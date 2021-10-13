#include <stdio.h>
#include <math.h>
#include <string.h>
#include "EllipseFit.h"

struct Pixel{int r, c;};

///-----------------------------------------------------------
/// Allocate a 2D array of noRows rows and noColumns columns
///
static double **AllocateMatrix(int noRows, int noColumns){
  double **m = new double *[noRows];


  for (int i=0; i<noRows; i++){
	  m[i] = new double[noColumns];
	  memset(m[i], 0, sizeof(double)*noColumns);
  } // end-for

  return m;
} //end=AllocateMatrix

///-----------------------------------------------------------
/// Deallocates a 2D array of noRows rows and noColumns columns
///
static void DeallocateMatrix(double **m, int noRows){
  for (int i=0; i<noRows; i++) delete m[i];
  delete m;
} //end-DeallocateMatrix

///-----------------------------------------------------------
/// This function performs the meat of the calculations for the
/// curve plotting.  Note that it is not a matrix multiplier in the
/// pure sense.  The first matrix is the curve matrix (each curve type
/// has its own matrix), and the second matrix is the geometry matrix
/// (defined by the control points).  The result is returned in the
/// third matrix.
///
static void multMatrix(double **m, double **g, double **mg){
	// First clear the return array
	for (int i = 0; i < 4; i++)
		for (int j = 0; j < 2; j++)
			mg[i][j] = 0;

	// Perform the matrix math
	for (int i = 0; i < 4; i++)
		for (int j = 0; j < 2; j++)
			for (int k = 0; k < 4; k++)
				mg[i][j] = mg[i][j] + (m[i][k] * g[k][j]);
} //end-multMatrix


///-----------------------------------------------------------
/// ROTATE
///
static void ROTATE(double **a, int i, int j, int k, int l, double tau, double s){
	double g, h;
	g = a[i][j]; h = a[k][l]; a[i][j] = g - s * (h + g * tau);
	a[k][l] = h + s * (g - h * tau);
}

///-----------------------------------------------------------
/// jacobi
///
static void jacobi(double **a, int n, double d[], double **v, int nrot){
	int j, iq, ip, i;
	double tresh, theta, tau, t, sm, s, h, g, c;

	double *b = new double[n+1];
	double *z = new double[n+1];
	memset(b, 0, sizeof(double)*(n+1));
	memset(z, 0, sizeof(double)*(n+1));

	for (ip = 1; ip <= n; ip++)
	{
		for (iq = 1; iq <= n; iq++) v[ip][iq] = 0.0;
		v[ip][ip] = 1.0;
	}
	for (ip = 1; ip <= n; ip++)
	{
		b[ip] = d[ip] = a[ip][ip];
		z[ip] = 0.0;
	}
	nrot = 0;
	for (i = 1; i <= 50; i++)
	{
		sm = 0.0;
		for (ip = 1; ip <= n - 1; ip++)
		{
			for (iq = ip + 1; iq <= n; iq++)
				sm += fabs(a[ip][iq]);
		}
		if (sm == 0.0)
		{
      delete b;
      delete z;
			return;
		}
		if (i < 4)
			tresh = 0.2 * sm / (n * n);
		else
			tresh = 0.0;
		for (ip = 1; ip <= n - 1; ip++)
		{
			for (iq = ip + 1; iq <= n; iq++)
			{
				g = 100.0 * fabs(a[ip][iq]);
//				if (i > 4 && fabs(d[ip]) + g == fabs(d[ip])
//				&& fabs(d[iq]) + g == fabs(d[iq]))

				if (i > 4 && g == 0.0)
					a[ip][iq] = 0.0;
				else if (fabs(a[ip][iq]) > tresh)
				{
					h = d[iq] - d[ip];
					if (g == 0.0)
						t = (a[ip][iq]) / h;
					else
					{
						theta = 0.5 * h / (a[ip][iq]);
						t = 1.0 / (fabs(theta) + sqrt(1.0 + theta * theta));
						if (theta < 0.0) t = -t;
					}
					c = 1.0 / sqrt(1 + t * t);
					s = t * c;
					tau = s / (1.0 + c);
					h = t * a[ip][iq];
					z[ip] -= h;
					z[iq] += h;
					d[ip] -= h;
					d[iq] += h;
					a[ip][iq] = 0.0;
									for (j = 1; j <= ip - 1; j++)
					{
						ROTATE(a, j, ip, j, iq, tau, s);
					}
					for (j = ip + 1; j <= iq - 1; j++)
					{
						ROTATE(a, ip, j, j, iq, tau, s);
					}
					for (j = iq + 1; j <= n; j++)
					{
						ROTATE(a, ip, j, iq, j, tau, s);
					}
					for (j = 1; j <= n; j++)
					{
						ROTATE(v, j, ip, j, iq, tau, s);
					}
					++nrot;
				}
			}
		}
		for (ip = 1; ip <= n; ip++)
		{
			b[ip] += z[ip];
			d[ip] = b[ip];
			z[ip] = 0.0;
		}
	}
	//printf("Too many iterations in routine JACOBI");
	delete b;
	delete z;
} //end-jacobi

///-----------------------------------------------------------
/// Perform the Cholesky decomposition
/// Return the lower triangular L  such that L*L'=A
///
static void choldc(double **a, int n, double **l){
	int i, j, k;
	double sum;
	double *p = new double[n+1];
	memset(p, 0, sizeof(double)*(n+1));


	for (i = 1; i <= n; i++)
	{
		for (j = i; j <= n; j++)
		{
			for (sum = a[i][j], k = i - 1; k >= 1; k--) sum -= a[i][k] * a[j][k];
			if (i == j)
			{
				if (sum <= 0.0)
				// printf("\nA is not poitive definite!");
				{ }
				else
					p[i] = sqrt(sum);
			}
			else
			{
				a[j][i] = sum / p[i];
			}
		}
	}
  for (i = 1; i <= n; i++){
    for (j = i; j <= n; j++){
			if (i == j)
				l[i][i] = p[i];
			else
			{
				l[j][i] = a[j][i];
				l[i][j] = 0.0;
			}
    } //end-for-inner
  } // end-for-outer

  delete p;
} //end-choldc

/********************************************************************/
/**    Calcola la inversa della matrice  B mettendo il risultato   **/
/**    in InvB . Il metodo usato per l'inversione e' quello di     **/
/**    Gauss-Jordan.   N e' l'ordine della matrice .               **/
/**    ritorna 0 se l'inversione  corretta altrimenti ritorna     **/
/**    SINGULAR .                                                  **/
/********************************************************************/
static int inverse(double **TB, double **InvB, int N){
	int k, i, j, p, q;
	double mult;
	double D, temp;
	double maxpivot;
	int npivot;
	double **B = AllocateMatrix(N+1, N+2);
	double **A = AllocateMatrix(N+1, 2*N+2);
	double **C = AllocateMatrix(N+1, N+1);
	double eps = 10e-20;

	for (k = 1; k <= N; k++)
		for (j = 1; j <= N; j++)
			B[k][j] = TB[k][j];

	for (k = 1; k <= N; k++)
	{
		for (j = 1; j <= N + 1; j++)
			A[k][j] = B[k][j];
		for (j = N + 2; j <= 2 * N + 1; j++)
			A[k][j] = (double)0;
		A[k][k - 1 + N + 2] = (double)1;
	}
	for (k = 1; k <= N; k++)
	{
		maxpivot = fabs((double)A[k][k]);
		npivot = k;
		for (i = k; i <= N; i++)
			if (maxpivot < fabs((double)A[i][k]))
			{
				maxpivot = fabs((double)A[i][k]);
				npivot = i;
			}
		if (maxpivot >= eps)
		{
			if (npivot != k)
				for (j = k; j <= 2 * N + 1; j++)
				{
					temp = A[npivot][j];
					A[npivot][j] = A[k][j];
					A[k][j] = temp;
				};
			D = A[k][k];
			for (j = 2 * N + 1; j >= k; j--)
				A[k][j] = A[k][j] / D;
			for (i = 1; i <= N; i++)
			{
				if (i != k)
				{
					mult = A[i][k];
					for (j = 2 * N + 1; j >= k; j--)
						A[i][j] = A[i][j] - mult * A[k][j];
				}
			}
		}
		else
		{  // printf("\n The matrix may be singular !!") ;

	    DeallocateMatrix(B, N+1);
	    DeallocateMatrix(A, N+1);
	    DeallocateMatrix(C, N+1);

			return (-1);
		} //end-else
	}
	/**   Copia il risultato nella matrice InvB  ***/
	for (k = 1, p = 1; k <= N; k++, p++)
		for (j = N + 2, q = 1; j <= 2 * N + 1; j++, q++)
			InvB[p][q] = A[k][j];

	DeallocateMatrix(B, N+1);
	DeallocateMatrix(A, N+1);
	DeallocateMatrix(C, N+1);

	return (0);
} //end-inverse


///-----------------------------------------------------------
/// AperB
///
static void AperB(double **_A, double **_B, double **_res,
                		   int _righA, int _colA, int _righB, int _colB){
	int p, q, l;
	for (p = 1; p <= _righA; p++)
		for (q = 1; q <= _colB; q++)
		{
			_res[p][q] = 0.0;
			for (l = 1; l <= _colA; l++)
				_res[p][q] = _res[p][q] + _A[p][l] * _B[l][q];
		}
} //end-AperB

///-----------------------------------------------------------
/// A_TperB
///
static void A_TperB(double **_A, double **_B, double **_res,
		                     int _righA, int _colA, int _righB, int _colB){
	int p, q, l;
	for (p = 1; p <= _colA; p++)
		for (q = 1; q <= _colB; q++)
		{
			_res[p][q] = 0.0;
			for (l = 1; l <= _righA; l++)
				_res[p][q] = _res[p][q] + _A[l][p] * _B[l][q];
		}
} //end-A_TperB

///-----------------------------------------------------------
/// AperB_T
///
static void AperB_T(double **_A, double **_B, double **_res,
		                     int _righA, int _colA, int _righB, int _colB){
	int p, q, l;
	for (p = 1; p <= _colA; p++)
		for (q = 1; q <= _colB; q++)
		{
			_res[p][q] = 0.0;
			for (l = 1; l <= _righA; l++)
				_res[p][q] = _res[p][q] + _A[p][l] * _B[q][l];
		}
} //end-AperB_T

///--------------------------------------------------------------------------------
/// Computes the ellipse equation that best fits the given set of points
/// Returns "true" on success, "false" on failure
///
bool EllipseFit(void *points, int noPoints, EllipseEquation *pResult, int mode){
  Pixel *pixels = (Pixel *)points;

	double **D = AllocateMatrix(noPoints+1, 7);
	double **S = AllocateMatrix(7, 7);
	double **Const = AllocateMatrix(7, 7);
	double **temp = AllocateMatrix(7, 7);
	double **L = AllocateMatrix(7, 7);
	double **C = AllocateMatrix(7, 7);

	double **invL = AllocateMatrix(7, 7);
	double *d = new double[7];
	double **V = AllocateMatrix(7, 7);
	double **sol = AllocateMatrix(7, 7);
	double tx, ty;
	int nrot = 0;

	memset(d, 0, sizeof(double)*7);


	switch (mode)
	{
		case (FPF):
			//fprintf(stderr, "EllipseFit: FPF mode");
			Const[1][3] = -2;
			Const[2][2] = 1;
			Const[3][1] = -2;
			break;
		case (BOOKSTEIN):
			//fprintf(stderr, "EllipseFit: BOOKSTEIN mode");
			Const[1][1] = 2;
			Const[2][2] = 1;
			Const[3][3] = 2;
	} //end-switch

	if (noPoints < 6)
		return false;

	// Now first fill design matrix
	for (int i = 1; i <= noPoints; i++){
		tx = (double)pixels[i-1].c;
		ty = (double)pixels[i-1].r;

		D[i][1] = tx * tx;
		D[i][2] = tx * ty;
		D[i][3] = ty * ty;
		D[i][4] = tx;
		D[i][5] = ty;
		D[i][6] = 1.0;
	} //end-for

	//pm(Const,"Constraint");
	// Now compute scatter matrix  S
	A_TperB(D, D, S, noPoints, 6, noPoints, 6);
	//pm(S,"Scatter");

	choldc(S, 6, L);
	//pm(L,"Cholesky");

	inverse(L, invL, 6);
	//pm(invL,"inverse");

	AperB_T(Const, invL, temp, 6, 6, 6, 6);
	AperB(invL, temp, C, 6, 6, 6, 6);
	//pm(C,"The C matrix");

	jacobi(C, 6, d, V, nrot);
	//pm(V,"The Eigenvectors");  /* OK */
	//pv(d,"The eigevalues");

	A_TperB(invL, V, sol, 6, 6, 6, 6);
	//pm(sol,"The GEV solution unnormalized");  /* SOl */

	// Now normalize them
	for (int j = 1; j <= 6; j++)  /* Scan columns */
	{
		double mod = 0.0;
		for (int i = 1; i <= 6; i++)
			mod += sol[i][j] * sol[i][j];
		for (int i = 1; i <= 6; i++)
			sol[i][j] /= sqrt(mod);
	}

	//pm(sol,"The GEV solution");  /* SOl */

	double zero = 10e-20;
	double minev = 10e+20;
	int solind = 0;
  int i;
	switch (mode)
	{
		case (BOOKSTEIN):  // smallest eigenvalue
			for (i = 1; i <= 6; i++)
				if (d[i] < minev && fabs(d[i]) > zero)
					solind = i;
			break;
		case (FPF):
			for (i = 1; i <= 6; i++)
				if (d[i] < 0 && fabs(d[i]) > zero)
					solind = i;
	}

#if 0
  fprintf(stderr, "SOLUTIONS: Selected: %d\n", solind);
  for (i=1; i<=6; i++){
    fprintf(stderr, "d[i]: %e, a: %.7lf, b: %.7lf, c: %.7lf, d: %.7lf, e: %.7lf, f: %.7lf\n",
            d[i], sol[1][i], sol[2][i], sol[3][i], sol[4][i], sol[5][i], sol[6][i]);
  } //end-for

#endif

  bool valid = true;
  if (solind == 0) valid = false;

	// Now fetch the right solution
	for (int j = 1; j <= 6; j++){
		pResult->coeff[j] = sol[j][solind];
	} //end-for

	DeallocateMatrix(D, noPoints+1);
	DeallocateMatrix(S, 7);
	DeallocateMatrix(Const, 7);
	DeallocateMatrix(temp, 7);
	DeallocateMatrix(L, 7);
	DeallocateMatrix(C, 7);
	DeallocateMatrix(invL, 7);
	delete d;
	DeallocateMatrix(V, 7);
	DeallocateMatrix(sol, 7);

  int len = (int)ComputeEllipsePerimeter(pResult->coeff);
  if (len<=0 || len > 50000) valid = false;

  return valid;
} // end-FitEllipse


///--------------------------------------------------------------------------------
/// Computes the ellipse equation that best fits the given set of points
/// Returns "true" on success, "false" on failure
///
bool EllipseFit(double *x, double *y, int noPoints, EllipseEquation *pResult, int mode){
	double **D = AllocateMatrix(noPoints+1, 7);
	double **S = AllocateMatrix(7, 7);
	double **Const = AllocateMatrix(7, 7);
	double **temp = AllocateMatrix(7, 7);
	double **L = AllocateMatrix(7, 7);
	double **C = AllocateMatrix(7, 7);

	double **invL = AllocateMatrix(7, 7);
	double *d = new double[7];
	double **V = AllocateMatrix(7, 7);
	double **sol = AllocateMatrix(7, 7);
	double tx, ty;
	int nrot = 0;

	memset(d, 0, sizeof(double)*7);


	switch (mode)
	{
		case (FPF):
			//fprintf(stderr, "EllipseFit: FPF mode");
			Const[1][3] = -2;
			Const[2][2] = 1;
			Const[3][1] = -2;
			break;
		case (BOOKSTEIN):
			//fprintf(stderr, "EllipseFit: BOOKSTEIN mode");
			Const[1][1] = 2;
			Const[2][2] = 1;
			Const[3][3] = 2;
	} //end-switch

	if (noPoints < 6)
		return false;

	// Now first fill design matrix
	for (int i = 1; i <= noPoints; i++){
		tx = x[i-1];
		ty = y[i-1];

		D[i][1] = tx * tx;
		D[i][2] = tx * ty;
		D[i][3] = ty * ty;
		D[i][4] = tx;
		D[i][5] = ty;
		D[i][6] = 1.0;
	} //end-for

	//pm(Const,"Constraint");
	// Now compute scatter matrix  S
	A_TperB(D, D, S, noPoints, 6, noPoints, 6);
	//pm(S,"Scatter");

	choldc(S, 6, L);
	//pm(L,"Cholesky");

	inverse(L, invL, 6);
	//pm(invL,"inverse");

	AperB_T(Const, invL, temp, 6, 6, 6, 6);
	AperB(invL, temp, C, 6, 6, 6, 6);
	//pm(C,"The C matrix");

	jacobi(C, 6, d, V, nrot);
	//pm(V,"The Eigenvectors");  /* OK */
	//pv(d,"The eigevalues");

	A_TperB(invL, V, sol, 6, 6, 6, 6);
	//pm(sol,"The GEV solution unnormalized");  /* SOl */

	// Now normalize them
	for (int j = 1; j <= 6; j++)  /* Scan columns */
	{
		double mod = 0.0;
		for (int i = 1; i <= 6; i++)
			mod += sol[i][j] * sol[i][j];
		for (int i = 1; i <= 6; i++)
			sol[i][j] /= sqrt(mod);
	}

	//pm(sol,"The GEV solution");  /* SOl */

	double zero = 10e-20;
	double minev = 10e+20;
	int solind = 0;
  int i;
	switch (mode)
	{
		case (BOOKSTEIN):  // smallest eigenvalue
			for (i = 1; i <= 6; i++)
				if (d[i] < minev && fabs(d[i]) > zero)
					solind = i;
			break;
		case (FPF):
			for (i = 1; i <= 6; i++)
				if (d[i] < 0 && fabs(d[i]) > zero)
					solind = i;
	}

#if 0
  fprintf(stderr, "SOLUTIONS: Selected: %d\n", solind);
  for (i=1; i<=6; i++){
    fprintf(stderr, "d[i]: %e, a: %.7lf, b: %.7lf, c: %.7lf, d: %.7lf, e: %.7lf, f: %.7lf\n",
            d[i], sol[1][i], sol[2][i], sol[3][i], sol[4][i], sol[5][i], sol[6][i]);
  } //end-for

#endif

  bool valid = true;
  if (solind == 0) valid = false;

  if (valid){
	  // Now fetch the right solution
	  for (int j = 1; j <= 6; j++){
		  pResult->coeff[j] = sol[j][solind];
	  } //end-for
  } //end-if

	DeallocateMatrix(D, noPoints+1);
	DeallocateMatrix(S, 7);
	DeallocateMatrix(Const, 7);
	DeallocateMatrix(temp, 7);
	DeallocateMatrix(L, 7);
	DeallocateMatrix(C, 7);
	DeallocateMatrix(invL, 7);
	delete d;
	DeallocateMatrix(V, 7);
	DeallocateMatrix(sol, 7);

  if (valid){
    int len = (int)ComputeEllipsePerimeter(pResult);
    if (len<=0 || len > 50000) valid = false;
  } //end-else

  return valid;
} // end-FitEllipse


/// ---------------------------------------------------------------------------
/// Given an ellipse equation, computes "noPoints" many consecutive points
/// on the ellipse periferi. These points can be used to draw the ellipse
/// noPoints must be an even number.
///
void ComputeEllipsePoints(double *pvec, double *px, double *py, int noPoints){
  if (noPoints%2) noPoints--;
	int npts = noPoints / 2;

	double **u = AllocateMatrix(3, npts+1);
	double **Aiu = AllocateMatrix(3, npts+1);
	double **L = AllocateMatrix(3, npts+1);
	double **B = AllocateMatrix(3, npts+1);
	double **Xpos = AllocateMatrix(3, npts+1);
	double **Xneg = AllocateMatrix(3, npts+1);
	double **ss1 = AllocateMatrix(3, npts+1);
	double **ss2 = AllocateMatrix(3, npts+1);
	double *lambda = new double[npts+1];
	double **uAiu = AllocateMatrix(3, npts+1);
	double **A = AllocateMatrix(3, 3);
	double **Ai = AllocateMatrix(3, 3);
	double **Aib = AllocateMatrix(3, 2);
	double **b = AllocateMatrix(3, 2);
	double **r1 = AllocateMatrix(2, 2);
	double Ao, Ax, Ay, Axx, Ayy, Axy;

	double pi = 3.14781;
	double theta;
	int i;
	int j;
	double kk;

	memset(lambda, 0, sizeof(double)*(npts+1));

	Ao = pvec[6];
	Ax = pvec[4];
	Ay = pvec[5];
	Axx = pvec[1];
	Ayy = pvec[3];
	Axy = pvec[2];

	A[1][1] = Axx; A[1][2] = Axy / 2;
	A[2][1] = Axy / 2; A[2][2] = Ayy;
	b[1][1] = Ax; b[2][1] = Ay;

	// Generate normals linspace
	for (i = 1, theta = 0.0; i <= npts; i++, theta += (pi / npts))
	{
		u[1][i] = cos(theta);
		u[2][i] = sin(theta);
	}

	inverse(A, Ai, 2);

	AperB(Ai, b, Aib, 2, 2, 2, 1);
	A_TperB(b, Aib, r1, 2, 1, 2, 1);
	r1[1][1] = r1[1][1] - 4 * Ao;

	AperB(Ai, u, Aiu, 2, 2, 2, npts);
	for (i = 1; i <= 2; i++)
		for (j = 1; j <= npts; j++)
			uAiu[i][j] = u[i][j] * Aiu[i][j];

	for (j = 1; j <= npts; j++)
	{
		if ((kk = (r1[1][1] / (uAiu[1][j] + uAiu[2][j]))) >= 0.0)
			lambda[j] = sqrt(kk);
		else
			lambda[j] = -1.0;
	}

	// Builds up B and L
	for (j = 1; j <= npts; j++)
		L[1][j] = L[2][j] = lambda[j];
	for (j = 1; j <= npts; j++)
	{
		B[1][j] = b[1][1];
		B[2][j] = b[2][1];
	}

	for (j = 1; j <= npts; j++)
	{
		ss1[1][j] = 0.5 * (L[1][j] * u[1][j] - B[1][j]);
		ss1[2][j] = 0.5 * (L[2][j] * u[2][j] - B[2][j]);
		ss2[1][j] = 0.5 * (-L[1][j] * u[1][j] - B[1][j]);
		ss2[2][j] = 0.5 * (-L[2][j] * u[2][j] - B[2][j]);
	}

	AperB(Ai, ss1, Xpos, 2, 2, 2, npts);
	AperB(Ai, ss2, Xneg, 2, 2, 2, npts);

	for (j = 1; j <= npts; j++)
	{
		if (lambda[j] == -1.0)
		{
			px[j-1] = -1;
			py[j-1] = -1;
			px[j-1 + npts] = -1;
			py[j-1 + npts] = -1;
		}
		else
		{
			px[j-1] = Xpos[1][j];
			py[j-1] = Xpos[2][j];
			px[j-1 + npts] = Xneg[1][j];
			py[j-1 + npts] = Xneg[2][j];
		}
	}

	DeallocateMatrix(u, 3);
	DeallocateMatrix(Aiu, 3);
	DeallocateMatrix(L, 3);
	DeallocateMatrix(B, 3);
	DeallocateMatrix(Xpos, 3);
	DeallocateMatrix(Xneg, 3);
	DeallocateMatrix(ss1, 3);
	DeallocateMatrix(ss2, 3);
	delete lambda;
	DeallocateMatrix(uAiu, 3);
	DeallocateMatrix(A, 3);
	DeallocateMatrix(Ai, 3);
	DeallocateMatrix(Aib, 3);
	DeallocateMatrix(b, 3);
	DeallocateMatrix(r1, 2);
} //end-ComputeEllipsePoints

/// ---------------------------------------------------------------------------
/// Given an ellipse equation, numerically computes the length of the perimeter of the ellipse
///
double ComputeEllipsePerimeter(double *pvec){
  /// Compute 8 points on the ellipse's perimeter
  /// Then find the center, x-axis length and y-axis length
  /// Finally, compute the length of the circumference of the ellipse
  /// and return that many points
#define NO_POINTS 16
  double x[NO_POINTS];
  double y[NO_POINTS];

  ComputeEllipsePoints(pvec, x, y, NO_POINTS);

  double xc = 0;
  double yc = 0;
  double xmin = 1e10;
  double xmax = 0;
  double ymin = 1e10;
  double ymax = 0;
  for (int n=0; n<NO_POINTS; n++){
    xc += x[n];
    yc += y[n];

    if      (x[n] < xmin) xmin = x[n];
    else if (x[n] > xmax) xmax = x[n];

    if      (y[n] < ymin) ymin = y[n];
    else if (y[n] > ymax) ymax = y[n];
  } //end-for

  xc /= NO_POINTS;
  yc /= NO_POINTS;
  double a = (xc-xmin + xmax-xc)/2;
  double b = (yc-ymin + ymax-yc)/2;
//  return 3.14159*sqrt(2*a*a + 2*b*b)*1.25;
  return 3.14159*sqrt(2*a*a + 2*b*b);
#undef NO_POINTS
} //end-ComputeEllipsePeripmeter


/// ---------------------------------------------------------------------------
/// Given an ellipse equation, computes the length of the perimeter of the ellipse
/// Calculates the ellipse perimeter wrt the Ramajunan II formula
///
double ComputeEllipsePerimeter(EllipseEquation *eq){
//	double A = coefs[0], B = coefs[1], C = coefs[2], D = coefs[3], E = coefs[4], F = coefs[5];

  double mult = 1;

	double A = eq->A()*mult;
  double B = eq->B()*mult;
  double C = eq->C()*mult;
  double D = eq->D()*mult;
  double E = eq->E()*mult;
  double F = eq->F()*mult;

	double A2, B2, C2, D2, E2, F2, theta;  //rotated coefficients
	double			   D3, E3, F3;             //ellipse form coefficients
	double cX, cY, a, b;		               //(cX,cY) center, a & b: semimajor & semiminor axes
	double h;							                 //h = (a-b)^2 / (a+b)^2
	bool rotation = false;

#define pi 3.14159265

	//Normalize coefficients
	B /= A; C /= A; D /= A; E /= A; F /= A; A /= A;

	if (B == 0) //Then not need to rotate the axes
	{
		A2 = A;
		B2 = B;
		C2 = C;
		D2 = D;
		E2 = E;
		F2 = F;

//		printf("No Rotation is applied\n\n");
	}
	else if (B != 0) //Rotate the axes
	{
		rotation = true;

		//Determine the rotation angle (in radians)
		theta = atan(B / (A - C)) / 2;

		//Compute the coefficients wrt the new coordinate system
		A2 = 0.5 * (A * (1 + cos(2 * theta) + B * sin(2 * theta) + C * (1 - cos ( 2 * theta))));

		B2 = (C - A) * sin(2 * theta) + B * cos(2 * theta); //B2 should turn to be zero?

		C2 = 0.5 * (A * (1 - cos(2 * theta) - B * sin(2 * theta) + C * (1 + cos ( 2 * theta))));

		D2 = D * cos(theta) + E * sin(theta);

		E2 = -D * sin(theta) + E * cos(theta);

		F2 = F;

//		printf("Rotation in degrees: %.2f\n\n", theta * 180 / pi);
	}

	//Transform the conic equation into the ellipse form
	D3 = D2 / A2; //normalize x term's coef
	//A3 = 1;     //A2 / A2

	E3 = E2 / C2; //normalize y term's coef
	//C3 = 1;     //C2 / C2

	cX = -(D3 / 2);	//center X
	cY = -(E3 / 2);	//center Y

	F3 = A2 * pow(cX, 2.0) + C2 * pow(cY, 2.0) - F2;

	//semimajor axis
	a = sqrt(F3 / A2);
	//semiminor axis
	b = sqrt(F3 / C2);

	//Center coordinates have to be re-transformed if rotation is applied!
	if(rotation)
	{
		double tmpX = cX, tmpY = cY;
		cX = tmpX * cos(theta) - tmpY * sin(theta);
		cY = tmpX * sin(theta) + tmpY * cos(theta);
	}

	///printf("F(x,y)= (x + %.2f) / %.3f + (y + %.2f) / %.3f = 1\n\n", -cX, a, -cY, b);
//	printf("Semi-major axis: a = %.3f\nSemi-minor axis: b = %.3f\n", a, b);
//	printf("Center of ellipse: (%.2f , %.2f)\n\n", cX, cY);

	//Perimeter Computation(s)
	h = pow((a - b), 2.0) / pow((a + b), 2.0);

	//Ramajunan I
	//double P1 = pi * (a + b) * (3 - sqrt(4 - h));
	///printf("Perimeter of the ellipse is %.5f (Ramajunan I)\n", P1);

	//Ramajunan II
	double P2 = pi * (a + b) * (1 + 3*h / (10 + sqrt(4 - 3*h)));
//	printf("Perimeter of the ellipse is %.5f (Ramajunan II)\n", P2);

	//lise kitabindaki formul
//	double P3 = 2 * pi * sqrt(0.5 * (a*a + b*b));
//	printf("Perimeter of the ellipse is %.5f (simple formula)\n", P3);

	return P2;
#undef pi
} //end-ComputeEllipsePerimeter


/// ---------------------------------------------------------------------------
/// Computes the center, major and minor axis lengths of an ellipse
///
void ComputeEllipseCenterAndAxisLengths(EllipseEquation *eq, double *pxc, double *pyc, double *pmajorAxisLength, double *pminorAxisLength){
//	double A = coefs[0], B = coefs[1], C = coefs[2], D = coefs[3], E = coefs[4], F = coefs[5];

  double mult = 1;

  double A = eq->A()*mult;
  double B = eq->B()*mult;
  double C = eq->C()*mult;
  double D = eq->D()*mult;
  double E = eq->E()*mult;
  double F = eq->F()*mult;

	double A2, B2, C2, D2, E2, F2, theta;  //rotated coefficients
	double			   D3, E3, F3;             //ellipse form coefficients
	double cX, cY, a, b;		               //(cX,cY) center, a & b: semimajor & semiminor axes
	bool rotation = false;

#define pi 3.14159265

	//Normalize coefficients
	B /= A; C /= A; D /= A; E /= A; F /= A; A /= A;

	if (B == 0) //Then not need to rotate the axes
	{
		A2 = A;
		B2 = B;
		C2 = C;
		D2 = D;
		E2 = E;
		F2 = F;

//		printf("No Rotation is applied\n\n");
	}
	else if (B != 0) //Rotate the axes
	{
		rotation = true;

		//Determine the rotation angle (in radians)
		theta = atan(B / (A - C)) / 2;

		//Compute the coefficients wrt the new coordinate system
		A2 = 0.5 * (A * (1 + cos(2 * theta) + B * sin(2 * theta) + C * (1 - cos ( 2 * theta))));

		B2 = (C - A) * sin(2 * theta) + B * cos(2 * theta); //B2 should turn to be zero?

		C2 = 0.5 * (A * (1 - cos(2 * theta) - B * sin(2 * theta) + C * (1 + cos ( 2 * theta))));

		D2 = D * cos(theta) + E * sin(theta);

		E2 = -D * sin(theta) + E * cos(theta);

		F2 = F;

//		printf("Rotation in degrees: %.2f\n\n", theta * 180 / pi);
	}

	//Transform the conic equation into the ellipse form
	D3 = D2 / A2; //normalize x term's coef
	//A3 = 1;     //A2 / A2

	E3 = E2 / C2; //normalize y term's coef
	//C3 = 1;     //C2 / C2

	cX = -(D3 / 2);	//center X
	cY = -(E3 / 2);	//center Y

	F3 = A2 * pow(cX, 2.0) + C2 * pow(cY, 2.0) - F2;

	//semimajor axis
	a = sqrt(F3 / A2);
	//semiminor axis
	b = sqrt(F3 / C2);

	//Center coordinates have to be re-transformed if rotation is applied!
	if(rotation)
	{
		double tmpX = cX, tmpY = cY;
		cX = tmpX * cos(theta) - tmpY * sin(theta);
		cY = tmpX * sin(theta) + tmpY * cos(theta);
	}

  *pxc = cX;
  *pyc = cY;

  if (a > b){
    *pmajorAxisLength = a;
    *pminorAxisLength = b;
  } else {
    *pmajorAxisLength = b;
    *pminorAxisLength = a;
  } //end-else

#undef pi
} //end-ComputeEllipseCenterAndAxisLengths

/// ---------------------------------------------------------------------------
/// Given an ellipse equation that was fitted to "noPoints" many points
/// stored in "px" & "py", computes the least square error
/// Slow but accurate method
///
double ComputeEllipseErrorSlow(EllipseEquation *eq, double *px, double *py, int noPoints){
  int len = (int)ComputeEllipsePerimeter(eq->coeff)*2;
  if (len%2) len--;

  double *x = new double[len];
  double *y = new double[len];

  /// Compute points on the periphery of the ellipse
  ComputeEllipsePoints(eq->coeff, x, y, len);

  // First, find the ellipse points that is closest to the first user point
  double error = 0;

  for (int i=0; i<noPoints; i++){
    double min = 1e10;

    int index = -1;
    for (int j=0; j<len; j++){
      double dx = px[i] - x[j];
      double dy = py[i] - y[j];
      double d = dx*dx + dy*dy;   // square of the distance
      if (d < min){ min = d; index=j;}
    } //end-for

    error += min;
  } //end-for

  error = sqrt(error/noPoints);

  delete x;
  delete y;

  return error;
} //end-ComputeEllipseErrorSlow

/// ---------------------------------------------------------------------------
/// Given an ellipse equation that was fitted to "noPoints" many points
/// stored in "px" & "py", computes the least square error.
/// Fast & accurate method
///
void ComputeDistances2Ellipse(EllipseEquation *eq, double *px, double *py, int noPoints, double *distances){
  double error = 0;

  double A = eq->A();
  double B = eq->B();
  double C = eq->C();
  double D = eq->D();
  double E = eq->E();
  double F = eq->F();

  double xc, yc, major, minor;
  ComputeEllipseCenterAndAxisLengths(eq, &xc, &yc, &major, &minor);

  for (int i=0; i<noPoints; i++){
    double dx = px[i]-xc;
    double dy = py[i]-yc;

    double min;
    double xs, ys;

    if (fabs(dx) > fabs(dy)){
      // The line equation is of the form: y = mx+n
      double m = dy/dx;
      double n = yc-m*xc;

      // a*x^2 + b*x + c
      double a = A + B*m + C*m*m;
      double b = B*n +2*C*m*n + D + E*m;
      double c = C*n*n + E*n + F;
      double det = b*b - 4*a*c;
      if (det<0) det = 0;
      double x1 = -(b+sqrt(det))/(2*a);
      double x2 = -(b-sqrt(det))/(2*a);

      double y1 = m*x1+n;
      double y2 = m*x2+n;

      dx = px[i]-x1;
      dy = py[i]-y1;
      double d1 = dx*dx + dy*dy;

      dx = px[i]-x2;
      dy = py[i]-y2;
      double d2 = dx*dx + dy*dy;

      if (d1<d2){min=d1; xs=x1; ys=y1;}
      else      {min=d2; xs=x2; ys=y2;}

    } else {
      // The line equation is of the form: x = my+n
      double m = dx/dy;
      double n = xc-m*yc;

      // a*y^2 + b*y + c
      double a = A*m*m + B*m + C;
      double b = 2*A*m*n + B*n + D*m + E;
      double c = A*n*n +D*n + F;
      double det = b*b - 4*a*c;
      if (det<0) det = 0;
      double y1 = -(b+sqrt(det))/(2*a);
      double y2 = -(b-sqrt(det))/(2*a);

      double x1 = m*y1+n;
      double x2 = m*y2+n;

      dx = px[i]-x1;
      dy = py[i]-y1;
      double d1 = dx*dx + dy*dy;

      dx = px[i]-x2;
      dy = py[i]-y2;
      double d2 = dx*dx + dy*dy;

      if (d1<d2){min=d1; xs=x1; ys=y1;}
      else      {min=d2; xs=x2; ys=y2;}
    } //end-else

    // Refine the search in the vicinity of (xs, ys)
	double delta = 0.1;
	double x = xs;
    while (1){
      x += delta;

      double a = C;
      double b = B*x + E;
      double c = A*x*x + D*x + F;
      double det = b*b - 4*a*c;
      if (det < 0) det = 0;

      double y1 = -(b+sqrt(det))/(2*a);
      double y2 = -(b-sqrt(det))/(2*a);

      dx = px[i]-x;
      dy = py[i]-y1;
      double d1 = dx*dx + dy*dy;

      dy = py[i]-y2;
      double d2 = dx*dx + dy*dy;

      if      (d1 <= min){min = d1;}
      else if (d2 <= min){min = d2;}
      else break;
    } //end-while

    x = xs;
    while (1){
      x -= delta;

      double a = C;
      double b = B*x + E;
      double c = A*x*x + D*x + F;
      double det = b*b - 4*a*c;
      if (det < 0) det = 0;

      double y1 = -(b+sqrt(det))/(2*a);
      double y2 = -(b-sqrt(det))/(2*a);

      dx = px[i]-x;
      dy = py[i]-y1;
      double d1 = dx*dx + dy*dy;

      dy = py[i]-y2;
      double d2 = dx*dx + dy*dy;

      if      (d1 <= min){min = d1;}
      else if (d2 <= min){min = d2;}
      else break;
    } //end-while

	distances[i] = min;
//	if (distances[i] > 4.0) distances[i] = 4.0;
  } //end-for
} //end-ComputeDistances2Ellipse

#if 0
/// -----------------------------------------------------------------------------------
/// 
double MyEllipseFitError(EllipseEquation* eq, double* px, double* py, int noPoints) {
	double error = 0;

	double A = eq->A();
	double B = eq->B();
	double C = eq->C();
	double D = eq->D();
	double E = eq->E();
	double F = eq->F();

	double xc, yc, major, minor;
	ComputeEllipseCenterAndAxisLengths(eq, &xc, &yc, &major, &minor);

	for (int i = 0; i < noPoints; i++) {
		double dx = px[i] - xc;
		double dy = py[i] - yc;

		double min;
		double xs, ys;

		if (fabs(dx) > fabs(dy)) {
			// The line equation is of the form: y = mx+n
			double m = dy / dx;
			double n = yc - m * xc;

			// a*x^2 + b*x + c
			double a = A + B * m + C * m * m;
			double b = B * n + 2 * C * m * n + D + E * m;
			double c = C * n * n + E * n + F;
			double det = b * b - 4 * a * c;
			if (det < 0) det = 0;
			double x1 = -(b + sqrt(det)) / (2 * a);
			double x2 = -(b - sqrt(det)) / (2 * a);

			double y1 = m * x1 + n;
			double y2 = m * x2 + n;

			dx = px[i] - x1;
			dy = py[i] - y1;
			double d1 = dx * dx + dy * dy;

			dx = px[i] - x2;
			dy = py[i] - y2;
			double d2 = dx * dx + dy * dy;

			if (d1 < d2) { min = d1; xs = x1; ys = y1; }
			else { min = d2; xs = x2; ys = y2; }

		}
		else {
			// The line equation is of the form: x = my+n
			double m = dx / dy;
			double n = xc - m * yc;

			// a*y^2 + b*y + c
			double a = A * m * m + B * m + C;
			double b = 2 * A * m * n + B * n + D * m + E;
			double c = A * n * n + D * n + F;
			double det = b * b - 4 * a * c;
			if (det < 0) det = 0;
			double y1 = -(b + sqrt(det)) / (2 * a);
			double y2 = -(b - sqrt(det)) / (2 * a);

			double x1 = m * y1 + n;
			double x2 = m * y2 + n;

			dx = px[i] - x1;
			dy = py[i] - y1;
			double d1 = dx * dx + dy * dy;

			dx = px[i] - x2;
			dy = py[i] - y2;
			double d2 = dx * dx + dy * dy;

			if (d1 < d2) { min = d1; xs = x1; ys = y1; }
			else { min = d2; xs = x2; ys = y2; }
		} //end-else

		error += min;
	} //end-for

	error = sqrt(error / noPoints);

	return error;
} //end-MyEllipseFitError
#endif