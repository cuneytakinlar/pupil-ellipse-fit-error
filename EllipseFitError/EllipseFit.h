#ifndef _ELLIPSE_FIT_H_
#define _ELLIPSE_FIT_H_

///----------------------------------------------------------
/// Ellipse Equation is of the form:
/// Ax^2 + Bxy + Cy^2 + Dx + Ey + F = 0
///
struct EllipseEquation {
  double coeff[7];  // coeff[1] = A

  EllipseEquation(){
    for (int i=0; i<7; i++) coeff[i] = 0;
  } //end-EllipseEquation

  double A(){return coeff[1];}
  double B(){return coeff[2];}
  double C(){return coeff[3];}
  double D(){return coeff[4];}
  double E(){return coeff[5];}
  double F(){return coeff[6];}
};


//#ifdef __cplusplus
//extern "C" {
//#endif /* __cplusplus */


#define BOOKSTEIN 0       // method1
#define FPF       1       // method2

/// Computes the ellipse equation that best fits the given set of points
/// Returns "true" on success, "false" on failure
bool EllipseFit(void *points, int noPoints, EllipseEquation *pResult, int mode=FPF);

/// Computes the ellipse equation that best fits the given set of points
/// Returns "true" on success, "false" on failure
bool EllipseFit(double *x, double *y, int noPoints, EllipseEquation *pResult, int mode=FPF);

/// Given an ellipse equation, computes the length of the perimeter of the ellipse
/// This is my numeric solution, which gives incorrect results, but is OK for the time being
double ComputeEllipsePerimeter(double *pvec);

/// Given an ellipse equation, computes the length of the perimeter of the ellipse
double ComputeEllipsePerimeter(EllipseEquation *eq);

/// Computes the center, major and minor axis lengths of an ellipse
void ComputeEllipseCenterAndAxisLengths(EllipseEquation *eq, double *pxc, double *pyc, double *pmajorAxisLength, double *pminorAxisLength);

/// Given an ellipse equation, computes "noPoints" many consecutive points on the ellipse periferi. 
/// noPoints must be an even number.
void ComputeEllipsePoints(double *pvec, double *px, double *py, int noPoints);

/// Given an ellipse equation that was fitted to "noPoints" many points
/// stored in "px" & "py", computes the least square error
/// Slow but accurate method
double ComputeEllipseErrorSlow(EllipseEquation *eq, double *px, double *py, int noPoints);

/// Given an ellipse equation that was fitted to "noPoints" many points
/// stored in "px" & "py", computes the distance of each point to the nearest ellipse point
/// Fast & accurate method
void ComputeDistances2Ellipse(EllipseEquation *eq, double *px, double *py, int noPoints, double *distances);

/// -------------------------------------------------------------------------------
/// My Ellipse fit error used as the regularization term for UNet
/// 
double MyEllipseFitError(EllipseEquation* eq, double* px, double* py, int noPoints);

//#ifdef __cplusplus
//}
//#endif /* __cplusplus */


#endif