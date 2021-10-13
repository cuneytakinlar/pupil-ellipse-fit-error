#ifndef ELLIPSE_H
#define ELLIPSE_H 

#include "EDDefines.h"

//conflicts Pixel in the EDDefines
//typedef struct pixel{	int x;	int y; }Pixel;


//This Ellipse class is designed to work with the "Screen Coordinate System"
//which has the origin at the left-top corner.
//This is why the word "Pixel" is especially preferred instead of "Point".
//Thus, all points' y component are reflected with a convenient transform.

namespace Conic {
	class Ellipse
	{

	private:
		//Ellipse coefficients
		double A1, B1, C1, D1, E1, F1;	//ellipse coefficients
		double A2, B2, C2, D2, E2, F2;
		double A3, B3, C3, D3, E3, F3;
		double cX, cY;				//center of ellipse
		double a, b, eccentricity;	//semimajor and semiminor axes, and eccentricity (a/b)
		double perimeter, alpha;	//perimeter and spanned angle (the ratio)			
		double rotation;			//rotation angle of the ellipse in radians
		double fitError;			//average fit error (per point)
		double rmsError;			//root mean square fitting error
		double normRmsErr;			//root mean square fitting error

		//variables that are being used multiple times
		double a2_b2; //(a^2 - b^2) is frequently used in point distance

		//points...
		int noPoints;
		double* xPoints;	//x coordinate of points in double
		double* yPoints;	//y coordinate of points in double
		double* estRadians;	//angle (in radians) estimations for 

		Pixel* fitPoints;	//points that ellipse will be fit
		Pixel* drawPoints;	//points that are the polygonal approximation of the ellipse
		Pixel* closePoints;	//closest points to the fit points on the ellipse curve (with the same order, respectively)

		//initialize the parameters for safety
		void InitParams();

		//One step of NR method: f / fPrime
		inline double NewtonRaphsonIteration(double theta, double x, double y);

		//variable steps NR estimation for theta of the closest ellipse point
		//One iteration of NR method: f - fPrime
		double NewtonRaphsonThetaEstimation(double theta, double x, double y);

		//returns distance between two points
		inline double DistanceBwPoints(Pixel p1, Pixel p2);

		inline double DistanceBwPoints(double x1, double y1, double x2, double y2);

		//returns the SQUARED distance between two points
		inline double SquaredDistanceBwPoints(Pixel p1, Pixel p2);

		inline double SquaredDistanceBwPoints(double x1, double y1, double x2, double y2);

		//Computes and returns of a pixels' distance to the ellipse
		double GetDistance(Pixel point, double& tmpEst);

		//overload for double version
		double GetDistance(double pX, double pY, double& estimation);

		//Computes and returns of a pixels' squared distance to the ellipse
		double GetSquaredDistance(Pixel point, double& tmpEst);

		//overload for double version
		double GetSquaredDistance(double pX, double pY, double& estimation);

	public:
		Ellipse(double coefs[6]);

		Ellipse(Pixel* points, int numPoints);

		Ellipse(double* pX, double* pY, int numPoints);

		~Ellipse();

		int GetNoPoints();

		Pixel GetCenter();

		double GetCenterX();

		double GetCenterY();

		double GetRotation();

		double GetSemiMajorAxis();

		double GetSemiMinorAxis();

		double GetEccentricity();

		//Gets the coefficients of the ellipse in conic form
		void GetCoefficients(double*);

		//Draws an ellipse with the desired resolution i.t.o. points
		//resolution may be selected wrt the perimeter
		Pixel* DrawEllipse(int resolution);

		double GetPerimeter();

		double GetSpannedAngleRatio();

		//Computes and returns the fitting error of the points that
		//are used to generate the ellipse equation
		double GetRmsFittingError();

		double GetRmsFittingError(Pixel* pixels, int noPixels);

		//Rms fitting error per point
		//(GetRmsFittingError / noPoints) x 100
		double GetNormalizedRmsFitErr();

		//Computes and returns the fitting error of the points that
		//are used to generate the ellipse equation
		double GetAverageFittingError();

		//Computes and returns the closest points to the fitting points 
		//GetFittingError has to be invoked before calling this fcn
		Pixel* GetClosestPoints();

		//Gets the closest point to the test point and returns the distance
		double GetClosestPointAndDistance(Pixel test, Pixel& closest);

		//overload for the double version
		double GetClosestPointAndDistance(double testX, double testY, Pixel& closest);
	};
}
#endif