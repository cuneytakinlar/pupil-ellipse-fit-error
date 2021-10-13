#ifndef EXTSEGMENTMAP_H
#define EXTSEGMENTMAP_H

#include <stdlib.h>

#include "Ellipse.h"

//enum SegState {SHORT = 0, CANDIDATE, CURVY, NEAR_CIRCULAR};

class ExtSegment
{
public:
	int noPixels;
	double distance;	//aperture between first and last pixels
	double gradEntropy;	
	//SegState segState;

	double *gradients;  //former int*
	Pixel *pixels;		//pointer to the pixels in the ED
	Conic::Ellipse *ellipse;	//used if only the segment is complete pupil 

	int noCorners;
	int *corners;		//index array for corners

public:
	ExtSegment()
	{
		noPixels = 0;
		noCorners = 0;
		distance = 0;	//aperture between first and last pixels
		gradEntropy = 0;

		ellipse = NULL;
		gradients = NULL;
		corners = NULL;
	}

	~ExtSegment()
	{
		delete ellipse;
		delete[] gradients;
		delete[] corners;
	}
};

class ExtSegmentMap
{
public:
	int noSegments;

	ExtSegment* pSegments = NULL;
	//Pixel *pixels;	//uses the same memory within the ED

public:
	ExtSegmentMap(int noSeg)
	{
		noSegments = noSeg;

		pSegments = new ExtSegment[noSegments];
	}

	~ExtSegmentMap()
	{
		if (pSegments != NULL) delete pSegments;
	}
};


#endif