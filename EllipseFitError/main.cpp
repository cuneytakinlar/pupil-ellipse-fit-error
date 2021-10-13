#pragma warning(disable:4996)  // For _CRT_SECURE_NO_WARNINGS

#include <iostream>
#include <stdio.h>
#include <opencv.hpp>
#include <cxcore.hpp>
#include <imgcodecs.hpp>
#include <highgui.hpp>

#include "EDDefines.h"
#include "EdgeDrawing.h"
#include "ExtSegmentMap.h"
#include "Ellipse.h"
#include "EllipseFit.h"

using namespace cv;
using namespace std;
using namespace Conic;

#ifdef _DEBUG
#pragma comment(lib, "C:/opencv-3.4.0/build/x64/vc15/lib/opencv_world340d.lib")
#else
#pragma comment(lib, "C:/opencv-3.4.0/build/x64/vc15/lib/opencv_world340.lib")
#endif

/// Compute the distance between ALL detected blob pixels & the GT ellipse
extern "C" __declspec(dllexport) double EllipseFitError(int, int, uint8_t output[], uint8_t gt[], double distances[]);

/// --------------------------------------------------------------------------------
/// Compute the distance between ALL detected blob pixels & the GT ellipse
/// 
double EllipseFitError(int width, int height, uint8_t output[], uint8_t gt[], double distances[]) {
	cv::Mat input(height, width, CV_8UC1);

	// Find the boundary of the GT pupil		
	unsigned char* ptr = input.ptr(0, 0);
	for (int i = 0; i < height; i++) {
		for (int j = 0; j < width; j++) {
			ptr[i * width + j] = gt[i * width + j];
		} //end-for
	} //end-for

	EDImage* ed1 = new EDImage(input.cols, input.rows);
	ed1->image = input.data;
	DetectEdgesByEdgeDrawing(ed1);
	int noGtSegments = ed1->map->noSegments;

	// Find the boundary of the pupil output by the model
	for (int i = 0; i < height; i++) {
		for (int j = 0; j < width; j++) {
			ptr[i * width + j] = output[i * width + j];
		} //end-for
	} //end-for

	EDImage* ed2 = new EDImage(input.cols, input.rows);
	ed2->image = input.data;
	DetectEdgesByEdgeDrawing(ed2);
	int noOutputSegments = ed2->map->noSegments;

	// If both have no pupil, error is 0
	if (noGtSegments == 0 && noOutputSegments == 0) {
		delete ed1;
		delete ed2;
		return 0;
	} // end-if

	// If there is no pupil in GT, but there is a pupil in output, return max error
	// Similarly, if there is a pupil in GT, but there is no pupil in output, return max error
	if (noGtSegments == 0 || noOutputSegments == 0) {
		delete ed1;
		delete ed2;
		distances[0] = 1.0;
		return 1;
	} // end-if

	// Take the pixels in GT pupil and fit an ellipse to these pixels
	int noPixels = ed1->map->segments[0].noPixels;
	double* px = new double[noPixels];
	double* py = new double[noPixels];
	for (int i = 0; i < noPixels; i++) {
		px[i] = ed1->map->segments[0].pixels[i].c;
		py[i] = ed1->map->segments[0].pixels[i].r;
	} // end-for

	EllipseEquation eq;
	bool ret = EllipseFit(px, py, noPixels, &eq);

	delete[] px;
	delete[] py;

	// If ellipse fit fails: This should not happen
	if (ret == false) {
		delete ed1;
		delete ed2;
		distances[0] = 1.0;
		return 1;
	} //end=if

	// Now, take the longest segment to be the pupil segment in the output
	int longestSegID = 0;
	for (int i = 1; i < noOutputSegments; i++) {
		if (ed2->map->segments[i].noPixels > ed2->map->segments[longestSegID].noPixels) {
			longestSegID = i;
		} //end-if
	} //end-for

	// Take the pixels that make up the detected pupil boundary
	noPixels = ed2->map->segments[longestSegID].noPixels;
	px = new double[noPixels];
	py = new double[noPixels];
	for (int i = 0; i < noPixels; i++) {
		px[i] = ed2->map->segments[longestSegID].pixels[i].c;
		py[i] = ed2->map->segments[longestSegID].pixels[i].r;
	} // end-for

	// Finally, compute the average distance between the detected pupil boundary pixels & the GT ellipse
	ComputeDistances2Ellipse(&eq, px, py, noPixels, distances);

	delete[] px;
	delete[] py;

	delete ed1;
	delete ed2;

	return noPixels;
} //end-EllipseFitError

#if 0
int main() {
	for (int i = 1; i < 50; i++) {
		char filename[100];
		sprintf(filename, "D:/Home/DLWS/PythonDll&Usage/CustomLossDLL/mytestimages/out/%05d.png", i);
		Mat eyeImage = imread(filename, IMREAD_GRAYSCALE);

		sprintf(filename, "D:/Home/DLWS/PythonDll&Usage/CustomLossDLL/mytestimages/gt/%05d.png", i);
		Mat gtImage = imread(filename, IMREAD_GRAYSCALE);

		double* distances = new double[eyeImage.rows * eyeImage.cols];
		int noPixels = EllipseFitError(eyeImage.cols, eyeImage.rows, eyeImage.data, gtImage.data, distances);

		double error = 0;
		for (int k = 0; k < noPixels; k++) error += distances[k];
		error = sqrt(error/ noPixels);
		delete[] distances;

		printf("i: %2d, noPixels: %4d, Error: %.5lf\n", i, noPixels, error);
	} //end-for

	return 0;
} //end-main
#endif