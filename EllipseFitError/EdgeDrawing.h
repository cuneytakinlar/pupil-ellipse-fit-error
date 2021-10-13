#ifndef _EDGE_DRAWING_H_
#define _EDGE_DRAWING_H_

#include <stdlib.h>
#include "EDDefines.h"

///-----------------------------------------------------------------------
/// EDImage class declaration below
///
struct EDImage {
public:
  int width, height;          // Width and height of the image

  unsigned char *image;       // Image being manipulated
  unsigned char *lowPassImg;  // Image obtained by Gauss filter
  unsigned char *edgeImg;     // Edge image
  unsigned char *dirImg;      // Direction image
  short *gradImg;             // Gradient image
  short *gXImage;             // Vertical gradient image
  short *gYImage;             // Horizontal gradient image

  // parameters to be used during edge detection
  GradientOperator GRADIENT_OPERATOR;    // Operator to be used for gradient computation -- SOBEL by default
  int GRADIENT_THRESH;        // Gradient threshold - Pixels with gradient values lower than this are ignored
  int ANCHOR_THRESH;          // For a pixel to be an anchor, its gradient must be greater than both of its neighbors by this amount
  int SCAN_INTERVAL;          // Used in the determination of the anchor points. Specifies how many rows/columns we need to skip
                              // during a search for anchor points - default value is 1
  int MIN_PATH_LEN;           // An edge segment is accepted only if it has at least MIN_PATH_LEN pixels

  EdgeMap *map;               // Edge map obtained during edge detection
 
  // Temporary variables used by JoinAnchorPoints
  Pixel *tmpPixels;           // Temporary pixels array
  void *tmpChains;            // Temporary chains array
  void *tmpStack;             // Temporary stack array

  Pixel *pixels;              // Pixels array (all edge pixels are stored here)

  // Used during anchor computation & sorting
  bool useSortedAnchors;      // Sort the anchors & start with the one having the max. gradient value during linking. False by default.
  int maxGradValue;           // Value of the maximum gradient
  int *C;                     // We use counting sort. C is the count array
  int *A;                     // A hold the offsets of the anchors in sorted order 
  int noAnchors;              // Total # of anchors
  
  bool validateSegments;      // Should be validate edge segments? False by defaul

  // Statistic information
  int noImages;               // Total number of images processed
  double gaussTime;           // Time spent on gauss filter
  double gradientTime;        // Time spent on gradient computation
  double anchorTime;          // Time spent on anchor point computation
  double joinTime;            // Joining anchors time
  double validationTime;      // Time spent on segment validation

public:
  /// Constructor
  EDImage(int _width, int _height);

  /// Destructor
  ~EDImage();
};

/// Detect Edges by Edge Drawing
void DetectEdgesByEdgeDrawing(EDImage *ed);

/// Detect Edges by Edge Drawing & validate them by DMM algorithm using the ParameterFree ED
void DetectEdgesByEDParameterFree(EDImage *ed);

/// Splits segments into multiple segments if they contain closed objects in them
/// Returns the new number of segments in the edgemap
int SplitSegments(EdgeMap *map);

/// Dump edge segments to file
void DumpEdgeSegments2File(EdgeMap *map, char *fname);

#endif