#pragma warning(disable:4996)  // For _CRT_SECURE_NO_WARNINGS

#include <stdio.h> 
#include <windows.h>
#include <math.h>
#include <limits.h>

#include <cvconfig.h>
#include <highgui.h>
#include <cxcore.h>

#include "EdgeDrawing.h"
#include "Timer.h"

/// Special defines
#define EDGE_VERTICAL   1
#define EDGE_HORIZONTAL 2

#define TEMP_PIXEL    253
#define ANCHOR_PIXEL  254
#define EDGE_PIXEL    255

#define LEFT  1
#define RIGHT 2
#define UP    3
#define DOWN  4

struct StackNode {
  int r, c;   // starting pixel
  int parent; // parent chain (-1 if no parent)
  int dir;    // direction where you are supposed to go
};

// Used during Edge Linking
struct TmpChain {
  short dir;                   // Direction of the chain
  unsigned short len;          // # of pixels in the chain
  short parent;                // Parent of this node (-1 if no parent)
  short children[2];           // Children of this node (-1 if no children)
  Pixel *pixels;               // Pointer to the beginning of the pixels array
};

static void ValidateSegments(EDImage *ed);
static void ExtractNewSegments(EDImage *ed);

///===========================================================================
/// EDImage Constructor
///
EDImage::EDImage(int _width, int _height){
  width = _width;
  height = _height;

  lowPassImg = new unsigned char[width*height];
  dirImg = new unsigned char[width*height];
  edgeImg = new unsigned char[width*height];

  gradImg = new short[width*height];
  gXImage = new short[width*height];
  gYImage = new short[width*height];

  map = new EdgeMap(width*height/50);   // # of segments depends on the image size

  // Temporary arrays
  tmpPixels = new Pixel[width*height];   
  tmpChains = new TmpChain[width*height];
  tmpStack = new StackNode[width*height];

  // This is where the result is returned
  pixels = new Pixel[width*height];

  // Parameter initializations
//  GRADIENT_OPERATOR = SOBEL_OPERATOR;
  GRADIENT_OPERATOR = PREWITT_OPERATOR;
  GRADIENT_THRESH = 24;
  ANCHOR_THRESH = 0;
  SCAN_INTERVAL = 1;
  MIN_PATH_LEN = 10;

  noImages = 0;
  gaussTime = 0;
  gradientTime = 0;
  anchorTime = 0;
  joinTime = 0;
  validationTime = 0;

  // Sorted anchor computation stuff
  useSortedAnchors = false;
  maxGradValue = 32*256;
  C = NULL;
  A = NULL;

  // Should be validate segments after detection?
  validateSegments = false;
} //end-EDImage

///----------------------------------------------
/// EDImage Destructor
///
EDImage::~EDImage(){
  delete lowPassImg;
  delete dirImg;
  delete edgeImg;

  delete gradImg;
  delete gXImage;
  delete gYImage;

  delete map;

  delete tmpPixels;
  delete tmpChains;
  delete tmpStack;

  delete pixels;

  if (C) delete C;
  if (A) delete A;
} //end-~EDImage

///------------------------------------------------------------------------------------------
/// Perform Gauss filter on "src" and store the result in "dst"
///
static void GaussFilter(unsigned char *src, unsigned char *dst, int width, int height){
#if 0
  for (int j=0; j<width; j++){
    dst[j] = src[j];                    // row=0
    dst[width+j] = src[width+j];        // row=1

    dst[(height-2)*width+j] = src[(height-2)*width+j];  // row=height-2
    dst[(height-1)*width+j] = src[(height-1)*width+j];  // row=height-1
  } //end-for

  for (int i=2; i<height-2; i++){
    dst[i*width] = src[i*width];        // column=0
    dst[i*width+1] = src[i*width+1];    // column=1

    dst[i*width+width-2] = src[i*width+width-2];  // column=width-2
    dst[i*width+width-1] = src[i*width+width-1];  // column=width-1
  } //end-for
#else
  memcpy(dst, src, width*height);
#endif

#if 0
  // 5x5 kernel
  // { 1, 4, 7, 4, 1}
  // { 4, 16, 26, 16, 4}
  // { 7, 26, 41, 26, 7}
  // { 4, 16, 26, 16, 4}
  // { 1, 4, 7, 4, 1}
  for (int i=2; i<height-2; i++){
    for (int j=2; j<width-2; j++){
      dst[i*width+j] = 
        (src[(i-2)*width+j-2] + 4*src[(i-2)*width+j-1] + 7*src[(i-2)*width+j] + 4*src[(i-2)*width+j+1] + src[(i-2)*width+j+2] +
         4*src[(i-1)*width+j-2] + 16*src[(i-1)*width+j-1] + 26*src[(i-1)*width+j] + 16*src[(i-1)*width+j+1] + 4*src[(i-1)*width+j+2] +
         7*src[i*width+j-2] + 26*src[i*width+j-1] + 41*src[i*width+j] + 26*src[i*width+j+1] + 7*src[i*width+j+2] +
         4*src[(i+1)*width+j-2] + 16*src[(i+1)*width+j-1] + 26*src[(i+1)*width+j] + 16*src[(i+1)*width+j+1] + 4*src[(i+1)*width+j+2] +
         src[(i+2)*width+j-2] + 4*src[(i+2)*width+j-1] + 7*src[(i+2)*width+j] + 4*src[(i+2)*width+j+1] + src[(i+2)*width+j+2] + 137)/273;      
    } //end-for
  } //end-for

#else
  // Another 5x5 kernel
  //  {2, 4, 5, 4, 2},
  //  {4, 9, 12, 9, 4},
  //  {5, 12, 15, 12, 5},
  //  {4, 9, 12, 9, 4},
  //  {2, 4, 5, 4, 2}
  for (int i=2; i<height-2; i++){
    for (int j=2; j<width-2; j++){
      dst[i*width+j] =  
        (2*src[(i-2)*width+j-2] + 4*src[(i-2)*width+j-1] + 5*src[(i-2)*width+j] + 4*src[(i-2)*width+j+1] + 2*src[(i-2)*width+j+2] +
         4*src[(i-1)*width+j-2] + 9*src[(i-1)*width+j-1] + 12*src[(i-1)*width+j] + 9*src[(i-1)*width+j+1] + 4*src[(i-1)*width+j+2] +
         5*src[i*width+j-2] + 12*src[i*width+j-1] + 15*src[i*width+j] + 12*src[i*width+j+1] + 5*src[i*width+j+2] +
         4*src[(i+1)*width+j-2] + 9*src[(i+1)*width+j-1] + 12*src[(i+1)*width+j] + 9*src[(i+1)*width+j+1] + 4*src[(i+1)*width+j+2] +
         2*src[(i+2)*width+j-2] + 4*src[(i+2)*width+j-1] + 5*src[(i+2)*width+j] + 4*src[(i+2)*width+j+1] + 2*src[(i+2)*width+j+2] + 80) / 159;
    } //end-for
  } //end-for
#endif
} //end-GaussFilter


//*************************************************************************************
// GRADIENT OPERATORS
//*************************************************************************************
///---------------------------------------------------------
/// LSD Operator
///
static void LSDOperator(EDImage *ed){
  int width = ed->width;
  int height = ed->height;

  unsigned char *lowPassImg = ed->lowPassImg;
  unsigned char *dirImg = ed->dirImg;
  short *gradImg = ed->gradImg;

  int GRADIENT_THRESH = ed->GRADIENT_THRESH;

  // Initialize gradient image for row = 0, row = height-1, column=0, column=width-1                   
  for (int j=0; j<width; j++){gradImg[j] = gradImg[(height-1)*width+j] = SHRT_MIN;}
  for (int i=1; i<height-1; i++){gradImg[i*width] = gradImg[(i+1)*width-1] = SHRT_MIN;}

  // Compute the gradient image and edge directions for the rest of the pixels
  for (int i=1; i<height-1; i++){
    for (int j=1; j<width-1; j++){
      // LSD gradient computation
      // A B
      // C D
      // gx = (B-A) + (D-C)
      // gy = (C-A) + (D-B)
      //
      // To make this faster: 
      // com1 = (D-A)
      // com2 = (B-C)
      // Then: gx = com1 + com2 = (D-A) + (B-C) = (B-A) + (D-C)
      //       gy = com1 - com2 = (D-A) - (B-C) = (C-A) + (D-B)
      // 
      int com1 = lowPassImg[(i+1)*width+j+1] - lowPassImg[i*width+j];
      int com2 = lowPassImg[i*width+j+1] - lowPassImg[(i+1)*width+j];

      int gx = abs(com1 + com2);
      int gy = abs(com1 - com2);

      int sum = gx + gy;      
      int index = i*width+j;
      gradImg[index] = sum;

      if (sum >= GRADIENT_THRESH){
        if (gx >= gy) dirImg[index] = EDGE_VERTICAL;
        else          dirImg[index] = EDGE_HORIZONTAL;
      } //end-if

    } // end-for
  } // end-for
} //end-LSDOperator

///---------------------------------------------------------
/// Prewitt Operator
///
static void PrewittOperator(EDImage *ed){
  int width = ed->width;
  int height = ed->height;

  unsigned char *lowPassImg = ed->lowPassImg;
  unsigned char *dirImg = ed->dirImg;
  short *gradImg = ed->gradImg;

  int GRADIENT_THRESH = ed->GRADIENT_THRESH;

  // Initialize gradient image for row = 0, row = height-1, column=0, column=width-1                   
  for (int j=0; j<width; j++){gradImg[j] = gradImg[(height-1)*width+j] = SHRT_MIN;}
  for (int i=1; i<height-1; i++){gradImg[i*width] = gradImg[(i+1)*width-1] = SHRT_MIN;}

  // Compute the gradient image and edge directions for the rest of the pixels
  for (int i=1; i<height-1; i++){
    for (int j=1; j<width-1; j++){
      // Prewitt Operator in horizontal and vertical direction
      // A B C
      // D x E
      // F G H
      // gx = (C-A) + (E-D) + (H-F)
      // gy = (F-A) + (G-B) + (H-C)
      //
      // To make this faster: 
      // com1 = (H-A)
      // com2 = (C-F)
      // Then: gx = com1 + com2 + (E-D) = (H-A) + (C-F) + (E-D) = (C-A) + (E-D) + (H-F)
      //       gy = com1 - com2 + (G-B) = (H-A) - (C-F) + (G-B) = (F-A) + (G-B) + (H-C)
      // 
      int com1 = lowPassImg[(i+1)*width+j+1] - lowPassImg[(i-1)*width+j-1];
      int com2 = lowPassImg[(i-1)*width+j+1] - lowPassImg[(i+1)*width+j-1];

      int gx = abs(com1 + com2 + (lowPassImg[i*width+j+1] - lowPassImg[i*width+j-1]));
      int gy = abs(com1 - com2 + (lowPassImg[(i+1)*width+j] - lowPassImg[(i-1)*width+j]));

      int sum = gx + gy;      
      int index = i*width+j;
      gradImg[index] = sum;

      if (sum >= GRADIENT_THRESH){
        if (gx >= gy) dirImg[index] = EDGE_VERTICAL;
        else          dirImg[index] = EDGE_HORIZONTAL;
      } //end-if

    } // end-for
  } // end-for
} //end-PrewittOperator

///------------------------------------------------------------------------------------
/// Performs Sobel Operator on the image
///
static void SobelOperator(EDImage *ed){
  int width = ed->width;
  int height = ed->height;

  unsigned char *lowPassImg = ed->lowPassImg;
  unsigned char *dirImg = ed->dirImg;
  short *gradImg = ed->gradImg;
  short *gXImage = ed->gXImage;
  short *gYImage = ed->gYImage;

  int GRADIENT_THRESH = ed->GRADIENT_THRESH;

  // Initialize gradient image for row = 0, row = height-1, column=0, column=width-1                   
  for (int j=0; j<width; j++){gradImg[j] = gradImg[(height-1)*width+j] = SHRT_MIN;}
  for (int i=1; i<height-1; i++){gradImg[i*width] = gradImg[(i+1)*width-1] = SHRT_MIN;}

  // Compute the gradient image and edge directions for the rest of the pixels
  for (int i=1; i<height-1; i++){
    for (int j=1; j<width-1; j++){
      // Sobel Operator in horizontal and vertical direction
      // Faster method below
      // A B C
      // D x E
      // F G H
      // gx = (C-A) + 2*(E-D) + (H-F)
      // gy = (F-A) + 2*(G-B) + (H-C)
      //
      // To make this faster: 
      // com1 = (H-A)
      // com2 = (C-F)
      // Then: gx = com1 + com2 + 2*(E-D) = (H-A) + (C-F) + 2*(E-D) = (C-A) + 2*(E-D) + (H-F)
      //       gy = com1 - com2 + 2*(G-B) = (H-A) - (C-F) + 2*(G-B) = (F-A) + 2*(G-B) + (H-C)
      // 
//#define ORIGINAL
#ifdef ORIGINAL
	  int com1 = lowPassImg[(i+1)*width+j+1] - lowPassImg[(i-1)*width+j-1];
      int com2 = lowPassImg[(i-1)*width+j+1] - lowPassImg[(i+1)*width+j-1];

      int gx = abs(com1 + com2 + 2*(lowPassImg[i*width+j+1] - lowPassImg[i*width+j-1]));
      int gy = abs(com1 - com2 + 2*(lowPassImg[(i+1)*width+j] - lowPassImg[(i-1)*width+j]));

      int sum = gx + gy;      
      int index = i*width+j;
      gradImg[index] = sum;

      if (sum >= GRADIENT_THRESH){
        if (gx >= gy) dirImg[index] = EDGE_VERTICAL;
        else          dirImg[index] = EDGE_HORIZONTAL;
	  }
#else

	  int com1 = lowPassImg[(i+1)*width+j+1] - lowPassImg[(i-1)*width+j-1];
      int com2 = lowPassImg[(i-1)*width+j+1] - lowPassImg[(i+1)*width+j-1];

      int gx = com1 + com2 + 2*(lowPassImg[i*width+j+1] - lowPassImg[i*width+j-1]);
      int gy = com1 - com2 + 2*(lowPassImg[(i+1)*width+j] - lowPassImg[(i-1)*width+j]);

	  int index = i*width+j;
      int sum = abs(gx) + abs(gy);      

      gradImg[index] = sum;
	  gXImage[index] = gx; 
	  gYImage[index] = gy;

      if (sum >= GRADIENT_THRESH){
        if (abs(gx) >= abs(gy)) dirImg[index] = EDGE_VERTICAL;
        else          dirImg[index] = EDGE_HORIZONTAL;
      } //end-if
#endif

    } // end-for
  } // end-for
} //end-SobelOperator

///------------------------------------------------------------------------------------
/// Performs Scharr Operator on the image
///
static void ScharrOperator(EDImage *ed){
  int width = ed->width;
  int height = ed->height;

  unsigned char *lowPassImg = ed->lowPassImg;
  unsigned char *dirImg = ed->dirImg;
  short *gradImg = ed->gradImg;

  int GRADIENT_THRESH = ed->GRADIENT_THRESH;

  // Initialize gradient image for row = 0, row = height-1, column=0, column=width-1                   
  for (int j=0; j<width; j++){gradImg[j] = gradImg[(height-1)*width+j] = SHRT_MIN;}
  for (int i=1; i<height-1; i++){gradImg[i*width] = gradImg[(i+1)*width-1] = SHRT_MIN;}

  // Compute the gradient image and edge directions for the rest of the pixels
  for (int i=1; i<height-1; i++){
    for (int j=1; j<width-1; j++){
      // Scharr Operator in horizontal and vertical direction
      // A B C
      // D x E
      // F G H
      // gx = 3*(C-A) + 10*(E-D) + 3*(H-F)
      // gy = 3*(F-A) + 10*(G-B) + 3*(H-C)
      //
      // To make this faster: 
      // com1 = (H-A)
      // com2 = (C-F)
      // Then: gx = 3*(com1 + com2) + 10*(E-D) = 3*(H-A) + 3*(C-F) + 10*(E-D) = 3*(C-A) + 10*(E-D) + 3*(H-F)
      //       gy = 3*(com1 - com2) + 10*(G-B) = 3*(H-A) - 3*(C-F) + 10*(G-B) = 3*(F-A) + 10*(G-B) + 3*(H-C)
      // 
      int com1 = lowPassImg[(i+1)*width+j+1] - lowPassImg[(i-1)*width+j-1];
      int com2 = lowPassImg[(i-1)*width+j+1] - lowPassImg[(i+1)*width+j-1];

      int gx = abs(3*(com1 + com2) + 10*(lowPassImg[i*width+j+1] - lowPassImg[i*width+j-1]));
      int gy = abs(3*(com1 - com2) + 10*(lowPassImg[(i+1)*width+j] - lowPassImg[(i-1)*width+j]));

      int sum = gx + gy;      
      int index = i*width+j;
      gradImg[index] = sum;

      if (sum >= GRADIENT_THRESH){
        if (gx >= gy) dirImg[index] = EDGE_VERTICAL;
        else          dirImg[index] = EDGE_HORIZONTAL;
      } //end-if

    } // end-for
  } // end-for
} //end-ScharrOperator

///-----------------------------------------------------------------------------------
/// Compute anchor points
///
static void ComputeAnchorPoints(EDImage *ed){
  int width = ed->width;
  int height = ed->height;

  unsigned char *lowPassImg = ed->lowPassImg;
  unsigned char *dirImg = ed->dirImg;
  unsigned char *edgeImg = ed->edgeImg;
  short *gradImg = ed->gradImg;

  int ANCHOR_THRESH = ed->ANCHOR_THRESH;
  int GRADIENT_THRESH = ed->GRADIENT_THRESH;
  int SCAN_INTERVAL = ed->SCAN_INTERVAL;

  memset(edgeImg, 0, width*height);

  for (int i=2; i<height-2; i++){
    int start = 2;
    int inc = 1;
    if (i%SCAN_INTERVAL != 0){start=SCAN_INTERVAL; inc=SCAN_INTERVAL;}

   for (int j=start; j<width-2; j+=inc){
      if (gradImg[i*width+j] < GRADIENT_THRESH) continue;

      if (dirImg[i*width+j] == EDGE_VERTICAL){
        // vertical edge
        int diff1 = gradImg[i*width+j] - gradImg[i*width+j-1];
        int diff2 = gradImg[i*width+j] - gradImg[i*width+j+1];
        if (diff1 >= ANCHOR_THRESH && diff2 >= ANCHOR_THRESH) edgeImg[i*width+j] = ANCHOR_PIXEL;
        
      } else {
        // horizontal edge
        int diff1 = gradImg[i*width+j] - gradImg[(i-1)*width+j];
        int diff2 = gradImg[i*width+j] - gradImg[(i+1)*width+j];
        if (diff1 >= ANCHOR_THRESH && diff2 >= ANCHOR_THRESH) edgeImg[i*width+j] = ANCHOR_PIXEL;
      } // end-else
    } //end-for-inner
  } //end-for-outer
} //end-ComputeAnchorPoints


///-----------------------------------------------------------------------------------
/// Compute anchor points2
///
//#include "ImageIO.h"
static void ComputeAnchorPoints2(EDImage *ed){
  int width = ed->width;
  int height = ed->height;

  unsigned char *cannyImg;
/*
  if (ReadImagePGM("D:/VisionWS/Test/OutputImages/CannyEdgeImg.pgm", (char **)&cannyImg, &width, &height) == 0){
    printf("Failed opening D:/VisionWS/Test/OutputImages/CannyEdgeImg.pgm\n");
    exit(0);
  } //end-if
*/
  unsigned char *edgeImg = ed->edgeImg;
  memset(edgeImg, 0, width*height);
  short *gradImg = ed->gradImg;

  int GRADIENT_THRESH = ed->GRADIENT_THRESH;

  for (int i=2; i<height-2; i++){
   for (int j=2; j<width-2; j++){
      if (gradImg[i*width+j] < GRADIENT_THRESH) continue;
  
      if (cannyImg[i*width+j]) edgeImg[i*width+j] = ANCHOR_PIXEL;
    } //end-for-inner
  } //end-for-outer
} //end-ComputeAnchorPoints2

///--------------------------------------------------------------------------
/// Computes anchor offsets & sorts them
/// 
static void SortAnchorsByGradValue(EDImage *ed){
  int width = ed->width;
  int height = ed->height;

  if (ed->C) delete ed->C;
  int SIZE = ed->maxGradValue;  
  int *C = ed->C = new int[SIZE];
  memset(C, 0, sizeof(int)*SIZE);  

  // Count the number of grad values
  for (int i=1; i<height-1; i++){
    for (int j=1; j<width-1; j++){
      if (ed->edgeImg[i*width+j] != ANCHOR_PIXEL) continue;

      int grad = ed->gradImg[i*width+j];
      C[grad]++;
    } //end-for
  } //end-for  

  // Compute indices
  for (int i=1; i<SIZE; i++) C[i] += C[i-1];

  int noAnchors = ed->noAnchors = C[SIZE-1];
  if (ed->A) delete ed->A;
  int *A = ed->A = new int[noAnchors];
  memset(A, 0, sizeof(int)*noAnchors);

  for (int i=1; i<height-1; i++){
    for (int j=1; j<width-1; j++){
      if (ed->edgeImg[i*width+j] != ANCHOR_PIXEL) continue;

      int grad = ed->gradImg[i*width+j];
      int index = --C[grad];
      A[index] = i*width+j;    // anchor's offset 
    } //end-for
  } //end-for  
} //end-SortAnchorsByGradValue

///--------------------------------------------------------------------------
/// Computes the length of the longest chain
///
static int LongestChain(TmpChain *chains, int root){
  if (root == -1 || chains[root].len == 0) return 0;

  int len0 = 0;
  if (chains[root].children[0] != -1) len0 = LongestChain(chains, chains[root].children[0]);

  int len1 = 0;
  if (chains[root].children[1] != -1) len1 = LongestChain(chains, chains[root].children[1]);

  int max = 0;

  if (len0 >= len1){
    max = len0;
    chains[root].children[1] = -1;

  } else {
    max = len1;
    chains[root].children[0] = -1;
  } //end-else

  return chains[root].len + max;
} //end-LongestChain

///-----------------------------------------------------------
/// Retrieves the chain nos from the tree
///
static int RetrieveChainNos(TmpChain *chains, int root, int chainNos[]){
  int count = 0;

  while (root != -1){
    chainNos[count] = root;
    count++;

    if (chains[root].children[0] != -1) root = chains[root].children[0];
    else                                root = chains[root].children[1];
  } //end-while

  return count;
} //end-RetrieveChainNos

///-----------------------------------------------------------------------------------
/// Join anchor points and compute segment chains at the same time
///
static void JoinAnchorPoints(EDImage *ed){
  int width = ed->width;
  int height = ed->height;

  unsigned char *dirImg = ed->dirImg;
  unsigned char *edgeImg = ed->edgeImg;
  short *gradImg = ed->gradImg;
  int *chainNos = new int[10000];

  int GRADIENT_THRESH = ed->GRADIENT_THRESH;
  int MIN_PATH_LEN = ed->MIN_PATH_LEN;

  Pixel *pixels = ed->tmpPixels;
  StackNode *stack = (StackNode *)ed->tmpStack;
  TmpChain *chains = (TmpChain *)ed->tmpChains;

  EdgeMap *map = ed->map;

  int noSegments = 0;  
  int totalPixels = 0;
  
  for (int i=2; i<height-2; i++){
    for (int j=2; j<width-2; j++){
      if (edgeImg[i*width+j] != ANCHOR_PIXEL) continue;   

      chains[0].len = 0;
      chains[0].parent = -1;
      chains[0].dir = 0;
      chains[0].children[0] = chains[0].children[1] = -1;
      chains[0].pixels = NULL;

      int noChains = 1;
      int len = 0;
      int duplicatePixelCount = 0;

      int top = -1;  // top of the stack 

      if (dirImg[i*width+j] == EDGE_VERTICAL){
        stack[++top].r = i;
        stack[top].c = j;
        stack[top].dir = DOWN;
        stack[top].parent = 0;

        stack[++top].r = i;
        stack[top].c = j;
        stack[top].dir = UP;
        stack[top].parent = 0;

      } else {
        stack[++top].r = i;
        stack[top].c = j;
        stack[top].dir = RIGHT;
        stack[top].parent = 0;

        stack[++top].r = i;
        stack[top].c = j;
        stack[top].dir = LEFT;
        stack[top].parent = 0;
      } //end-else

      // While the stack is not empty
StartOfWhile:
      while (top >= 0){
        int r = stack[top].r;
        int c = stack[top].c;
        int dir = stack[top].dir;
        int parent = stack[top].parent;
        top--;

        if (edgeImg[r*width+c] != EDGE_PIXEL) duplicatePixelCount++;

        
        chains[noChains].dir = dir;   // traversal direction
        chains[noChains].parent = parent;
        chains[noChains].children[0] = chains[noChains].children[1] = -1;

        int chainLen = 0;
        chains[noChains].pixels = &pixels[len];
        
        pixels[len].r = r;
        pixels[len].c = c;
        len++;
        chainLen++;
        
        if (dir == LEFT){
          while (dirImg[r*width+c] == EDGE_HORIZONTAL){
            edgeImg[r*width+c] = EDGE_PIXEL;

            // The edge is horizontal. Look LEFT
            //
            //   A
            //   B x 
            //   C 
            //
            // cleanup up & down pixels
            if (edgeImg[(r-1)*width+c] == ANCHOR_PIXEL) edgeImg[(r-1)*width+c] = 0;
            if (edgeImg[(r+1)*width+c] == ANCHOR_PIXEL) edgeImg[(r+1)*width+c] = 0;

            // Look if there is an edge pixel in the neighbors
            if (edgeImg[r*width+c-1] >= ANCHOR_PIXEL){c--;}
            else if (edgeImg[(r-1)*width+c-1] >= ANCHOR_PIXEL){r--; c--;}
            else if (edgeImg[(r+1)*width+c-1] >= ANCHOR_PIXEL){r++; c--;}
            else {
              // else -- follow max. pixel to the LEFT
              int A = gradImg[(r-1)*width+c-1];
              int B = gradImg[r*width+c-1];
              int C = gradImg[(r+1)*width+c-1];

              if (A > B){
                if (A > C) r--; 
                else       r++;
              } else  if (C > B) r++;
              c--;
            } //end-else

            if (edgeImg[r*width+c] == EDGE_PIXEL || gradImg[r*width+c] < GRADIENT_THRESH){
              if (chainLen > 0){
                chains[noChains].len = chainLen;
                if (dir == LEFT || dir == UP) chains[parent].children[0] = noChains;
                else                          chains[parent].children[1] = noChains;
                noChains++;
              } // end-if
              goto StartOfWhile;
            } //end-else

            pixels[len].r = r;
            pixels[len].c = c;
            len++;
            chainLen++;
          } //end-while

          stack[++top].r = r;
          stack[top].c = c;
          stack[top].dir = DOWN;
          stack[top].parent = noChains;
    
          stack[++top].r = r;
          stack[top].c = c;
          stack[top].dir = UP;
          stack[top].parent = noChains;

          len--;
          chainLen--;

          chains[noChains].len = chainLen;
          if (dir == LEFT || dir == UP) chains[parent].children[0] = noChains;
          else                          chains[parent].children[1] = noChains;
          noChains++;

        } else if (dir == RIGHT){
          while (dirImg[r*width+c] == EDGE_HORIZONTAL){
            edgeImg[r*width+c] = EDGE_PIXEL;

            // The edge is horizontal. Look RIGHT
            //
            //     A
            //   x B
            //     C
            //
            // cleanup up&down pixels
            if (edgeImg[(r+1)*width+c] == ANCHOR_PIXEL) edgeImg[(r+1)*width+c] = 0;
            if (edgeImg[(r-1)*width+c] == ANCHOR_PIXEL) edgeImg[(r-1)*width+c] = 0;

            // Look if there is an edge pixel in the neighbors
            if (edgeImg[r*width+c+1] >= ANCHOR_PIXEL){c++;}
            else if (edgeImg[(r+1)*width+c+1] >= ANCHOR_PIXEL){r++; c++;}
            else if (edgeImg[(r-1)*width+c+1] >= ANCHOR_PIXEL){r--; c++;}
            else {
              // else -- follow max. pixel to the RIGHT
              int A = gradImg[(r-1)*width+c+1];
              int B = gradImg[r*width+c+1];
              int C = gradImg[(r+1)*width+c+1];

              if (A > B){
                if (A > C) r--;       // A
                else       r++;       // C
              } else if (C > B) r++;  // C
              c++;
            } //end-else

            if (edgeImg[r*width+c] == EDGE_PIXEL || gradImg[r*width+c] < GRADIENT_THRESH){
              if (chainLen > 0){
                chains[noChains].len = chainLen;
                if (dir == LEFT || dir == UP) chains[parent].children[0] = noChains;
                else                          chains[parent].children[1] = noChains;
                noChains++;
              } // end-if
              goto StartOfWhile;
            } //end-else

            pixels[len].r = r;
            pixels[len].c = c;
            len++;
            chainLen++;
          } //end-while

          stack[++top].r = r;
          stack[top].c = c;
          stack[top].dir = DOWN;  // Go down
          stack[top].parent = noChains;

          stack[++top].r = r;
          stack[top].c = c;
          stack[top].dir = UP;   // Go up
          stack[top].parent = noChains;

          len--;
          chainLen--;

          chains[noChains].len = chainLen;
          if (dir == LEFT || dir == UP) chains[parent].children[0] = noChains;
          else                          chains[parent].children[1] = noChains;
          noChains++;

        } else if (dir == UP){
          while (dirImg[r*width+c] == EDGE_VERTICAL){
            edgeImg[r*width+c] = EDGE_PIXEL;

            // The edge is vertical. Look UP
            //
            //   A B C
            //     x
            //
            // Cleanup left & right pixels
            if (edgeImg[r*width+c-1] == ANCHOR_PIXEL) edgeImg[r*width+c-1] = 0;
            if (edgeImg[r*width+c+1] == ANCHOR_PIXEL) edgeImg[r*width+c+1] = 0;

            // Look if there is an edge pixel in the neighbors
            if (edgeImg[(r-1)*width+c] >= ANCHOR_PIXEL){r--;}
            else if (edgeImg[(r-1)*width+c-1] >= ANCHOR_PIXEL){r--; c--;}
            else if (edgeImg[(r-1)*width+c+1] >= ANCHOR_PIXEL){r--; c++;}
            else {
              // else -- follow the max. pixel UP
              int A = gradImg[(r-1)*width+c-1];
              int B = gradImg[(r-1)*width+c];
              int C = gradImg[(r-1)*width+c+1];

              if (A > B){
                if (A > C) c--;
                else       c++;
              } else if (C > B) c++;
              r--;
            } //end-else

            if (edgeImg[r*width+c] == EDGE_PIXEL || gradImg[r*width+c] < GRADIENT_THRESH){
              if (chainLen > 0){
                chains[noChains].len = chainLen;
                if (dir == LEFT || dir == UP) chains[parent].children[0] = noChains;
                else                          chains[parent].children[1] = noChains;
                noChains++;
              } // end-if
              goto StartOfWhile;
            } //end-else

            pixels[len].r = r;
            pixels[len].c = c;
            len++;
            chainLen++;
          } //end-while

          stack[++top].r = r;
          stack[top].c = c;
          stack[top].dir = RIGHT;
          stack[top].parent = noChains;

          stack[++top].r = r;
          stack[top].c = c;
          stack[top].dir = LEFT;
          stack[top].parent = noChains;

          len--;
          chainLen--;

          chains[noChains].len = chainLen;
          if (dir == LEFT || dir == UP) chains[parent].children[0] = noChains;
          else                          chains[parent].children[1] = noChains;
          noChains++;

        } else { // dir == DOWN
          while (dirImg[r*width+c] == EDGE_VERTICAL){
            edgeImg[r*width+c] = EDGE_PIXEL;

            // The edge is vertical
            //
            //     x
            //   A B C
            //
            // cleanup side pixels
            if (edgeImg[r*width+c+1] == ANCHOR_PIXEL) edgeImg[r*width+c+1] = 0;
            if (edgeImg[r*width+c-1] == ANCHOR_PIXEL) edgeImg[r*width+c-1] = 0;

            // Look if there is an edge pixel in the neighbors
            if (edgeImg[(r+1)*width+c] >= ANCHOR_PIXEL){r++;}
            else if (edgeImg[(r+1)*width+c+1] >= ANCHOR_PIXEL){r++; c++;}
            else if (edgeImg[(r+1)*width+c-1] >= ANCHOR_PIXEL){r++; c--;}
            else {
              // else -- follow the max. pixel DOWN
              int A = gradImg[(r+1)*width+c-1];
              int B = gradImg[(r+1)*width+c];
              int C = gradImg[(r+1)*width+c+1];

              if (A > B){
                if (A > C) c--;       // A
                else       c++;       // C
              } else if (C > B) c++;  // C
              r++;                    
            } //end-else

            if (edgeImg[r*width+c] == EDGE_PIXEL || gradImg[r*width+c] < GRADIENT_THRESH){
              if (chainLen > 0){
                chains[noChains].len = chainLen;
                if (dir == LEFT || dir == UP) chains[parent].children[0] = noChains;
                else                          chains[parent].children[1] = noChains;
                noChains++;
              } // end-if
              goto StartOfWhile;
            } //end-else

            pixels[len].r = r;
            pixels[len].c = c;
            len++;
            chainLen++;
          } //end-while

          stack[++top].r = r;
          stack[top].c = c;
          stack[top].dir = RIGHT;
          stack[top].parent = noChains;

          stack[++top].r = r;
          stack[top].c = c;
          stack[top].dir = LEFT;
          stack[top].parent = noChains;

          len--;
          chainLen--;

          chains[noChains].len = chainLen;
          if (dir == LEFT || dir == UP) chains[parent].children[0] = noChains;
          else                          chains[parent].children[1] = noChains;
          noChains++;
        } //end-else

      } //end-while

      if (len-duplicatePixelCount < MIN_PATH_LEN){
        for (int k=0; k<len; k++){
          edgeImg[pixels[k].r*width+pixels[k].c] = 0;          
        } //end-for

      } else {
        map->segments[noSegments].pixels = &ed->pixels[totalPixels];

        int totalLen = LongestChain(chains, chains[0].children[1]);
        int noSegmentPixels = 0;

        if (totalLen > 0){
          // Retrieve the chainNos
          int count = RetrieveChainNos(chains, chains[0].children[1], chainNos);

          // Copy these pixels in the reverse order
          for (int k=count-1; k>=0; k--){
            int chainNo = chainNos[k];

#if 1
            /* See if we can erase some pixels from the last chain. This is for cleanup */
            int fr = chains[chainNo].pixels[chains[chainNo].len-1].r;
            int fc = chains[chainNo].pixels[chains[chainNo].len-1].c;

            int index = noSegmentPixels-2;
            while (index >=0){
              int dr = abs(fr-map->segments[noSegments].pixels[index].r);
              int dc = abs(fc-map->segments[noSegments].pixels[index].c);

              if (dr<=1 && dc<=1){
                // neighbors. Erase last pixel
                noSegmentPixels--;
                index--;
              } else break;
            } //end-while

            if (chains[chainNo].len > 1){
              fr = chains[chainNo].pixels[chains[chainNo].len-2].r;
              fc = chains[chainNo].pixels[chains[chainNo].len-2].c;

              int dr = abs(fr-map->segments[noSegments].pixels[noSegmentPixels-1].r);
              int dc = abs(fc-map->segments[noSegments].pixels[noSegmentPixels-1].c);

              if (dr<=1 && dc<=1) chains[chainNo].len--;
            } //end-if
#endif

            for (int l=chains[chainNo].len-1; l>=0; l--){
              map->segments[noSegments].pixels[noSegmentPixels++] = chains[chainNo].pixels[l];
            } //end-for

            chains[chainNo].len = 0;  // Mark as copied
          } //end-for
        } //end-if

        totalLen = LongestChain(chains, chains[0].children[0]);
        if (totalLen > 1){
          // Retrieve the chainNos
          int count = RetrieveChainNos(chains, chains[0].children[0], chainNos);

          // Copy these chains in the forward direction. Skip the first pixel of the first chain
          // due to repetition with the last pixel of the previous chain
          int lastChainNo = chainNos[0];
          chains[lastChainNo].pixels++;
          chains[lastChainNo].len--;

          for (int k=0; k<count; k++){
            int chainNo = chainNos[k];

#if 1
            /* See if we can erase some pixels from the last chain. This is for cleanup */
            int fr = chains[chainNo].pixels[0].r;
            int fc = chains[chainNo].pixels[0].c;

            int index = noSegmentPixels-2;
            while (index >=0){
              int dr = abs(fr-map->segments[noSegments].pixels[index].r);
              int dc = abs(fc-map->segments[noSegments].pixels[index].c);

              if (dr<=1 && dc<=1){
                // neighbors. Erase last pixel
                noSegmentPixels--;
                index--;
              } else break;
            } //end-while

            int startIndex = 0;
            int chainLen = chains[chainNo].len;
            if (chainLen > 1){
              int fr = chains[chainNo].pixels[1].r;
              int fc = chains[chainNo].pixels[1].c;

              int dr = abs(fr-map->segments[noSegments].pixels[noSegmentPixels-1].r);
              int dc = abs(fc-map->segments[noSegments].pixels[noSegmentPixels-1].c);

              if (dr<=1 && dc<=1){startIndex = 1;}
            } //end-if
#endif

            /* Start a new chain & copy pixels from the new chain */
            for (int l=startIndex; l<chains[chainNo].len; l++){
              map->segments[noSegments].pixels[noSegmentPixels++] = chains[chainNo].pixels[l];
            } //end-for

            chains[chainNo].len = 0;  // Mark as copied
          } //end-for
        } //end-if

        map->segments[noSegments].noPixels = noSegmentPixels;
        totalPixels += noSegmentPixels;

        // See if the first pixel can be cleaned up
        int fr = map->segments[noSegments].pixels[1].r;
        int fc = map->segments[noSegments].pixels[1].c;

        int dr = abs(fr-map->segments[noSegments].pixels[noSegmentPixels-1].r);
        int dc = abs(fc-map->segments[noSegments].pixels[noSegmentPixels-1].c);

        if (dr<=1 && dc<=1){
          map->segments[noSegments].pixels++;
          map->segments[noSegments].noPixels--;
        } //end-if

        noSegments++;

        // Copy the rest of the long chains here
        for (int k=2; k<noChains; k++){
          if (chains[k].len < 2) continue;

          totalLen = LongestChain(chains, k);

          // If long enough, copy
//          if (totalLen >= 12){
          if (totalLen >= 10){
            map->segments[noSegments].pixels = &ed->pixels[totalPixels];
    
            // Retrieve the chainNos
            int count = RetrieveChainNos(chains, k, chainNos);

            // Copy the pixels
            noSegmentPixels = 0;
            for (int k=0; k<count; k++){
              int chainNo = chainNos[k];

#if 1
              /* See if we can erase some pixels from the last chain. This is for cleanup */
              int fr = chains[chainNo].pixels[0].r;
              int fc = chains[chainNo].pixels[0].c;

              int index = noSegmentPixels-2;
              while (index >=0){
                int dr = abs(fr-map->segments[noSegments].pixels[index].r);
                int dc = abs(fc-map->segments[noSegments].pixels[index].c);

                if (dr<=1 && dc<=1){
                  // neighbors. Erase last pixel
                  noSegmentPixels--;
                  index--;
                } else break;
              } //end-while

              int startIndex = 0;
              int chainLen = chains[chainNo].len;
              if (chainLen > 1){
                int fr = chains[chainNo].pixels[1].r;
                int fc = chains[chainNo].pixels[1].c;

                int dr = abs(fr-map->segments[noSegments].pixels[noSegmentPixels-1].r);
                int dc = abs(fc-map->segments[noSegments].pixels[noSegmentPixels-1].c);

                if (dr<=1 && dc<=1){startIndex = 1;}
              } //end-if
#endif
              /* Start a new chain & copy pixels from the new chain */
              for (int l=startIndex; l<chains[chainNo].len; l++){
                map->segments[noSegments].pixels[noSegmentPixels++] = chains[chainNo].pixels[l];
              } //end-for

              chains[chainNo].len = 0;  // Mark as copied
            } //end-for

            map->segments[noSegments].noPixels = noSegmentPixels;
            totalPixels += noSegmentPixels;        

            noSegments++;
          } //end-if          
        } //end-for

      } //end-else

    } //end-for-inner
  } //end-for-outer

  map->noSegments = noSegments;

  delete chainNos;
} //end-JoinAnchorPoints

///-----------------------------------------------------------------------------------
/// Join anchors starting with the anchor having the maximum gradient value.
/// To do this, we need to first sort the anchors
///
static void JoinAnchorPointsUsingSortedAnchors(EDImage *ed){
  int width = ed->width;
  int height = ed->height;

  int GRADIENT_THRESH = ed->GRADIENT_THRESH;
  int MIN_PATH_LEN = ed->MIN_PATH_LEN;

  unsigned char *dirImg = ed->dirImg;
  unsigned char *edgeImg = ed->edgeImg;
  short *gradImg = ed->gradImg;
  int *chainNos = new int[10000];

  Pixel *pixels = ed->tmpPixels;
  StackNode *stack = (StackNode *)ed->tmpStack;
  TmpChain *chains = (TmpChain *)ed->tmpChains;

  EdgeMap *map = ed->map;

  // First Sort the anchors
  SortAnchorsByGradValue(ed);

  int noSegments = 0;  
  int totalPixels = 0;

  for (int k=ed->noAnchors-1; k>=0; k--){
    int pixelOffset = ed->A[k];

    int i = pixelOffset/width;
    int j = pixelOffset % width;

    if (edgeImg[i*width+j] != ANCHOR_PIXEL) continue;   

      chains[0].len = 0;
      chains[0].parent = -1;
      chains[0].dir = 0;
      chains[0].children[0] = chains[0].children[1] = -1;
      chains[0].pixels = NULL;

      int noChains = 1;
      int len = 0;
      int duplicatePixelCount = 0;

      int top = -1;  // top of the stack 

      if (dirImg[i*width+j] == EDGE_VERTICAL){
        stack[++top].r = i;
        stack[top].c = j;
        stack[top].dir = DOWN;
        stack[top].parent = 0;

        stack[++top].r = i;
        stack[top].c = j;
        stack[top].dir = UP;
        stack[top].parent = 0;

      } else {
        stack[++top].r = i;
        stack[top].c = j;
        stack[top].dir = RIGHT;
        stack[top].parent = 0;

        stack[++top].r = i;
        stack[top].c = j;
        stack[top].dir = LEFT;
        stack[top].parent = 0;
      } //end-else

      // While the stack is not empty
StartOfWhile:
      while (top >= 0){
        int r = stack[top].r;
        int c = stack[top].c;
        int dir = stack[top].dir;
        int parent = stack[top].parent;
        top--;

        if (edgeImg[r*width+c] != EDGE_PIXEL) duplicatePixelCount++;

        
        chains[noChains].dir = dir;   // traversal direction
        chains[noChains].parent = parent;
        chains[noChains].children[0] = chains[noChains].children[1] = -1;

        int chainLen = 0;
        chains[noChains].pixels = &pixels[len];
        
        pixels[len].r = r;
        pixels[len].c = c;
        len++;
        chainLen++;
        
        if (dir == LEFT){
          while (dirImg[r*width+c] == EDGE_HORIZONTAL){
            edgeImg[r*width+c] = EDGE_PIXEL;

            // The edge is horizontal. Look LEFT
            //
            //   A
            //   B x 
            //   C 
            //
            // cleanup up & down pixels
            if (edgeImg[(r-1)*width+c] == ANCHOR_PIXEL) edgeImg[(r-1)*width+c] = 0;
            if (edgeImg[(r+1)*width+c] == ANCHOR_PIXEL) edgeImg[(r+1)*width+c] = 0;

            // Look if there is an edge pixel in the neighbors
            if (edgeImg[r*width+c-1] >= ANCHOR_PIXEL){c--;}
            else if (edgeImg[(r-1)*width+c-1] >= ANCHOR_PIXEL){r--; c--;}
            else if (edgeImg[(r+1)*width+c-1] >= ANCHOR_PIXEL){r++; c--;}
            else {
              // else -- follow max. pixel to the LEFT
              int A = gradImg[(r-1)*width+c-1];
              int B = gradImg[r*width+c-1];
              int C = gradImg[(r+1)*width+c-1];

              if (A > B){
                if (A > C) r--; 
                else       r++;
              } else  if (C > B) r++;
              c--;
            } //end-else

            if (edgeImg[r*width+c] == EDGE_PIXEL || gradImg[r*width+c] < GRADIENT_THRESH){
              if (chainLen > 0){
                chains[noChains].len = chainLen;
                if (dir == LEFT || dir == UP) chains[parent].children[0] = noChains;
                else                          chains[parent].children[1] = noChains;
                noChains++;
              } // end-if
              goto StartOfWhile;
            } //end-else

            pixels[len].r = r;
            pixels[len].c = c;
            len++;
            chainLen++;
          } //end-while

          stack[++top].r = r;
          stack[top].c = c;
          stack[top].dir = DOWN;
          stack[top].parent = noChains;
    
          stack[++top].r = r;
          stack[top].c = c;
          stack[top].dir = UP;
          stack[top].parent = noChains;

          len--;
          chainLen--;

          chains[noChains].len = chainLen;
          if (dir == LEFT || dir == UP) chains[parent].children[0] = noChains;
          else                          chains[parent].children[1] = noChains;
          noChains++;

        } else if (dir == RIGHT){
          while (dirImg[r*width+c] == EDGE_HORIZONTAL){
            edgeImg[r*width+c] = EDGE_PIXEL;

            // The edge is horizontal. Look RIGHT
            //
            //     A
            //   x B
            //     C
            //
            // cleanup up&down pixels
            if (edgeImg[(r+1)*width+c] == ANCHOR_PIXEL) edgeImg[(r+1)*width+c] = 0;
            if (edgeImg[(r-1)*width+c] == ANCHOR_PIXEL) edgeImg[(r-1)*width+c] = 0;

            // Look if there is an edge pixel in the neighbors
            if (edgeImg[r*width+c+1] >= ANCHOR_PIXEL){c++;}
            else if (edgeImg[(r+1)*width+c+1] >= ANCHOR_PIXEL){r++; c++;}
            else if (edgeImg[(r-1)*width+c+1] >= ANCHOR_PIXEL){r--; c++;}
            else {
              // else -- follow max. pixel to the RIGHT
              int A = gradImg[(r-1)*width+c+1];
              int B = gradImg[r*width+c+1];
              int C = gradImg[(r+1)*width+c+1];

              if (A > B){
                if (A > C) r--;       // A
                else       r++;       // C
              } else if (C > B) r++;  // C
              c++;
            } //end-else

            if (edgeImg[r*width+c] == EDGE_PIXEL || gradImg[r*width+c] < GRADIENT_THRESH){
              if (chainLen > 0){
                chains[noChains].len = chainLen;
                if (dir == LEFT || dir == UP) chains[parent].children[0] = noChains;
                else                          chains[parent].children[1] = noChains;
                noChains++;
              } // end-if
              goto StartOfWhile;
            } //end-else

            pixels[len].r = r;
            pixels[len].c = c;
            len++;
            chainLen++;
          } //end-while

          stack[++top].r = r;
          stack[top].c = c;
          stack[top].dir = DOWN;  // Go down
          stack[top].parent = noChains;

          stack[++top].r = r;
          stack[top].c = c;
          stack[top].dir = UP;   // Go up
          stack[top].parent = noChains;

          len--;
          chainLen--;

          chains[noChains].len = chainLen;
          if (dir == LEFT || dir == UP) chains[parent].children[0] = noChains;
          else                          chains[parent].children[1] = noChains;
          noChains++;

        } else if (dir == UP){
          while (dirImg[r*width+c] == EDGE_VERTICAL){
            edgeImg[r*width+c] = EDGE_PIXEL;

            // The edge is vertical. Look UP
            //
            //   A B C
            //     x
            //
            // Cleanup left & right pixels
            if (edgeImg[r*width+c-1] == ANCHOR_PIXEL) edgeImg[r*width+c-1] = 0;
            if (edgeImg[r*width+c+1] == ANCHOR_PIXEL) edgeImg[r*width+c+1] = 0;

            // Look if there is an edge pixel in the neighbors
            if (edgeImg[(r-1)*width+c] >= ANCHOR_PIXEL){r--;}
            else if (edgeImg[(r-1)*width+c-1] >= ANCHOR_PIXEL){r--; c--;}
            else if (edgeImg[(r-1)*width+c+1] >= ANCHOR_PIXEL){r--; c++;}
            else {
              // else -- follow the max. pixel UP
              int A = gradImg[(r-1)*width+c-1];
              int B = gradImg[(r-1)*width+c];
              int C = gradImg[(r-1)*width+c+1];

              if (A > B){
                if (A > C) c--;
                else       c++;
              } else if (C > B) c++;
              r--;
            } //end-else

            if (edgeImg[r*width+c] == EDGE_PIXEL || gradImg[r*width+c] < GRADIENT_THRESH){
              if (chainLen > 0){
                chains[noChains].len = chainLen;
                if (dir == LEFT || dir == UP) chains[parent].children[0] = noChains;
                else                          chains[parent].children[1] = noChains;
                noChains++;
              } // end-if
              goto StartOfWhile;
            } //end-else

            pixels[len].r = r;
            pixels[len].c = c;
            len++;
            chainLen++;
          } //end-while

          stack[++top].r = r;
          stack[top].c = c;
          stack[top].dir = RIGHT;
          stack[top].parent = noChains;

          stack[++top].r = r;
          stack[top].c = c;
          stack[top].dir = LEFT;
          stack[top].parent = noChains;

          len--;
          chainLen--;

          chains[noChains].len = chainLen;
          if (dir == LEFT || dir == UP) chains[parent].children[0] = noChains;
          else                          chains[parent].children[1] = noChains;
          noChains++;

        } else { // dir == DOWN
          while (dirImg[r*width+c] == EDGE_VERTICAL){
            edgeImg[r*width+c] = EDGE_PIXEL;

            // The edge is vertical
            //
            //     x
            //   A B C
            //
            // cleanup side pixels
            if (edgeImg[r*width+c+1] == ANCHOR_PIXEL) edgeImg[r*width+c+1] = 0;
            if (edgeImg[r*width+c-1] == ANCHOR_PIXEL) edgeImg[r*width+c-1] = 0;

            // Look if there is an edge pixel in the neighbors
            if (edgeImg[(r+1)*width+c] >= ANCHOR_PIXEL){r++;}
            else if (edgeImg[(r+1)*width+c+1] >= ANCHOR_PIXEL){r++; c++;}
            else if (edgeImg[(r+1)*width+c-1] >= ANCHOR_PIXEL){r++; c--;}
            else {
              // else -- follow the max. pixel DOWN
              int A = gradImg[(r+1)*width+c-1];
              int B = gradImg[(r+1)*width+c];
              int C = gradImg[(r+1)*width+c+1];

              if (A > B){
                if (A > C) c--;       // A
                else       c++;       // C
              } else if (C > B) c++;  // C
              r++;                    
            } //end-else

            if (edgeImg[r*width+c] == EDGE_PIXEL || gradImg[r*width+c] < GRADIENT_THRESH){
              if (chainLen > 0){
                chains[noChains].len = chainLen;
                if (dir == LEFT || dir == UP) chains[parent].children[0] = noChains;
                else                          chains[parent].children[1] = noChains;
                noChains++;
              } // end-if
              goto StartOfWhile;
            } //end-else

            pixels[len].r = r;
            pixels[len].c = c;
            len++;
            chainLen++;
          } //end-while

          stack[++top].r = r;
          stack[top].c = c;
          stack[top].dir = RIGHT;
          stack[top].parent = noChains;

          stack[++top].r = r;
          stack[top].c = c;
          stack[top].dir = LEFT;
          stack[top].parent = noChains;

          len--;
          chainLen--;

          chains[noChains].len = chainLen;
          if (dir == LEFT || dir == UP) chains[parent].children[0] = noChains;
          else                          chains[parent].children[1] = noChains;
          noChains++;
        } //end-else

      } //end-while

      if (len-duplicatePixelCount < MIN_PATH_LEN){
        for (int k=0; k<len; k++){
          edgeImg[pixels[k].r*width+pixels[k].c] = 0;          
        } //end-for

      } else {
        map->segments[noSegments].pixels = &ed->pixels[totalPixels];

        int totalLen = LongestChain(chains, chains[0].children[1]);
        int noSegmentPixels = 0;

        if (totalLen > 0){
          // Retrieve the chainNos
          int count = RetrieveChainNos(chains, chains[0].children[1], chainNos);

          // Copy these pixels in the reverse order
          for (int k=count-1; k>=0; k--){
            int chainNo = chainNos[k];

#if 1
            /* See if we can erase some pixels from the last chain. This is for cleanup */
            int fr = chains[chainNo].pixels[chains[chainNo].len-1].r;
            int fc = chains[chainNo].pixels[chains[chainNo].len-1].c;

            int index = noSegmentPixels-2;
            while (index >=0){
              int dr = abs(fr-map->segments[noSegments].pixels[index].r);
              int dc = abs(fc-map->segments[noSegments].pixels[index].c);

              if (dr<=1 && dc<=1){
                // neighbors. Erase last pixel
                noSegmentPixels--;
                index--;
              } else break;
            } //end-while

            if (chains[chainNo].len > 1){
              fr = chains[chainNo].pixels[chains[chainNo].len-2].r;
              fc = chains[chainNo].pixels[chains[chainNo].len-2].c;

              int dr = abs(fr-map->segments[noSegments].pixels[noSegmentPixels-1].r);
              int dc = abs(fc-map->segments[noSegments].pixels[noSegmentPixels-1].c);

              if (dr<=1 && dc<=1) chains[chainNo].len--;
            } //end-if
#endif

            for (int l=chains[chainNo].len-1; l>=0; l--){
              map->segments[noSegments].pixels[noSegmentPixels++] = chains[chainNo].pixels[l];
            } //end-for

            chains[chainNo].len = 0;  // Mark as copied
          } //end-for
        } //end-if

        totalLen = LongestChain(chains, chains[0].children[0]);
        if (totalLen > 1){
          // Retrieve the chainNos
          int count = RetrieveChainNos(chains, chains[0].children[0], chainNos);

          // Copy these chains in the forward direction. Skip the first pixel of the first chain
          // due to repetition with the last pixel of the previous chain
          int lastChainNo = chainNos[0];
          chains[lastChainNo].pixels++;
          chains[lastChainNo].len--;

          for (int k=0; k<count; k++){
            int chainNo = chainNos[k];

#if 1
            /* See if we can erase some pixels from the last chain. This is for cleanup */
            int fr = chains[chainNo].pixels[0].r;
            int fc = chains[chainNo].pixels[0].c;

            int index = noSegmentPixels-2;
            while (index >=0){
              int dr = abs(fr-map->segments[noSegments].pixels[index].r);
              int dc = abs(fc-map->segments[noSegments].pixels[index].c);

              if (dr<=1 && dc<=1){
                // neighbors. Erase last pixel
                noSegmentPixels--;
                index--;
              } else break;
            } //end-while

            int startIndex = 0;
            int chainLen = chains[chainNo].len;
            if (chainLen > 1){
              int fr = chains[chainNo].pixels[1].r;
              int fc = chains[chainNo].pixels[1].c;

              int dr = abs(fr-map->segments[noSegments].pixels[noSegmentPixels-1].r);
              int dc = abs(fc-map->segments[noSegments].pixels[noSegmentPixels-1].c);

              if (dr<=1 && dc<=1){startIndex = 1;}
            } //end-if
#endif

            /* Start a new chain & copy pixels from the new chain */
            for (int l=startIndex; l<chains[chainNo].len; l++){
              map->segments[noSegments].pixels[noSegmentPixels++] = chains[chainNo].pixels[l];
            } //end-for

            chains[chainNo].len = 0;  // Mark as copied
          } //end-for
        } //end-if

        map->segments[noSegments].noPixels = noSegmentPixels;
        totalPixels += noSegmentPixels;

        // See if the first pixel can be cleaned up
        int fr = map->segments[noSegments].pixels[1].r;
        int fc = map->segments[noSegments].pixels[1].c;

        int dr = abs(fr-map->segments[noSegments].pixels[noSegmentPixels-1].r);
        int dc = abs(fc-map->segments[noSegments].pixels[noSegmentPixels-1].c);

        if (dr<=1 && dc<=1){
          map->segments[noSegments].pixels++;
          map->segments[noSegments].noPixels--;
        } //end-if

        noSegments++;

        // Copy the rest of the long chains here
        for (int k=2; k<noChains; k++){
          if (chains[k].len < 2) continue;

          totalLen = LongestChain(chains, k);

          // If long enough, copy
//          if (totalLen >= 12){
          if (totalLen >= 10){
            map->segments[noSegments].pixels = &ed->pixels[totalPixels];
    
            // Retrieve the chainNos
            int count = RetrieveChainNos(chains, k, chainNos);

            // Copy the pixels
            noSegmentPixels = 0;
            for (int k=0; k<count; k++){
              int chainNo = chainNos[k];

#if 1
              /* See if we can erase some pixels from the last chain. This is for cleanup */
              int fr = chains[chainNo].pixels[0].r;
              int fc = chains[chainNo].pixels[0].c;

              int index = noSegmentPixels-2;
              while (index >=0){
                int dr = abs(fr-map->segments[noSegments].pixels[index].r);
                int dc = abs(fc-map->segments[noSegments].pixels[index].c);

                if (dr<=1 && dc<=1){
                  // neighbors. Erase last pixel
                  noSegmentPixels--;
                  index--;
                } else break;
              } //end-while

              int startIndex = 0;
              int chainLen = chains[chainNo].len;
              if (chainLen > 1){
                int fr = chains[chainNo].pixels[1].r;
                int fc = chains[chainNo].pixels[1].c;

                int dr = abs(fr-map->segments[noSegments].pixels[noSegmentPixels-1].r);
                int dc = abs(fc-map->segments[noSegments].pixels[noSegmentPixels-1].c);

                if (dr<=1 && dc<=1){startIndex = 1;}
              } //end-if
#endif
              /* Start a new chain & copy pixels from the new chain */
              for (int l=startIndex; l<chains[chainNo].len; l++){
                map->segments[noSegments].pixels[noSegmentPixels++] = chains[chainNo].pixels[l];
              } //end-for

              chains[chainNo].len = 0;  // Mark as copied
            } //end-for

            map->segments[noSegments].noPixels = noSegmentPixels;
            totalPixels += noSegmentPixels;        

            noSegments++;
          } //end-if          
        } //end-for

      } //end-else

  } //end-for-outer

  map->noSegments = noSegments;

  delete chainNos;
} //end-JoinAnchorPointsUsingSortedAnchors

///----------------------------------------------------------------------------------------------
/// Splits segments into multiple segments if they contain closed objects in them
/// Returns the new number of segments in the edgemap
/// 
int SplitSegments(EdgeMap *map){
  int noSegments = map->noSegments;

  for (int i=0; i<map->noSegments; i++){
    int sr = map->segments[i].pixels[0].r;
    int sc = map->segments[i].pixels[0].c;

    int er = map->segments[i].pixels[map->segments[i].noPixels-1].r;
    int ec = map->segments[i].pixels[map->segments[i].noPixels-1].c;

    bool withinNewSegment = false;
    int len = 3;
    int j = 3;
    int noPixels = map->segments[i].noPixels;
    while (j<noPixels-3){
      int r = map->segments[i].pixels[j].r;
      int c = map->segments[i].pixels[j].c;    
      len++;
      j++;

      int dr1 = abs(r-sr);
      int dc1 = abs(c-sc);

      int dr2 = abs(r-er);
      int dc2 = abs(c-ec);
      if ((dr1<=1 &&dc1<=1) || (dr2<=1 && dc2<=1)){
        if (withinNewSegment == false) map->segments[i].noPixels = len;
        else {map->segments[noSegments].noPixels = len; noSegments++;}

        map->segments[noSegments].pixels = &map->segments[i].pixels[j];

        withinNewSegment = true;
        len = 3;
        j += 3;
      } //end-if
    } //end-while

    if (j < noPixels) len += noPixels-j;
    else              len -= j-noPixels;

    if (len > 1 && withinNewSegment){map->segments[noSegments].noPixels = len; noSegments++;}
  } //end-for

  return noSegments;
} //end-SplitChains

///----------------------------------------------------------------------------------------------
/// Detect edges by edge drawing method
/// 
void DetectEdgesByEdgeDrawing(EDImage *ed){
  TimerClass timer;

  ed->noImages++;

  int width = ed->width;
  int height = ed->height;    

  unsigned char *image = ed->image;
  unsigned char *lowPassImg = ed->lowPassImg;

  /*------------ GAUSSIAN SMOOTHING -------------------*/
  timer.StartTiming();

  // Do Gauss low pass filter
#define CV_GAUSS_FILTER
#ifdef CV_GAUSS_FILTER
  /// OpenCV requires the image width to be a multiple of 4.
  /// If that is not the case, we have to handle it differently
  if (width % 4 == 0){
    IplImage *img1 = cvCreateImageHeader(cvSize(width, height),IPL_DEPTH_8U, 1);
    IplImage *img2 = cvCreateImageHeader(cvSize(width, height),IPL_DEPTH_8U, 1);
    img1->imageData = (char *)image;
    img2->imageData = (char *)lowPassImg;

    cvSmooth(img1, img2, CV_GAUSSIAN, 5);   // Gauss filter with kernel size 5

    cvReleaseImageHeader(&img1);
    cvReleaseImageHeader(&img2);

  } else {
    int w = width+4-width%4;

    IplImage *img1 = cvCreateImage(cvSize(w, height),IPL_DEPTH_8U, 1);
    IplImage *img2 = cvCreateImage(cvSize(w, height),IPL_DEPTH_8U, 1);

    char *p = img1->imageData;
    for (int i=0; i<height; i++){
      memcpy(p, &image[i*width], width);
      for (int j=width; j<w; j++) p[j]=p[width-1];
      p += w;
    } // end-for

    cvSmooth(img1, img2, CV_GAUSSIAN, 5);   // Gauss filter with kernel size 5

    p = img2->imageData;
    for (int i=0; i<height; i++){
      memcpy(&lowPassImg[i*width], p, width);
      p += w;
    } // end-for

    cvReleaseImage(&img1);
    cvReleaseImage(&img2);
  } //end-else
#else
  GaussFilter(image, lowPassImg, width, height);
#endif

  // Check parameters for sanity
  if (ed->GRADIENT_THRESH < 1) ed->GRADIENT_THRESH = 1;
  if (ed->ANCHOR_THRESH < 1) ed->ANCHOR_THRESH = 1;
  if (ed->SCAN_INTERVAL < 1) ed->SCAN_INTERVAL = 1;
  if (ed->MIN_PATH_LEN < 5) ed->MIN_PATH_LEN = 5;

  timer.StopTiming();
  ed->gaussTime += timer.SecondsBetween()*1e3;

  /*------------ COMPUTE GRADIENT & EDGE DIRECTION MAPS -------------------*/
  timer.StartTiming();

  switch (ed->GRADIENT_OPERATOR){
    case LSD_OPERATOR:     LSDOperator(ed);     break;
    case PREWITT_OPERATOR: PrewittOperator(ed); break;
    case SOBEL_OPERATOR:   SobelOperator(ed);   break;
    case SCHARR_OPERATOR:  ScharrOperator(ed);  break;
  } //end-switch

  timer.StopTiming();
  ed->gradientTime += timer.SecondsBetween()*1e3;

  /*------------ COMPUTE ANCHORS -------------------*/
  timer.StartTiming();

  ComputeAnchorPoints(ed);

  timer.StopTiming();
  ed->anchorTime += timer.SecondsBetween()*1e3;

  /*------------ JOIN ANCHORS -------------------*/
  timer.StartTiming();
    
  if (ed->useSortedAnchors) 
    JoinAnchorPointsUsingSortedAnchors(ed);
  else 
    JoinAnchorPoints(ed);

  timer.StopTiming();
  ed->joinTime += timer.SecondsBetween()*1e3;


  /*------------ VALIDATE SEGMENTS -------------------*/
  timer.StartTiming();
  if (ed->validateSegments){
    ValidateSegments(ed);    // Compute minimum over all pixels, but use chainLen/2 for NFA computation (suggested by DMM) 
    ExtractNewSegments(ed);  // Extract the new valid segments
  } //end-if

  timer.StopTiming();
  ed->validationTime += timer.SecondsBetween()*1e3;
} //DetectEdgesByEdgeDrawing

///================================ BELOW CODE IS USED FOR SEGMENT VALIDATION ===================
///---------------------------------------------------------------------------
/// Prewitt gradient map computation during segment validation
///
static void ComputePrewitt(EDImage *ed, double *H){
  int width = ed->width;
  int height = ed->height;

  // Compute gradients & gradient statistics
  int *grads = new int[ed->maxGradValue];
  memset(grads, 0, sizeof(int)*ed->maxGradValue);

  for (int i=1; i<height-1; i++){
    for (int j=1; j<width-1; j++){
      // Prewitt Operator in horizontal and vertical direction
      // A B C
      // D x E
      // F G H
      // gx = (C-A) + (E-D) + (H-F)
      // gy = (F-A) + (G-B) + (H-C)
      //
      // To make this faster: 
      // com1 = (H-A)
      // com2 = (C-F)
      // Then: gx = com1 + com2 + (E-D) = (H-A) + (C-F) + (E-D) = (C-A) + (E-D) + (H-F)
      //       gy = com1 - com2 + (G-B) = (H-A) - (C-F) + (G-B) = (F-A) + (G-B) + (H-C)
      // 
      int com1 = ed->image[(i+1)*width+j+1] - ed->image[(i-1)*width+j-1];
      int com2 = ed->image[(i-1)*width+j+1] - ed->image[(i+1)*width+j-1];

      int gx = abs(com1 + com2 + (ed->image[i*width+j+1] - ed->image[i*width+j-1]));
      int gy = abs(com1 - com2 + (ed->image[(i+1)*width+j] - ed->image[(i-1)*width+j]));
      int g = gx + gy;      

      ed->gradImg[i*width+j] = g;
      grads[g]++;
    } // end-for
  } //end-for

  // Compute probability function H
  int size = (width-2)*(height-2);
//  size -= grads[0];
  for (int i=ed->maxGradValue-1; i>0; i--) grads[i-1] += grads[i];
  for (int i=0; i<ed->maxGradValue; i++) H[i] = (double)grads[i]/((double)size);

  delete grads;
} //end-ComputePrewitt

///---------------------------------------------------------------------------
/// Number of false alarms code as suggested by Desolneux, Moisan and Morel (DMM)
///
#define EPSILON 1.0
static double NFA(int np, double prob, int len){
  double nfa = np;
  for (int i=0; i<len && nfa > EPSILON; i++) nfa*=prob;
//  for (int i=0; i<len; i++) nfa*=prob;

  return nfa;
} //end-NFA

///----------------------------------------------------------------------------------
/// Resursive validation using half of the pixels as suggested by DMM algorithm
/// We take pixels at Nyquist distance, i.e., 2 (as suggested by DMM)
///
static void TestSegment(EDImage *ed, int i, int index1, int index2, int np, double *H){
  int width = ed->width;

  int chainLen = index2-index1+1;
  if (chainLen < 5) return;

  /// Test from index1 to index2. If OK, then we are done. Otherwise, split into two and 
  /// recursively test the left & right halves

  // First find the min. gradient along the segment
  EdgeMap *map = ed->map;
  int minGrad = ed->maxGradValue;
  int minGradIndex;
  for (int k=index1; k<=index2; k++){
    int r = map->segments[i].pixels[k].r;
    int c = map->segments[i].pixels[k].c;
    if (ed->gradImg[r*width+c] < minGrad){minGrad = ed->gradImg[r*width+c]; minGradIndex = k;}
  } //end-for

  // Compute nfa
  double nfa = NFA(np, H[minGrad], chainLen/2);

  if (nfa <= EPSILON){
    for (int k=index1; k<=index2; k++){
      int r = map->segments[i].pixels[k].r;
      int c = map->segments[i].pixels[k].c;

      ed->edgeImg[r*width+c] = 255;
    } //end-for

    return;
  } //end-if  

  // Split into two halves. We divide at the point where the gradient is the minimum
  int end = minGradIndex-1;
  while (end > index1){
    int r = map->segments[i].pixels[end].r;
    int c = map->segments[i].pixels[end].c;

    if (ed->gradImg[r*width+c] <= minGrad) end--;
    else break;
  } //end-while

  int start = minGradIndex+1;
  while (start < index2){
    int r = map->segments[i].pixels[start].r;
    int c = map->segments[i].pixels[start].c;

    if (ed->gradImg[r*width+c] <= minGrad) start++;
    else break;
  } //end-while

  TestSegment(ed, i, index1, end, np, H);
  TestSegment(ed, i, start, index2, np, H);
} //end-TestSegment

static void ValidateSegments(EDImage *ed){  
  EdgeMap *map = ed->map;
  int width = ed->width;
  int height = ed->height;

  memset(ed->edgeImg, 0, width*height);

  double *H = new double[ed->maxGradValue];  
  ComputePrewitt(ed, H);
  
  // Compute np: # of segment pieces
  int np = 0;
  for (int i=0; i<map->noSegments; i++){
    int len = map->segments[i].noPixels;
    np += (len*(len-1))/2;
  } //end-for

  // Validate segments
  for (int i=0; i<map->noSegments; i++){
    TestSegment(ed, i, 0, map->segments[i].noPixels-1, np, H);
  } //end-for

  delete H;
} //end-ValidateSegments

///----------------------------------------------------------------------------------------------
/// After the validation of the edge segments, extracts the valid ones
/// In other words, updates the valid segments' pixel arrays and their lengths
/// 
static void ExtractNewSegments(EDImage *ed){
  /// Extract new segments
  int width = ed->width;
  EdgeMap *map = ed->map;
  EdgeSegment *segments = &map->segments[ed->map->noSegments];
  int noSegments = 0;

  for (int i=0; i<map->noSegments; i++){
    int start = 0;
    while (start < map->segments[i].noPixels){

      while (start < map->segments[i].noPixels){
        int r = map->segments[i].pixels[start].r;
        int c = map->segments[i].pixels[start].c;

        if (ed->edgeImg[r*width+c]) break;
        start++;
      } //end-while
  
      int end = start+1;
      while (end < map->segments[i].noPixels){
        int r = map->segments[i].pixels[end].r;
        int c = map->segments[i].pixels[end].c;

        if (ed->edgeImg[r*width+c] == 0) break;
        end++;
      } //end-while

      int len = end-start;
      if (len >= 8){
        // A new segment. Accepted only only long enough (whatever that means)
        segments[noSegments].pixels = &map->segments[i].pixels[start];
        segments[noSegments].noPixels = len;
        noSegments++;
      } //end-else

      start = end+1;
    } //end-while
  } //end-for

  // Copy to ed
  for (int i=0; i<noSegments; i++) map->segments[i] = segments[i];
  map->noSegments = noSegments;
} //end-ExractNewSegments

///--------------------------------------------------------------
/// Detect Edges by Edge Drawing & validate them by DMM algorithm
/// Can use any of the well-known operators
///
void DetectEdgesByEDParameterFree(EDImage *ed){
  TimerClass timer;

  if (ed->GRADIENT_OPERATOR == LSD_OPERATOR){
    ed->GRADIENT_OPERATOR = LSD_OPERATOR;
    ed->GRADIENT_THRESH = 7;
    ed->ANCHOR_THRESH = 1;
    ed->SCAN_INTERVAL = 1;
    ed->MIN_PATH_LEN = 10;

  } else if (ed->GRADIENT_OPERATOR == SOBEL_OPERATOR){
    ed->GRADIENT_OPERATOR = SOBEL_OPERATOR;
    ed->GRADIENT_THRESH = 17;
    ed->ANCHOR_THRESH = 1;
    ed->SCAN_INTERVAL = 1;
    ed->MIN_PATH_LEN = 10;

  } else if (ed->GRADIENT_OPERATOR == SCHARR_OPERATOR){
    ed->GRADIENT_OPERATOR = SCHARR_OPERATOR;
    ed->GRADIENT_THRESH = 68;
    ed->ANCHOR_THRESH = 1;
    ed->SCAN_INTERVAL = 1;
    ed->MIN_PATH_LEN = 10;

//  } else if (ed->GRADIENT_OPERATOR == PREWITT_OPERATOR){
  } else {
    // By default we use the Prewitt operator for parameter-free ED
    ed->GRADIENT_OPERATOR = PREWITT_OPERATOR;
    ed->GRADIENT_THRESH = 13;  
    ed->ANCHOR_THRESH = 1;     
    ed->SCAN_INTERVAL = 1;
    ed->MIN_PATH_LEN = 10;
  } //end-else

  ed->useSortedAnchors = true;
  ed->validateSegments = true;
  DetectEdgesByEdgeDrawing(ed);
} //end-DetectEdgesByEDParameterFree

///-------------------------------------------------------------
/// Dump edge segments to file: Just the pixels, not the chains
///
void DumpEdgeSegments2File(EdgeMap *map, char *fname){
  FILE *fp = fopen(fname, "w");

  for (int i=0; i<map->noSegments; i++){
    fprintf(fp, "\n+===SEGMENT: %3d. NoPixels: %d ====+\n", i, map->segments[i].noPixels);

    fprintf(fp, "+------+------+-------+\n");
    fprintf(fp, "|  Row |Column|Segment|\n"); 
    fprintf(fp, "+------+------+-------+\n");

    for (int j=0; j<map->segments[i].noPixels; j++){      
      fprintf(fp, "|%6d|%6d|%7d|\n", map->segments[i].pixels[j].r,  map->segments[i].pixels[j].c, i);     
    } //end-for

    fprintf(fp, "+------+------+-------+\n");
  } //end-for

  fclose(fp);
} //end-DumpEdgeSegments2File2

void SaveFrameRAW(char *sFileName, char *m_buffer, int width, int height){
  FILE *fp = fopen( sFileName, "wb" );
  fwrite(m_buffer, 1, width*height, fp);

  fclose( fp );
} //end-SaveFrameRAW

//dll entry
//extern int EdgeDrawingDetectorV(unsigned char *srcImg,  unsigned char *edgeImg, int width, int height, int DIFF_THRESH, 
//						   int SOBEL_THRESH, int MIN_PATH_LEN, int DETAIL_RATIO, myEdgeMap *myMap)
//{
//	double total = 0.0;
//	TimerClass timer;
//	timer.StartTiming();
//
//	int angles[3][3] = {{135, 90, 45}, {180, -1, 0}, {225, 270, 315}};  
//	
//	EDImage *ed = new EDImage(width, height);
//	ed->image = srcImg;
//	//ed->computeChains = true;
//
//	ed->ANCHOR_THRESH = DIFF_THRESH;
//	ed->GRADIENT_THRESH = SOBEL_THRESH;
//	ed->MIN_PATH_LEN = MIN_PATH_LEN;
//	ed->SCAN_INTERVAL = DETAIL_RATIO;
//
//	DetectEdgesByEdgeDrawing(ed);
//
//	//SaveFrameRAW("EdgeDrawing.raw", (char *)ed->edgeImg, width, height);
//
//	int pixCount = 0;
//	myMap->noOfSegments = ed->map->noSegments;
//
//	for(int i=0 ; i < ed->map->noSegments ; i++)
//	{
//		myMap->segments[i].noOfPixels = ed->map->segments[i].noPixels; 
//		myMap->segments[i].pixels = (myPixel*) myMap->pixels + pixCount;
//		
//		for(int j=0 ; j < ed->map->segments[i].noPixels ; j++)
//		{
//			int x = ed->map->segments[i].pixels[j].c;
//			int y = ed->map->segments[i].pixels[j].r;
//
//			myMap->segments[i].pixels[j].c = x; //ed->map->segments[i].pixels[j].c;
//			myMap->segments[i].pixels[j].r = y; //ed->map->segments[i].pixels[j].r;
//
//#define ED_CONTOURS			
//#ifdef ED_CONTOURS
//
//			int index = y*width+x; //x*width+y;
//			//with fast arctan
//			//myMap->segments[i].pixels[j].angle = (int)myAtan2(ed->gYImage[index], ed->gXImage[index]);
//
//			//with atan function
//			double ratioXY = (double)ed->gXImage[index] / ed->gYImage[index];
//			myMap->segments[i].pixels[j].angle = (int)(atan(ratioXY) * 180.0 / 3.14159);
//
//#else
//			if(j < ed->map->segments[i].noPixels - 1)
//			{
//				int x = ed->map->segments[i].pixels[j].c;		//pixels[k].c;
//				int y = ed->map->segments[i].pixels[j].r;		//pixels[k].r;
//				int nextX = ed->map->segments[i].pixels[j+1].c;	//pixels[k+1].c;
//				int nextY = ed->map->segments[i].pixels[j+1].r; //pixels[k+1].r;    
//			    
//				int yDiff = y-nextY;
//				int xDiff = x-nextX;
//				if (abs(yDiff) >= 2 || abs(xDiff) >= 2) continue;  // not neigbors. skip
//				
//				myMap->segments[i].pixels[j].angle = angles[yDiff+1][xDiff+1];
//			}
//#endif
//			pixCount++;
//		}
//		
//	}
//
//	timer.StopTiming();
//	total += timer.SecondsBetween()*1e3;
//
//	delete ed;
//	return (int)total;
//}