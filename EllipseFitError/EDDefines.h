#ifndef _ED_DEFINES_H_
#define _ED_DEFINES_H_

enum GradientOperator {LSD_OPERATOR=101, PREWITT_OPERATOR=102, SOBEL_OPERATOR=103, SCHARR_OPERATOR=104};

struct Pixel {int r, c;};

struct EdgeSegment {
	Pixel *pixels;
    int noPixels;
};

struct EdgeMap {
public:
  EdgeSegment *segments;
  int noSegments;

public:
  // constructor
  EdgeMap(int capacity){
    segments = new EdgeSegment[capacity];
    noSegments = 0;
  } //end-EdgeMap

  // Destructor
  ~EdgeMap(){delete segments;}
};

#endif