#ifndef PCANOISELEVELESTIMATOR_H_
#define PCANOISELEVELESTIMATOR_H_

//================================================================================

// Noise variance estimation.
// Arguments:
//   1. image_data - pointer to the memory block, which contains the input image.
//                   The image is assumed to be stored in raster order,
//                   i.e. pixel (x,y) is stored at image_data[x + image_w*y]
//   2. image_w    - the image width
//   3. image_h    - the image height
// Return value: the noise variance.
double EstimateNoiseVariance( const double* image_data, int image_w, int image_h );


#endif /*PCANOISELEVELESTIMATOR_H_*/
