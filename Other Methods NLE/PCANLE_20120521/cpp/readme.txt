Program for noise level estimation by principal component analysis.
C++ version.
Author: Stanislav Pyatykh

The algorithm is described in the following article:
S. Pyatykh, J. Hesser, and L. Zheng, "Image Noise Level Estimation by Principal Component Analysis", IEEE Transactions on Image Processing, Volume: 22, Issue: 2, Pages: 687 - 699, February 2013
DOI: 10.1109/TIP.2012.2221728

Contents:
CMakeLists.txt - project file used to build the program using CMake
demo.cpp - demonstration program (contains the entry point)
PCANoiseLevelEstimator.h/.cpp - the noise variance estimation algorithm
EigenvalueComputation.h/.cpp - eigenvalue computation algorithm
Param.h - noise variance estimation algorithm parameters
Vector.h - class for representation of vectors and matrices
