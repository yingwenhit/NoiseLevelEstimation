CMAKE_MINIMUM_REQUIRED( VERSION 2.4 )
PROJECT( PCANLE_demo )

SET( CMAKE_BUILD_TYPE Release )

ADD_EXECUTABLE( PCANLE_demo
	demo.cpp
	EigenvalueComputation.cpp
	PCANoiseLevelEstimator.cpp
)
