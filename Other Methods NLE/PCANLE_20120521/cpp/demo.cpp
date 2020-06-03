#include "PCANoiseLevelEstimator.h"

#include <stdio.h>
#include <math.h>


//================================================================================

static bool
ReadP5PGM( const char* file_name, double** image_data, int* image_w, int* image_h )
{
    bool success = false;
    
    FILE* file = fopen( file_name, "rb" );

    if( file )
    {
        char P5_header[32] = {0};
        int max_value = 0;
        fscanf( file, "%s %i %i %i\n", P5_header, image_w, image_h, &max_value );
        
        unsigned pixel_count = (*image_w) * (*image_h);
        
        *image_data = new double[pixel_count];

        for( unsigned i = 0; i < pixel_count; ++i )
        {
            unsigned char pixel = 0;
            fread( &pixel, 1, 1, file );
            (*image_data)[i] = pixel;
        }
        
        fclose( file );
        success = true;
    }
    
    return success;
}


//================================================================================

int
main()
{
    const char* file_name = "../data/cameraman-std=5.pgm";
    
    double* image_data = 0;
    int image_w = 0;
    int image_h = 0;
    
    if( ReadP5PGM(file_name, &image_data, &image_w, &image_h) )
    {
        double var = EstimateNoiseVariance( image_data, image_w, image_h );
        double std = sqrt(var);
        
        printf( "expected noise standard deviation = 5.00\n" );
        printf( "computed noise standard deviation = %.2f\n", std );
        
        delete [] image_data;
    }
    else
    {
        printf( "failed to read input file \"%s\"\n", file_name );
    }
    
    return 0;
}
