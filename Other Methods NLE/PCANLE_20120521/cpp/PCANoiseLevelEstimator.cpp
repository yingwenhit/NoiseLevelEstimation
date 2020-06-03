#include "PCANoiseLevelEstimator.h"

#include "Param.h"
#include "Vector.h"
#include "EigenvalueComputation.h"

#include <math.h>
#include <stdio.h>
#include <algorithm>


//================================================================================

inline int
Round( double x )
{
    return x > 0.0 ? int(x + 0.5) : int(x - 0.5);
}


//================================================================================

inline int
Clamp( int value, int begin_value, int end_value )
{
    if( value < begin_value )
        return begin_value;
        
    if( value > end_value )
        return end_value;
        
    return value;
}


//================================================================================

struct BlockInfo
{
    double      variance;
    unsigned    offset;
};

inline bool
operator < ( const BlockInfo& a, const BlockInfo& b )
{
    return a.variance < b.variance;
}


//================================================================================

static unsigned
ComputeBlockInfo( const double* image_data, int image_w, int image_h, BlockInfo* block_info )
{
    unsigned block_count = 0;

    for( int y = 0; y < image_h - M2; ++y )
    {
        for( int x = 0; x < image_w - M1; ++x )
        {
            double sum1 = 0.0;
            double sum2 = 0.0;
            unsigned clipped_pixel_count = 0;
            
            for( int by = y; by < y + M2; ++by )
            {
                for( int bx = x; bx < x + M1; ++bx )
                {
                    double val = image_data[by*image_w + bx];
                    sum1 += val;
                    sum2 += val*val;

                    if( val == 0.0 || val == 255.0 )
                        ++clipped_pixel_count;
                }
            }

            if( clipped_pixel_count <= MaxClippedPixelCount )
            {
                BlockInfo bi;
                bi.variance = (sum2 - sum1*sum1/M) / M;
                bi.offset = y*image_w + x;
                
                block_info[block_count] = bi;
                ++block_count;
            }
        }
    }

    return block_count;
}


//================================================================================

static unsigned
ComputeStatistics( const double* image_data, int image_w, int image_h, const BlockInfo* block_info, unsigned block_count, Vector* sum1, Matrix* sum2, unsigned* subset_size )
{
    unsigned subset_count = 0;
    
    for( double p = 1.0; p - MinLevel > -1e-6; p -= LevelStep )
    {
        double      q           = p - LevelStep - MinLevel > -1e-6 ? p - LevelStep : 0.0;
        unsigned    max_index   = block_count - 1;
        unsigned    beg         = Clamp( Round(q*max_index), 0, max_index );
        unsigned    end         = Clamp( Round(p*max_index), 0, max_index );

        Vector curr_sum1 = {{0}};
        Matrix curr_sum2 = {{0}};
        
        for( unsigned k = beg; k < end; ++k )
        {
            unsigned offset = block_info[k].offset;

            double block[M];
            for( int by = 0; by < M2; ++by )
                for( int bx = 0; bx < M1; ++bx )
                    block[by*M1+bx] = image_data[offset + by*image_w + bx];
            
            for( unsigned i = 0; i < M; ++i )
                curr_sum1[i] += block[i];

            for( unsigned i = 0; i < M; ++i )
                for( unsigned j = i; j < M; ++j )
                    curr_sum2[i*M+j] += block[i]*block[j];
        }
        
        sum1[subset_count] = curr_sum1;
        sum2[subset_count] = curr_sum2;
        subset_size[subset_count] = end - beg;

        ++subset_count;
    }
    
    for( unsigned i = subset_count - 1; i > 0; --i )
    {
        sum1[i-1] += sum1[i];
        sum2[i-1] += sum2[i];
        subset_size[i-1] += subset_size[i];
    }

    return subset_count;
}


//================================================================================

static double
ComputeUpperBound( const BlockInfo* block_info, unsigned block_count )
{
    unsigned max_index  = block_count - 1;
    unsigned index      = Clamp( Round(UpperBoundLevel*max_index), 0, max_index );
    return UpperBoundFactor * block_info[index].variance;
}


//================================================================================

static Vector
ApplyPCA( const Vector& sum1, const Matrix& sum2, unsigned subset_size )
{
    double mean[M];
    for( unsigned i = 0; i < M; ++i )
        mean[i] = sum1[i] / subset_size;

    Matrix cov_matrix;
    for( unsigned i = 0; i < M; ++i )
    {
        for( unsigned j = i; j < M; ++j )
        {
            cov_matrix[i*M+j] = sum2[i*M+j]/subset_size - mean[i]*mean[j];
            cov_matrix[j*M+i] = cov_matrix[i*M+j];
        }
    }
    
    return ComputeEigenvalues( &cov_matrix );
}


//================================================================================

static double
GetNextEstimate( const Vector* sum1, const Matrix* sum2, unsigned* subset_size, unsigned sum_count, double prev_estimate, double upper_bound )
{
    double var = 0.0;
    
    for( unsigned i = 0; i < sum_count; ++i )
    {
        Vector eigen_value = ApplyPCA( sum1[i], sum2[i], subset_size[i] );

        var = eigen_value[0];

        if( var < 1e-6 )
            break;

        double diff             = eigen_value[EigenValueCount-1] - eigen_value[0];
        double diff_threshold   = EigenValueDiffThreshold * prev_estimate / sqrt((double)subset_size[i]);
        
        if( diff < diff_threshold && var < upper_bound )
            break;
    }
    
    return var < upper_bound ? var : upper_bound;
}


//================================================================================

double
EstimateNoiseVariance( const double* image_data, int image_w, int image_h )
{
    BlockInfo*  block_info  = new BlockInfo[image_w*image_h];
    unsigned    block_count = ComputeBlockInfo( image_data, image_w, image_h, block_info );
    
    std::sort( block_info, block_info + block_count );
    
    Vector sum1[MaxSubsetCount];
    Matrix sum2[MaxSubsetCount];
    unsigned subset_size[MaxSubsetCount];
    unsigned subset_count = ComputeStatistics( image_data, image_w, image_h, block_info, block_count, sum1, sum2, subset_size );
    
    double upper_bound  = ComputeUpperBound( block_info, block_count );
    double prev_var     = 0.0;
    double var          = upper_bound;

    for( unsigned iter = 0; iter < 10; ++iter )
    {
        if( fabs(prev_var - var) < 1e-6 )
            break;
        
        prev_var = var;
        var = GetNextEstimate( sum1, sum2, subset_size, subset_count, var, upper_bound );
    }
    
    delete [] block_info;

    return var > 0.0 ? var : 0.0;
}

