
// Copyright (c) 2017-2024 King Abdullah University of Science and Technology,
// All rights reserved.
// ExaGeoStat is a software package, provided by King Abdullah University of Science and Technology (KAUST).

/**
 * @file JDFKernels.cu
 * @brief Implementation of DCMG kernel used in dmatrix generation.
 * @version 2.0.0
 * @author Mahmoud ElKarargy
 * @author Sameh Abdulah
 * @date 2025-02-04
**/

#include <cuda_runtime.h>
#include <cublas_v2.h>
#include <runtime/parsec/JDFKernels.h>

#define CHUNKSIZE 32

extern "C" {
/**
 * @brief 3D kernel to compute power-exponential covariance distance.
 * @param[in,out] apDescA        (double*) The output matrix on the GPU (already allocated).
 * @param[in]     aM, aN         Dimensions of the tile/block we are processing.
 * @param[in]     aM0, aN0       Global offsets (not always used, but kept for consistency).
 * @param[in]     apLocationX1, apLocationY1 Arrays of x, y for first location set (size aN).
 * @param[in]     apLocationX2, apLocationY2 Arrays of x, y for second location set (size aM).
 * @param[in]     aLocalTheta0, aLocalTheta1, aLocalTheta2    Kernel hyperparameters (e.g., sigma², correlation length, exponent).
 * @param[in]     aDistanceMetric 1 if you want some specialized distance formula (like "great circle"),
 *                                0 for Euclidian. Example below just does Euclidian in 3D.
 */
__global__ void DcmgArrayKernel(double *apDescA, int aM, int aN, int aM0, int aN0, double* apLocationX1, double* apLocationY1,
                                double* apLocationX2, double* apLocationY2, double aLocalTheta0, double aLocalTheta1,
                                double aLocalTheta2, int aDistanceMetric)
    {
    const int tx  = threadIdx.x;
    const int ty  = threadIdx.y;
    const int idx = blockIdx.x * blockDim.x + tx;
    const int idy = blockIdx.y * blockDim.y + ty;

    if(idx>=aM || idy >=aN){return;}

    //double x0, y0;
    double expr  = 0.0;
    double expr1 = 0.0;

    double sigma_square = aLocalTheta0;// * localtheta[0];

    expr = sqrt(pow((apLocationX2[idx] - apLocationX1[idy]), 2) +
            pow((apLocationY2[idx] - apLocationY1[idy]), 2));

    expr1 = pow(expr, aLocalTheta2);
    if(expr == 0)
        apDescA[idx + idy * aM] = sigma_square /*+ 1e-4*/;
    else
        apDescA[idx + idy * aM] = sigma_square *  exp(-(expr1/aLocalTheta1)); // power-exp kernel
}

void DcmgArray(double *apDescA, int aM, int aN, int aM0, int aN0, double* apLocationX1, double* apLocationY1,
               double* apLocationX2, double* apLocationY2, double *apLocalTheta, int aDistanceMetric, cudaStream_t aCudaStream) {

    int block_num_x= (aM+CHUNKSIZE-1)/CHUNKSIZE;
    int block_num_y= (aN+CHUNKSIZE-1)/CHUNKSIZE;
    dim3 dimBlock(CHUNKSIZE,CHUNKSIZE);
    dim3 dimGrid(block_num_x,block_num_y);

    DcmgArrayKernel<<<dimGrid,dimBlock, 0, aCudaStream>>>(apDescA, aM, aN, aM0, aN0, apLocationX1, apLocationY1,
               apLocationX2, apLocationY2,  apLocalTheta[0],apLocalTheta[1],apLocalTheta[2], aDistanceMetric);

}

/**
 * @brief 3D kernel to compute power-exponential covariance distance.
 * @param[in,out] apDescA        (double*) The output matrix on the GPU (already allocated).
 * @param[in]     aM, aN         Dimensions of the tile/block we are processing.
 * @param[in]     aM0, aN0       Global offsets (not always used, but kept for consistency).
 * @param[in]     apLocationX1, apLocationY1, apLocationZ1    Arrays of x, y, z for first location set (size aN).
 * @param[in]     apLocationX2, apLocationY2, apLocationZ2    Arrays of x, y, z for second location set (size aM).
 * @param[in]     aLocalTheta0, aLocalTheta1, aLocalTheta2    Kernel hyperparameters (e.g., sigma², correlation length, exponent).
 * @param[in]     aDistanceMetric 1 if you want some specialized distance formula (like "great circle"),
 *                                0 for Euclidian. Example below just does Euclidian in 3D.
 */
__global__ void DcmgArray3DKernel(double *apDescA, int aM, int aN, int aM0, int aN0, const double* apLocationX1,
                                  const double* apLocationY1, const double* apLocationZ1, const double* apLocationX2,
                                  const double* apLocationY2, const double* apLocationZ2, double aLocalTheta0,
                                  double aLocalTheta1, double aLocalTheta2, int aDistanceMetric)
{
    int tx  = threadIdx.x;
    int ty  = threadIdx.y;
    int idx = blockIdx.x * blockDim.x + tx;
    int idy = blockIdx.y * blockDim.y + ty;

    if (idx >= aM || idy >= aN) {
        return;
    }

    // 1. Compute 3D distance between point idx (in X2/Y2/Z2) and point idy (in X1/Y1/Z1).
    double dx = apLocationX2[idx] - apLocationX1[idy];
    double dy = apLocationY2[idx] - apLocationY1[idy];
    double dz = apLocationZ2[idx] - apLocationZ1[idy];

    // For Euclidian distance
    double expr  = sqrt(dx * dx + dy * dy + dz * dz);

    // 2. Apply the power-exponential formula
    double sigma_sq = aLocalTheta0;  // Could be apLocalTheta[0]
    double expr1    = pow(expr, aLocalTheta2);

    // If distance is zero, put the variance on the diagonal
    if (expr == 0.0) {
        apDescA[idx + idy * aM] = sigma_sq; // + a small jitter if desired
    } else {
        apDescA[idx + idy * aM]
            = sigma_sq * exp( -(expr1 / aLocalTheta1) );
    }
}

void DcmgArray3D(double *apDescA, int aM, int aN, int aM0, int aN0, double* apLocationX1, double* apLocationY1,
                 double* apLocationZ1,  double* apLocationX2, double* apLocationY2, double* apLocationZ2, double *apLocalTheta,
                 int aDistanceMetric, cudaStream_t aCudaStream){

    int block_num_x = (aM + CHUNKSIZE - 1) / CHUNKSIZE;
    int block_num_y = (aN + CHUNKSIZE - 1) / CHUNKSIZE;

    dim3 dimBlock(CHUNKSIZE, CHUNKSIZE);
    dim3 dimGrid(block_num_x, block_num_y);

    // Launch the 3D kernel
    DcmgArray3DKernel<<<dimGrid, dimBlock, 0, aCudaStream>>>(
        apDescA,
        aM, aN, aM0, aN0,
        apLocationX1, apLocationY1, apLocationZ1,
        apLocationX2, apLocationY2, apLocationZ2,
        apLocalTheta[0], apLocalTheta[1], apLocalTheta[2],
        aDistanceMetric
    );
}
}