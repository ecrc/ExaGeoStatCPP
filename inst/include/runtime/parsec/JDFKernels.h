
// Copyright (c) 2017-2024 King Abdullah University of Science and Technology,
// All rights reserved.
// ExaGeoStat is a software package, provided by King Abdullah University of Science and Technology (KAUST).

/**
 * @file JDFKernels.h
 * @brief Declaration of DCMG kernel used in dmatrix generation.
 * @version 2.0.0
 * @author Mahmoud ElKarargy
 * @author Sameh Abdulah
 * @date 2025-02-04
**/

#ifdef __cplusplus
extern "C" {
#endif
    /**
     * @brief Computes a 2D power-exponential covariance matrix on the GPU.  This function launches a CUDA kernel to fill
     * the matrix apDescA with power-exponential covariance values computed from 2D locations.
     * @param[out] apDescA Pointer to the (device) array of size  aM *  aN where the results will be stored.
     * @param[in] aM  Number of rows (e.g., in the  x dimension).
     * @param[in] aN Number of columns (e.g., in the  y dimension).
     * @param[in] aM0 Global offset in the row dimension (often 0, used if the tile starts partway through the matrix).
     * @param[in] aN0 Global offset in the column dimension (similar to  aM0).
     * @param[in] apLocationX1 Pointer to the array of  x coordinates for the first location set (size =  aN).
     * @param[in] apLocationY1 Pointer to the array of  y coordinates for the first location set (size =  aN).
     * @param[in] apLocationX2 Pointer to the array of  x coordinates for the second location set (size =  aM).
     * @param[in] apLocationY2 Pointer to the array of  y coordinates for the second location set (size =  aM).
     * @param[in] apLocalTheta Array of hyperparameters for the power-exponential kernel.
     * @param[in] aDistanceMetric Indicator for which distance metric to use (e.g., 0 for standard Euclidian,1 for great-circle, etc.).
     * @param[in] aCudaStream CUDA stream on which to launch the kernel.
     * @return void
     */
void DcmgArray(double *apDescA, int aM, int aN, int aM0, int aN0, double* apLocationX1, double* apLocationY1,
               double* apLocationX2, double* apLocationY2, double *apLocalTheta, int aDistanceMetric, cudaStream_t aCudaStream);

    /**
     * @brief Computes a 3D power-exponential covariance matrix on the GPU. Similar to DcmgArray, but extended to 3D coordinates.
     * @param[out] apDescA Pointer to the (device) array of size aM * aN where the results will be stored.
     * @param[in] aM Number of rows (e.g.,  x dimension).
     * @param[in] aN Number of columns (e.g.,  y dimension).
     * @param[in] aM0 Global offset in the row dimension (tile offset).
     * @param[in] aN0 Global offset in the column dimension (tile offset).
     * @param[in] apLocationX1 Pointer to the array of  x coordinates for the first location set (size =  aN).
     * @param[in] apLocationY1 Pointer to the array of  y coordinates for the first location set (size =  aN).
     * @param[in] apLocationZ1 Pointer to the array of  z coordinates for the first location set (size =  aN).
     * @param[in] apLocationX2 Pointer to the array of  x coordinates for the second location set (size =  aM).
     * @param[in] apLocationY2 Pointer to the array of  y coordinates for the second location set (size =  aM).
     * @param[in] apLocationZ2 Pointer to the array of  z coordinates for the second location set (size =  aM).
     * @param[in] apLocalTheta Array of hyperparameters for the power-exponential kernel.
     * @param[in] aDistanceMetric Indicator for which distance metric to use (e.g., 0 for Euclidian).
     * @param[in] aCudaStream CUDA stream on which to launch the kernel.
     * @return void
     */
void DcmgArray3D(double *apDescA, int aM, int aN, int aM0, int aN0, double* apLocationX1, double* apLocationY1,
                 double* apLocationZ1,  double* apLocationX2, double* apLocationY2, double* apLocationZ2, double *apLocalTheta,
                 int aDistanceMetric, cudaStream_t aCudaStream);

#ifdef __cplusplus
}
#endif