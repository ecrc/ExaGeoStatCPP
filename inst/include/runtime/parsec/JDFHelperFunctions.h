/**
 * @file JDFHelperFunctions.h
 * @brief A header file for declarations of JDF helper functions.
 * @details Contains function prototypes for JDF operations, including computations,
 *          transformations, file I/O, and MPI-related data handling for matrix and data structures.
 * @version 2.0.0
 * @author Mahmoud ElKarargy
 * @date 2024-10-20
**/

#include <runtime/parsec/ParsecHeader.h>
#include <complex.h>

#ifdef USE_CUDA
#include <runtime/parsec/GPUHelperFunctions.h>
#endif

/**
 * @brief Calculates a unique single index from given dimensions.
 * @param[in] aN The row index.
 * @param[in] aM The column index.
 * @return The calculated single index.
 *
 */
int CalculateSingleIndex(int aN, int aM);

/**
 * @brief Sums the elements of a double-precision data matrix.
 * @param[in] apData Pointer to the data matrix.
 * @param[in] aColumn The number of columns in the matrix.
 * @param[in] aRow The number of rows in the matrix.
 * @return The sum of matrix elements as a double.
 *
 */
double SumDoubleData(double *apData, int aColumn, int aRow);

/**
 * @brief Sums the elements of a complex double-precision data matrix.
 * @param[in] apData Pointer to the complex data matrix.
 * @param[in] aColumn The number of columns in the matrix.
 * @param[in] aRow The number of rows in the matrix.
 * @return The sum of matrix elements as a complex double.
 * @return void
 *
 */
complex double SumComplexData(complex double *apData, int aColumn, int aRow);

/**
 * @brief Performs forward Spherical Harmonic Transform (SHT) calculations.
 * @param[in,out] apFlm Pointer to SHT coefficients.
 * @param[in] apF_data Pointer to spatial data for transformation.
 * @param[in] aFDataM Number of rows in spatial data.
 * @param[in] aFDataN Number of columns in spatial data.
 * @param[in] apEt1, apEt2, apEp, apSlmn, apIe, apIo, apP Arrays for intermediate calculations.
 * @param[in] aEt1M, aEt2M, aEpM, aEpN, aSlmnN, aSlmnM, aIeM, aIeN, aIoM, aIoN Dimensions of respective arrays.
 * @param[in] apD, apGmtheta_r, apFmnm, apTmp1, apTmp2 Additional intermediate matrices.
 * @param[in] aL Maximum degree for the SHT.
 * @return void
 *
 */
void ForwardSHTHelper(double *apFlm, complex double *apF_data, int aFDataM, int aFDataN,
                      complex double *apEt1, int aEt1M, complex double *apEt2, int aEt2M,
                      complex double *apEp, int aEpM, int aEpN, complex double *apSlmn,
                      int aSlmnN, int aSlmnM, complex double *apIe, int aIeM, int aIeN,
                      complex double *apIo, int aIoM, int aIoN, complex double *apP,
                      int aPM, int aPN, complex double *apD, complex double *apGmtheta_r,
                      complex double *apFmnm, complex double *apTmp1,
                      complex double *apTmp2, int aL);

/**
 * @brief Performs inverse Spherical Harmonic Transform (SHT) calculations.
 * @param[in] apFlm Pointer to coefficients from forward SHT.
 * @param[out] apF_spatial Pointer to result spatial data.
 * @param[in] apZlm, apSC, apSmt Arrays for intermediate calculations.
 * @param[in] aL Maximum degree for the SHT.
 * @return void
 *
 */
void InverseSHTHelper(double *apFlm, double *apF_spatial, double *apZlm, double *apSC,
                      double *apSmt, int aL);

/**
 * @brief Computes the sum of the elements within a matrix block.
 * @details This function iterates over a specified matrix block and accumulates
 *          the values of all elements into a single double-precision result.
 *          It also checks for NaN entries (not-a-number) and sets an indicator
 *          if any are encountered. The final sum is returned as a double.
 * @param[in] apA Pointer to the beginning of the matrix block data.
 * @param[in] aMb The number of rows in the matrix block.
 * @param[in] aNb The number of columns in the matrix block.
 * @param[in] aLda The leading dimension of the matrix.
 * @return The double-precision sum of all elements in the specified matrix block.
 *
 */
double ParsecMatrixSumCore(double *apA, int aMb, int aNb, int aLda);

#ifdef USE_CUDA

/**
 * @brief Performs forward Spherical Harmonic Transform (SHT) calculations on the GPU.
 * This function is used to compute the forward SHT using GPU acceleration.
 * @param[in] apFlm Pointer to coefficients used in the forward SHT.
 * @param[in,out] apF_data Pointer to GPU memory holding data for the SHT calculations.
 * @param[in] apEt1, apEt2, apEp Arrays for intermediate exponential terms in the SHT.
 * @param[in] apSlmn Pointer to an array for storing spherical harmonics coefficients.
 * @param[in] apIe, apIo Pointers to arrays for storing intermediate SHT results.
 * @param[in] apP Pointer to an array of precomputed SHT coefficients.
 * @param[in] apD Pointer to the differential operator array for SHT.
 * @param[in] apGmtheta_r Pointer to the array of spherical harmonics multipliers.
 * @param[in] apFmnm Pointer to an array for frequency modulation factors.
 * @param[in,out] apTmp1, apTmp2 Temporary buffers for intermediate GPU calculations.
 * @param[in] aL Maximum degree for the SHT.
 * @param[in] apCudaDevice Pointer to the CUDA device module used for the SHT.
 * @param[in] apCudaStream Pointer to the CUDA execution stream for asynchronous operations.
 * @param[in] apWorkSpace Pointer to the gpu workspace.
 * @return void
 *
 */
void ForwardSHTGPUCore(double *apFlm, cuDoubleComplex *apF_data, cuDoubleComplex *apEt1,
                       cuDoubleComplex *apEt2, cuDoubleComplex *apEp, cuDoubleComplex *apSlmn,
                       cuDoubleComplex *apIe, cuDoubleComplex *apIo, cuDoubleComplex *apP,
                       cuDoubleComplex *apD, parsec_device_cuda_module_t *apCudaDevice,
                       parsec_gpu_task_t *apGpuTask, parsec_cuda_exec_stream_t *apCudaStream,
                       void* apWorkSpace, int aL);
#endif