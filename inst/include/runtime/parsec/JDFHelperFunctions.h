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
