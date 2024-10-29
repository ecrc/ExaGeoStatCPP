// Copyright (c) 2017-2024 King Abdullah University of Science and Technology,
// All rights reserved.
// ExaGeoStat is a software package, provided by King Abdullah University of Science and Technology (KAUST).

/**
 * @file JobDescriptionFormat.h
 * @brief A header file for parsec generated functions from jdf files.
 * @version 2.0.0
 * @author Mahmoud ElKarargy
 * @date 2024-10-08
**/

#include <runtime/parsec/ParsecHeader.h>

/**
 * @brief Reads a CSV file into a matrix description.
 * @details This function reads data from a CSV file and populates the specified matrix description.
 * @param[in] apContext Pointer to the parsec context.
 * @param[in, out] apDesc Pointer to the matrix block cyclic descriptor.
 * @param[in] aMB Number of rows in the block.
 * @param[in] aNB Number of columns in the block.
 * @param[in] aNodes Number of nodes.
 * @param[in] aTimeSlot Time slot for the data.
 * @param[in] apFilename Filename of the CSV file.
 * @param[in] aRank Rank of the current process.
 * @param[in] aVerbose Verbosity level for output.
 * @param[in] aGpus Number of GPUs available.
 * @return 0 on success, negative value on error.
 *
 */
int ReadCSV(parsec_context_t *apContext, parsec_matrix_block_cyclic_t *apDesc, int aMB, int aNB, int aNodes,
            int aTimeSlot, char *apFilename, int aRank, int aVerbose, int aGpus);

/**
 * @brief Reads a specific time slot from a CSV file.
 * @details This function extracts data from a specified time slot in a CSV file and updates the matrix description.
 * @param[in] apContext Pointer to the parsec context.
 * @param[in, out] apDesc Pointer to the matrix block cyclic descriptor.
 * @param[in] aMB Number of rows in the block.
 * @param[in] aNB Number of columns in the block.
 * @param[in] aNodes Number of nodes.
 * @param[in] aTimeSlot Time slot for the data.
 * @param[in] apFilename Filename of the CSV file.
 * @param[in] aRank Rank of the current process.
 * @param[in] aVerbose Verbosity level for output.
 * @param[in] aGpus Number of GPUs available.
 * @return 0 on success, negative value on error.
 *
 */
int ReadCSVTimeSlot(parsec_context_t *apContext, parsec_matrix_block_cyclic_t *apDesc, int aMB, int aNB, int aNodes,
            int aTimeSlot, char *apFilename, int aRank, int aVerbose, int aGpus);

/**
 * @brief Reads a complex CSV file for a specific time slot.
 * @details This function reads complex data from a CSV file corresponding to a specific time slot.
 * @param[in] apContext Pointer to the parsec context.
 * @param[in, out] apDesc Pointer to the matrix block cyclic descriptor.
 * @param[in] aMB Number of rows in the block.
 * @param[in] aNB Number of columns in the block.
 * @param[in] aNodes Number of nodes.
 * @param[in] aTimeSlot Time slot for the data.
 * @param[in] apFilename Filename of the CSV file.
 * @param[in] aRank Rank of the current process.
 * @param[in] aVerbose Verbosity level for output.
 * @param[in] aGpus Number of GPUs available.
 * @return 0 on success, negative value on error.
 *
 */
int ReadCSVToComplexTimeSlot(parsec_context_t *apContext, parsec_matrix_block_cyclic_t *apDesc, int aMB, int aNB, int aNodes,
            int aTimeSlot, char *apFilename, int aRank, int aVerbose, int aGpus);

/**
 * @brief Reads a complex CSV file.
 * @details This function reads complex data from a CSV file and updates the matrix description.
 * @param[in] apContext Pointer to the parsec context.
 * @param[in, out] apDesc Pointer to the matrix block cyclic descriptor.
 * @param[in] aMB Number of rows in the block.
 * @param[in] aNB Number of columns in the block.
 * @param[in] aNodes Number of nodes.
 * @param[in] aTimeSlot Time slot for the data.
 * @param[in] apFilename Filename of the CSV file.
 * @param[in] aRank Rank of the current process.
 * @param[in] aVerbose Verbosity level for output.
 * @param[in] aGpus Number of GPUs available.
 * @return 0 on success, negative value on error.
 *
 */
int ReadCSVComplex(parsec_context_t *apContext, parsec_matrix_block_cyclic_t *apDesc, int aMB, int aNB, int aNodes,
            int aTimeSlot, char *apFilename, int aRank, int aVerbose, int aGpus);

/**
 * @brief Reads a complex CSV file into a matrix description.
 * @details This function reads complex data from a CSV file and populates the specified matrix description.
 * @param[in] apContext Pointer to the parsec context.
 * @param[in, out] apDesc Pointer to the matrix block cyclic descriptor.
 * @param[in] aMB Number of rows in the block.
 * @param[in] aNB Number of columns in the block.
 * @param[in] aNodes Number of nodes.
 * @param[in] aTimeSlot Time slot for the data.
 * @param[in] apFilename Filename of the CSV file.
 * @param[in] aRank Rank of the current process.
 * @param[in] aVerbose Verbosity level for output.
 * @param[in] aGpus Number of GPUs available.
 * @return 0 on success, negative value on error.
 *
 */
int ReadCSVToComplex(parsec_context_t *apContext, parsec_matrix_block_cyclic_t *apDesc, int aMB, int aNB, int aNodes,
            int aTimeSlot, char *apFilename, int aRank, int aVerbose, int aGpus);

/**
 * @brief Performs forward spherical harmonic transform.
 * @details This function computes the forward spherical harmonic transform using the provided data descriptors.
 * @param[in] apContext Pointer to the parsec context.
 * @param[in] apFDataDesc Pointer to the data descriptor for the forward transform.
 * @param[in] apFLMDesc Pointer to the descriptor for the spherical harmonic coefficients.
 * @param[in] apFLMTDesc Pointer to the descriptor for the transformed spherical harmonic coefficients.
 * @param[in] apET1Desc Pointer to the first auxiliary descriptor.
 * @param[in] apET2Desc Pointer to the second auxiliary descriptor.
 * @param[in] apEPDesc Pointer to the endpoint descriptor.
 * @param[in] apSLMNDesc Pointer to the descriptor for Spherical Harmonic Matrix.
 * @param[in] apIEDesc Pointer to the input descriptor.
 * @param[in] apIODesc Pointer to the output descriptor.
 * @param[in] apPDesc Pointer to the parameter descriptor.
 * @param[in] apDDesc Pointer to the descriptor for intermediate data.
 * @param[in] aFDataM Number of rows in the forward data matrix.
 * @param[in] aEPN Number of endpoints.
 * @param[in] aET1M Number of rows in the first auxiliary matrix.
 * @param[in] aET2M Number of rows in the second auxiliary matrix.
 * @param[in] aPN Number of processes.
 * @param[in] aFlmM Number of rows in the spherical harmonic coefficients matrix.
 * @param[in] aFlmN Number of columns in the spherical harmonic coefficients matrix.
 * @param[in] aLSize Size of the spherical harmonic basis.
 * @return 0 on success, negative value on error.
 *
 */
int ForwardSHT(parsec_context_t *apContext, parsec_tiled_matrix_t *apFDataDesc, parsec_tiled_matrix_t *apFLMDesc,
               parsec_tiled_matrix_t *apFLMTDesc, parsec_tiled_matrix_t *apET1Desc,
               parsec_tiled_matrix_t *apET2Desc, parsec_tiled_matrix_t *apEPDesc,
               parsec_tiled_matrix_t *apSLMNDesc, parsec_tiled_matrix_t *apIEDesc,
               parsec_tiled_matrix_t *apIODesc, parsec_tiled_matrix_t *apPDesc,
               parsec_tiled_matrix_t *apDDesc, int aFDataM, int aEPN, int aET1M,
               int aET2M, int aPN, int aFlmM, int aFlmN, int aLSize);

/**
 * @brief Performs forward spherical harmonic transform with reshaping.
 * @details This function computes the forward spherical harmonic transform and reshapes the resulting data.
 * @param[in] apContext Pointer to the parsec context.
 * @param[in] aRank Rank of the current process.
 * @param[in] aVerbose Verbosity level for output.
 * @param[in] apFDataDesc Pointer to the data descriptor for the forward transform.
 * @param[in] apFLMDesc Pointer to the descriptor for the spherical harmonic coefficients.
 * @param[in] apFLMTDesc Pointer to the descriptor for the transformed spherical harmonic coefficients.
 * @param[in] apET1Desc Pointer to the first auxiliary descriptor.
 * @param[in] apET2Desc Pointer to the second auxiliary descriptor.
 * @param[in] apEPDesc Pointer to the endpoint descriptor.
 * @param[in] apSLMNDesc Pointer to the descriptor for Spherical Harmonic Matrix.
 * @param[in] apIEDesc Pointer to the input descriptor.
 * @param[in] apIODesc Pointer to the output descriptor.
 * @param[in] apPDesc Pointer to the parameter descriptor.
 * @param[in] apDDesc Pointer to the descriptor for intermediate data.
 * @param[in] aFDataM Number of rows in the forward data matrix.
 * @param[in] aEPN Number of endpoints.
 * @param[in] aET1M Number of rows in the first auxiliary matrix.
 * @param[in] aET2M Number of rows in the second auxiliary matrix.
 * @param[in] aPN Number of processes.
 * @param[in] aFlmM Number of rows in the spherical harmonic coefficients matrix.
 * @param[in] aFlmN Number of columns in the spherical harmonic coefficients matrix.
 * @param[in] aLSize Size of the spherical harmonic basis.
 * @return 0 on success, negative value on error.
 *
 */
int ForwardSHTReshape(parsec_context_t *apContext, int aRank, int aVerbose, parsec_tiled_matrix_t *apFDataDesc,
                      parsec_tiled_matrix_t *apFLMDesc, parsec_tiled_matrix_t *apFLMTDesc,
                      parsec_tiled_matrix_t *apET1Desc, parsec_tiled_matrix_t *apET2Desc,
                      parsec_tiled_matrix_t *apEPDesc, parsec_tiled_matrix_t *apSLMNDesc,
                      parsec_tiled_matrix_t *apIEDesc, parsec_tiled_matrix_t *apIODesc,
                      parsec_tiled_matrix_t *apPDesc, parsec_tiled_matrix_t *apDDesc,
                      parsec_tiled_matrix_t *apADesc, int aFDataM, int aEPN, int aET1M,
                      int aET2M, int aPN, int aFlmTNB, int aT, int aLSize, double aNormGlobal,
                      int aNT, int aUpperLower);

/**
 * @brief Computes the norm of a tiled matrix.
 * @details This function calculates the norm of a matrix represented by the given tiled matrix descriptor.
 * The result can be computed for the entire matrix or for specific tiles based on the provided parameters.
 * It can also handle different matrix formats (upper/lower triangular).
 * @param[in] apContext Pointer to the parsec context, which holds information about the execution environment.
 * @param[in] apADesc Pointer to the descriptor of the tiled matrix for which the norm is being computed.
 * @param[in] aNormGlobal The global norm value to be updated. It should be initialized before calling this function.
 * @param[in] aNT The number of tiles in the tiled matrix. This indicates how many submatrices will be considered.
 * @param[in] aUpperLower Specifies whether to calculate the norm for the upper or lower triangular part of the matrix.
 * @param[out] apNormTile Pointer to an array where the computed norms for each tile will be stored.
 * The size of this array should be at least `aNT` to hold norms for all tiles.
 * @return The computed matrix norm. Returns a non-negative value on success, and a negative value if an error occurs.
 *
 */
double GetMatrixNorm(parsec_context_t *apContext, parsec_tiled_matrix_t *apADesc, double aNormGlobal, int aNT,
                     int aUpperLower, double *apNormTile);

/**
 * @brief Performs inverse spherical harmonic transform.
 * @details This function computes the inverse spherical harmonic transform using the provided data descriptors.
 * @param[in] apContext Pointer to the parsec context.
 * @param[in] apFSpatialDesc Pointer to the descriptor for the spherical harmonic coefficients.
 * @param[in] apFLMTDesc Pointer to the descriptor for the transformed spherical harmonic coefficients.
 * @param[in] apZLMDesc Pointer to the second auxiliary descriptor.
 * @param[in] apSCDesc Pointer to the input descriptor.
 * @param[in] aLSize Size of the spherical harmonic basis.
 * @return 0 on success, negative value on error.
 *
 */
int InverseSHT(parsec_context_t *apContext, parsec_tiled_matrix_t *apFSpatialDesc, parsec_tiled_matrix_t *apFLMDesc,
               parsec_tiled_matrix_t *apZLMDesc, parsec_tiled_matrix_t *apSCDesc, int aLSize);


