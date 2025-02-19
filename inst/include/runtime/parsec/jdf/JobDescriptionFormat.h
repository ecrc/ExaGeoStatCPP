
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
int ReadCSVTimeSlot(parsec_context_t *apContext, parsec_matrix_block_cyclic_t *apDesc, int aMB, int aNB,
                    int aNodes, int aTimeSlot, char *apFilename, int aRank, int aVerbose, int aGpus);

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
int ReadCSVToComplexTimeSlot(parsec_context_t *apContext, parsec_matrix_block_cyclic_t *apDesc, int aMB,
                             int aNB, int aNodes, int aTimeSlot, char *apFilename, int aRank, int aVerbose, int aGpus);

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
int ReadCSVComplex(parsec_context_t *apContext, parsec_matrix_block_cyclic_t *apDesc, int aMB, int aNB,
                   int aNodes, int aTimeSlot, char *apFilename, int aRank, int aVerbose, int aGpus);

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
int ReadCSVToComplex(parsec_context_t *apContext, parsec_matrix_block_cyclic_t *apDesc, int aMB, int aNB,
                     int aNodes, int aTimeSlot, char *apFilename, int aRank, int aVerbose, int aGpus);

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
 * @brief Computes the element-wise difference between two matrices.
 * @details This function calculates the difference between corresponding elements
 * of two matrices, `apDescA` and `apDescB`, which are described by a block-cyclic distribution.
 * @param[in] apContext Pointer to the PaRSEC context (`parsec_context_t`) in which the computation will be performed.
 * @param[in] apDescA Pointer to the descriptor of the first matrix.
 * @param[in] apDescB Pointer to the descriptor of the second matrix.
 * @return Returns 0 on successful completion or an error otherwise.
 *
 */
int DifferenceDouble(parsec_context_t *apContext, parsec_matrix_block_cyclic_t *apDescA, parsec_matrix_block_cyclic_t *apDescB);

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
 * @param[in] apADesc Pointer to the descriptor for A matrix.
 * @param[in] aFDataM Number of rows in the forward data matrix.
 * @param[in] aEPN Number of endpoints.
 * @param[in] aET1M Number of rows in the first auxiliary matrix.
 * @param[in] aET2M Number of rows in the second auxiliary matrix.
 * @param[in] aPN Number of processes.
 * @param[in] aFlmTNB Block size parameter specific to the reshaping of spherical harmonic coefficients.
 * @param[in] aT Transform-specific parameter, often used to define temporal or frequency characteristics.
 * @param[in] aLSize Size of the spherical harmonic basis, determining the resolution of the transformation.
 * @param[out] apNormGlobal Pointer to an array storing global normalization values applied during the transform.
 * @param[in] aNT Number of transformations applied, defining the iteration count for computation.
 * @param[in] aUpperLower Specifies the range for the spherical harmonic transformation, either upper or lower spherical components.
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
                      int aET2M, int aPN, int aFlmTNB, int aT, int aLSize, double *apNormGlobal,
                      int aNT, int aUpperLower);

/**
 * @brief Computes the mean squared error between data and spatial descriptors.
 * @details This function calculates the mean squared error (MSE) between two matrix descriptors, which represent
 * the observed data and the spatial data for a given size of spherical harmonic basis.
 * @param[in] apContext Pointer to the parsec context.
 * @param[in] apFDataDesc Pointer to the matrix descriptor for the observed data.
 * @param[in] apFSpatialDesc Pointer to the matrix descriptor for the spatial data.
 * @param[in] aLSize Size of the spherical harmonic basis.
 * @return 0 on success, negative value on error.
 *
 */
int MeanSquaredError(parsec_context_t *apContext, parsec_matrix_block_cyclic_t* apFDataDesc,
                     parsec_matrix_block_cyclic_t* apFSpatialDesc, int aLSize);

/**
 * @brief Computes the norm of a tiled matrix.
 * @details This function calculates the norm of a matrix represented by the given tiled matrix descriptor.
 * The result can be computed for the entire matrix or for specific tiles based on the provided parameters.
 * @param[in] apContext Pointer to the parsec context, which holds information about the execution environment.
 * @param[in] apNormGlobal The global norm value to be updated.
 * @param[in] apADesc Pointer to the descriptor of the tiled matrix for which the norm is being computed.
 * @param[in] aNT The number of tiles in the tiled matrix. This indicates how many submatrices will be considered.
 * @param[in] aIsSymmetric flag to specify if the data matrix is symmetric or not.
 * @param[in] aUpperLower Specifies whether to calculate the norm for the upper or lower triangular part of the matrix.
 * @return void
 *
 */
void GetMatrixNorm(parsec_context_t *apContext, double *apNormGlobal, parsec_tiled_matrix_t *apADesc,
                   int aNT, int aUpperLower, int aIsSymmetric);

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

/**
 * @brief Compresses a matrix using a specified compression strategy and parameters.
 * @details This function applies matrix compression based on the parameters provided.
 * @param[in] apContex Pointer to the Parsec context, which manages task execution and dataflow.
 * @param[in] apNormGlobal The global norm value to be updated.
 * @param[in] aUpperLower Integer flag indicating which triangular part of the matrix to compress.
 * @param[in] aBandSizeDense Band size for the dense region of the matrix in the compression process.
 * @param[in] aNT Number of tiles in one dimension of the matrix.
 * @param[in] aMaxRank Maximum allowable rank for the compressed matrix. Controls the degree of compression.
 * @param[in] aN Total number of columns in the matrix, used to define the matrix dimensions.
 * @param[in] aAdaptiveDecision Flag indicating whether to enable adaptive compression.
 * @param[in] aTolerance Compression tolerance value, used to decide the accuracy of the compressed matrix.
 * @param[in] aSendFullTile Flag indicating whether to transmit full tiles as part of the compression,.
 * @param[in] aAutoBand Flag to enable automatic adjustment of the band size during compression.
 * @param[in] aGpus Number of GPUs to utilize in the compression.
 * @param[in] apHicmaData The HiCMA data struct of descriptors
 * @param[in] apParamsKernel The Starsh struct of kernels
 * @return void
 *
 */
void MatrixCompress(parsec_context_t *apContext, double *apNormGlobal, int aUpperLower, int aBandSizeDense, int aNT,
                    int aMaxRank, int aN, int aAdaptiveDecision, int aTolerance, int aSendFullTile, int aAutoBand,
                    int aGpus, hicma_parsec_data_t *apHicmaData, starsh_params_t *apParamsKernel);