
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

#include <runtime/parsec/ParsecHeader.hpp>

int ReadCSV(parsec_context_t *apContext, parsec_matrix_block_cyclic_t *apDesc, int aMB, int aNB, int aNodes,
            int aTimeSlot, char *apFilename, int aRank, int aVerbose, int aGpus);

int ReadCSVTimeSlot(parsec_context_t *apContext, parsec_matrix_block_cyclic_t *apDesc, int aMB, int aNB, int aNodes,
            int aTimeSlot, char *apFilename, int aRank, int aVerbose, int aGpus);

int ReadCSVToComplexTimeSlot(parsec_context_t *apContext, parsec_matrix_block_cyclic_t *apDesc, int aMB, int aNB, int aNodes,
            int aTimeSlot, char *apFilename, int aRank, int aVerbose, int aGpus);

int ReadCSVComplex(parsec_context_t *apContext, parsec_matrix_block_cyclic_t *apDesc, int aMB, int aNB, int aNodes,
            int aTimeSlot, char *apFilename, int aRank, int aVerbose, int aGpus);

int ReadCSVToComplex(parsec_context_t *apContext, parsec_matrix_block_cyclic_t *apDesc, int aMB, int aNB, int aNodes,
            int aTimeSlot, char *apFilename, int aRank, int aVerbose, int aGpus);

int ForwardSHT(parsec_context_t *apContext, parsec_tiled_matrix_t *apFDataDesc, parsec_tiled_matrix_t *apFLMDesc,
               parsec_tiled_matrix_t *apFLMTDesc, parsec_tiled_matrix_t *apET1Desc,
               parsec_tiled_matrix_t *apET2Desc, parsec_tiled_matrix_t *apEPDesc,
               parsec_tiled_matrix_t *apSLMNDesc, parsec_tiled_matrix_t *apIEDesc,
               parsec_tiled_matrix_t *apIODesc, parsec_tiled_matrix_t *apPDesc,
               parsec_tiled_matrix_t *apDDesc, int aFDataM, int aEPN, int aET1M,
               int aET2M, int aPN, int aFlmM, int aFlmN, int aLSize);

int ForwardSHTReshape(parsec_context_t *apContext, int aRank, int aVerbose, parsec_tiled_matrix_t *apFDataDesc,
                       parsec_tiled_matrix_t *apFLMDesc,  parsec_tiled_matrix_t *apFLMTDesc,  parsec_tiled_matrix_t *apET1Desc,
                       parsec_tiled_matrix_t *apET2Desc,  parsec_tiled_matrix_t *apEPDesc,  parsec_tiled_matrix_t *apSLMNDesc,
                       parsec_tiled_matrix_t *apIEDesc,  parsec_tiled_matrix_t *apIODesc,  parsec_tiled_matrix_t *apPDesc,
                       parsec_tiled_matrix_t *apDDesc,  parsec_tiled_matrix_t *apADesc, int aFDataM, int aEPN, int aET1M, int aET2M,
                      int aPN, int aFlmTNB, int aT, int aLSize, double aNormGlobal, int aNT, int aUpperLower);

double GetMatrixNorm(parsec_context_t *apContext, parsec_tiled_matrix_t *apADesc, double aNormGlobal, int aNT,
                     int aUpperLower, double *apNormTile);

int InverseSHT(parsec_context_t *apContext, parsec_tiled_matrix_t *apFSpatialDesc, parsec_tiled_matrix_t *apFLMDesc,
               parsec_tiled_matrix_t *apZLMDesc, parsec_tiled_matrix_t *apSCDesc, int aLSize);


