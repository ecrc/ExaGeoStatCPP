
// Copyright (c) 2017-2024 King Abdullah University of Science and Technology,
// All rights reserved.
// ExaGeoStat is a software package, provided by King Abdullah University of Science and Technology (KAUST).

/**
 * @file JDFHelperFunctions.h
 * @brief A header file for declarations of JDF helper functions.
 * @version 2.0.0
 * @author Mahmoud ElKarargy
 * @date 2024-10-20
**/

#include <complex.h>
#include <runtime/parsec/ParsecHeader.hpp>


int CalculateSingleIndex(int aN, int aM);

double SumDoubleData(double *apData, int aColumn, int aRow);

complex double SumComplexData(complex double *apData, int aColumn, int aRow);

void ForwardSHTHelper(double *apFlm, complex double *apF_data, int aFDataM, int aFDataN,
                      complex double *apEt1, int aEt1M, complex double *apEt2, int aEt2M,
                      complex double *apEp, int aEpM, int aEpN, complex double *apSlmn,
                      int aSlmnN, int aSlmnM, complex double *apIe, int aIeM, int aIeN,
                      complex double *apIo,  int aIoM, int aIoN, complex double *apP,
                      int aPM, int aPN, complex double *apD, complex double *apGmtheta_r,
                      complex double *apFmnm, complex double *apTmp1,
                      complex double *apTmp2, int aL);

void InverseSHTHelper(double *apFlm, double *apF_spatial, double *apZlm, double *apSC,
                      double *apSmt, int aL);

void WriteMatrixToBinaryFile(const char *apFilename, size_t aM, size_t aN, double *apDataMatrix);

void WriteData(parsec_context_t *apContext, const char *apFilename,
               parsec_tiled_matrix_t *apParsecMatrix,
               hicma_parsec_data_t *apData, int symm,
               int aProcessRank, int aNodes);

void ReadData(parsec_context_t *apContext, const char *aFileName,
              parsec_tiled_matrix_t *apParsecMatrix,
              hicma_parsec_data_t *apData, int aP,
              int aQ, int aSymm, int aRank, int aNodes);

size_t PARSECGetAddressCM(parsec_tiled_matrix_t *apDescA, int aM, int aN, int aP, int aQ);

static void CreateMPIDatatype2D(int aProcesses, int aMPIRank, int aM, int aN,
                                int aDistM, int aDistN, int aMb, int aNb,
                                int aPRow, int aPCol, MPI_Datatype aOldType,
                                MPI_Datatype *apNewType);

void MPIWriteDFile(const char *apFilename, int aMPIRank, double *apDataMatrix,
                   MPI_Count aSize, MPI_Datatype aOldType, MPI_Datatype aNewType);

void MPIReadDFile(const char *apFilename, int aMPIRank, double *apDataMatrix,
                  MPI_Count aSize, MPI_Datatype aOldType, MPI_Datatype aNewType);
