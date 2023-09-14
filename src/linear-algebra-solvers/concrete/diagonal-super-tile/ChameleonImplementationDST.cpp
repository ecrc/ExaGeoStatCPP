
// Copyright (c) 2017-2023 King Abdullah University of Science and Technology,
// All rights reserved.
// ExaGeoStat is a software package, provided by King Abdullah University of Science and Technology (KAUST).

/**
 * @file ChameleonImplementationDST.cpp
 * @brief Diagonal Super Tile implementation of linear algebra methods.
 * @version 1.0.0
 * @author Mahmoud ElKarargy
 * @author Sameh Abdulah
 * @date 2023-03-20
**/

#include <lapacke.h>

#include <linear-algebra-solvers/concrete/diagonal-super-tile/ChameleonImplementationDST.hpp>
#include <data-units/ExaGeoStatData.hpp>

using namespace std;

using namespace exageostat::linearAlgebra::diagonalSuperTile;
using namespace exageostat::common;
using namespace exageostat::kernels;
using namespace exageostat::dataunits;
using namespace exageostat::helpers;
using namespace exageostat::configurations;

template<typename T>
void ChameleonImplementationDST<T>::InitiateDescriptors(configurations::Configurations &aConfigurations,
                                                        dataunits::DescriptorData<T> &aDescriptorData,
                                                        T *apMeasurementsMatrix) {

    // Check for initialize the Chameleon context.
    if (!this->mpContext) {
        throw std::runtime_error(
                "ExaGeoStat hardware is not initialized, please use 'ExaGeoStat<double/float>::ExaGeoStatInitializeHardware(configurations)'.");
    }
    // Get the problem size and other configuration parameters
    int n = aConfigurations.GetProblemSize();
    int dts = aConfigurations.GetDenseTileSize();
    int p_grid = aConfigurations.GetPGrid();
    int q_grid = aConfigurations.GetQGrid();
    bool is_OOC = aConfigurations.GetIsOOC();

    // For distributed system and should be removed
    T *z_cpy = new T[n];
    T dot_product_value;
    T det_value;

    // Create a Chameleon sequence
    RUNTIME_sequence_t *pSequence;
    RUNTIME_request_t request[2] = {CHAMELEON_SUCCESS, CHAMELEON_SUCCESS};
    CHAMELEON_Sequence_Create(&pSequence);

    // Set the floating point precision based on the template type
    FloatPoint float_point;
    if (sizeof(T) == SIZE_OF_FLOAT) {
        float_point = EXAGEOSTAT_REAL_FLOAT;
    } else if (sizeof(T) == SIZE_OF_DOUBLE) {
        float_point = EXAGEOSTAT_REAL_DOUBLE;
    } else {
        throw runtime_error("Unsupported for now!");
    }

    aDescriptorData.SetDescriptor(common::CHAMELEON_DESCRIPTOR, DESCRIPTOR_C, is_OOC, nullptr, float_point, dts, dts,
                                  dts * dts, n, n, 0, 0, n, n, p_grid, q_grid);
    aDescriptorData.SetDescriptor(common::CHAMELEON_DESCRIPTOR, DESCRIPTOR_Z, is_OOC, nullptr, float_point, dts, dts,
                                  dts * dts, n, 1, 0, 0, n, 1, p_grid, q_grid);
    aDescriptorData.SetDescriptor(common::CHAMELEON_DESCRIPTOR, DESCRIPTOR_Z_COPY, is_OOC, apMeasurementsMatrix,
                                  float_point, dts, dts, dts * dts, n, 1, 0, 0, n, 1, p_grid, q_grid);
    aDescriptorData.SetDescriptor(common::CHAMELEON_DESCRIPTOR, DESCRIPTOR_DETERMINANT, is_OOC, &det_value, float_point,
                                  dts, dts, dts * dts, 1, 1, 0, 0, 1, 1, p_grid, q_grid);
    aDescriptorData.SetDescriptor(common::CHAMELEON_DESCRIPTOR, DESCRIPTOR_PRODUCT, is_OOC, &dot_product_value,
                                  float_point, dts, dts, dts * dts, 1, 1, 0, 0, 1, 1, p_grid, q_grid);

    if (float_point == EXAGEOSTAT_REAL_DOUBLE) {
        aDescriptorData.SetDescriptor(common::CHAMELEON_DESCRIPTOR, DESCRIPTOR_PRODUCT_1, is_OOC, &dot_product_value,
                                      float_point, dts, dts, dts * dts, 1, 1, 0, 0, 1, 1, p_grid, q_grid);
        aDescriptorData.SetDescriptor(common::CHAMELEON_DESCRIPTOR, DESCRIPTOR_PRODUCT_2, is_OOC, &dot_product_value,
                                      float_point, dts, dts, dts * dts, 1, 1, 0, 0, 1, 1, p_grid, q_grid);
    }
    aDescriptorData.SetSequence(pSequence);
    aDescriptorData.SetRequest(request);

    //stop gsl error handler
    gsl_set_error_handler_off();
    delete[] z_cpy;
    aDescriptorData.SetIsDescriptorInitiated(true);
}

template<typename T>
void ChameleonImplementationDST<T>::InitiatePredictionDescriptors(
        configurations::Configurations &aConfigurations, dataunits::ExaGeoStatData<T> &aData) {
    throw std::runtime_error("unimplemented for now");
}

template<typename T>
void ChameleonImplementationDST<T>::CovarianceMatrixCodelet(dataunits::DescriptorData<T> *apDescriptorData,
                                                            void *apDescriptor, int &aTriangularPart,
                                                            Locations<T> *apLocation1, Locations<T> *apLocation2,
                                                            Locations<T> *apLocation3, T *aLocalTheta,
                                                            int aDistanceMetric, const std::string &aKernelName) {

    // Check for initialize the Chameleon context.
    if (!this->mpContext) {
        throw std::runtime_error(
                "ExaGeoStat hardware is not initialized, please use 'ExaGeoStat<double/float>::ExaGeoStatInitializeHardware(configurations)'.");
    }

    RUNTIME_option_t options;
    RUNTIME_options_init(&options, (CHAM_context_t *)
                                 this->mpContext, (RUNTIME_sequence_t *) apDescriptorData->GetSequence(),
                         (RUNTIME_request_t *) apDescriptorData->GetRequest());

    kernels::Kernel<T> *pKernel = exageostat::plugins::PluginRegistry<kernels::Kernel<T >>::Create(aKernelName);

    int tempmm, tempnn;
    auto *CHAM_apDescriptor = (CHAM_desc_t *) apDescriptor;
    CHAM_desc_t A = *CHAM_apDescriptor;

    struct starpu_codelet *cl = &this->cl_dcmg;
    int m, n, m0 = 0, n0 = 0;

    for (n = 0; n < A.nt; n++) {
        tempnn = n == A.nt - 1 ? A.n - n * A.nb : A.nb;
        if (aTriangularPart == ChamUpperLower) {
            m = 0;
        } else {
            m = A.m == A.n ? n : 0;
        }
        for (; m < A.mt; m++) {

            tempmm = m == A.mt - 1 ? A.m - m * A.mb : A.mb;
            m0 = m * A.mb;
            n0 = n * A.nb;

            // Register the data with StarPU
            starpu_insert_task(cl,
                               STARPU_VALUE, &tempmm, sizeof(int),
                               STARPU_VALUE, &tempnn, sizeof(int),
                               STARPU_VALUE, &m0, sizeof(int),
                               STARPU_VALUE, &n0, sizeof(int),
                               STARPU_W, (starpu_data_handle_t) RUNTIME_data_getaddr(CHAM_apDescriptor, m, n),
                               STARPU_VALUE, &apLocation1, sizeof(dataunits::Locations<T> *),
                               STARPU_VALUE, &apLocation2, sizeof(dataunits::Locations<T> *),
                               STARPU_VALUE, &apLocation3, sizeof(dataunits::Locations<T> *),
                               STARPU_VALUE, &aLocalTheta, sizeof(double *),
                               STARPU_VALUE, &aDistanceMetric, sizeof(int),
                               STARPU_VALUE, &pKernel, sizeof(kernels::Kernel<T> *),
                               0);

        }
    }
    RUNTIME_options_ws_free(&options);
    RUNTIME_options_finalize(&options, (CHAM_context_t *) this->mpContext);
    CHAMELEON_Sequence_Wait((RUNTIME_sequence_t *) apDescriptorData->GetSequence());
    delete pKernel;
}

template<typename T>
void ChameleonImplementationDST<T>::ExaGeoStatGaussianToNonTileAsync(
        dataunits::DescriptorData<T> *apDescriptorData, void *apDesc, T *apTheta) {
    throw std::runtime_error("unimplemented for now");
}


template<typename T>
void ChameleonImplementationDST<T>::GenerateObservationsVector(configurations::Configurations &aConfigurations,
                                                               dataunits::DescriptorData<T> *apDescriptorData,
                                                               BaseDescriptor aDescriptor, Locations<T> *apLocation1,
                                                               Locations<T> *apLocation2, Locations<T> *apLocation3,
                                                               int aDistanceMetric) {

    // Check for initialize the Chameleon context.
    if (!this->mpContext) {
        throw std::runtime_error(
                "ExaGeoStat hardware is not initialized, please use 'ExaGeoStat<double/float>::ExaGeoStatInitializeHardware(configurations)'.");
    }
    const int n = aConfigurations.GetProblemSize();
    int seed = aConfigurations.GetSeed();
    int iseed[4] = {seed, seed, seed, 1};
    auto *p_descriptor = aDescriptor.chameleon_desc;

    //nomral random generation of e -- ei~N(0, 1) to generate Z
    auto *n_rand = new T[n];
    LAPACKE_dlarnv(3, iseed, n, (double *) n_rand);

    //Generate the co-variance matrix C
    auto *theta = new T[aConfigurations.GetInitialTheta().size()];
    for (int i = 0; i < aConfigurations.GetInitialTheta().size(); i++) {
        theta[i] = aConfigurations.GetInitialTheta()[i];
    }

    VERBOSE("Initializing Covariance Matrix (Synthetic Dataset Generation Phase).....")
    int lower_upper = EXAGEOSTAT_LOWER;
    this->CovarianceMatrixCodelet(apDescriptorData, p_descriptor, lower_upper, apLocation1, apLocation2, apLocation3,
                                  theta, aDistanceMetric, aConfigurations.GetKernelName());
    delete[] theta;
    VERBOSE("Done.\n")

    //Copy n_rand to Z
    VERBOSE("Generate Normal Random Distribution Vector Z (Synthetic Dataset Generation Phase) .....")
    auto *CHAM_descriptorZ = apDescriptorData->GetDescriptor(CHAMELEON_DESCRIPTOR, DESCRIPTOR_Z).chameleon_desc;
    CopyDescriptorZ(apDescriptorData, CHAM_descriptorZ, n_rand);
    VERBOSE("Done.\n")

    //Cholesky factorization for the Co-variance matrix C
    VERBOSE("Cholesky factorization of Sigma (Synthetic Dataset Generation Phase) .....")
    int potential_failure = CHAMELEON_dpotrf_Tile(ChamLower, p_descriptor);
    FAILURE_LOGGER(potential_failure, "Factorization cannot be performed..\nThe matrix is not positive definite")
    VERBOSE("Done.\n")

    //Triangular matrix-matrix multiplication
    VERBOSE("Triangular matrix-matrix multiplication Z=L.e (Synthetic Dataset Generation Phase) .....")
    CHAMELEON_dtrmm_Tile(ChamLeft, ChamLower, ChamNoTrans, ChamNonUnit, 1, p_descriptor, CHAM_descriptorZ);

    VERBOSE("Done.\n")

    const int P = aConfigurations.GetP();
    if (aConfigurations.GetLogger()) {
        T *pMatrix;
        VERBOSE("Writing generated data to the disk (Synthetic Dataset Generation Phase) .....")
#ifdef CHAMELEON_USE_MPI
        pMatrix = new T[n];
        CHAMELEON_Tile_to_Lapack( *CHAM_descriptorZ, pMatrix, N);
        if ( CHAMELEON_My_Mpi_Rank() == 0 ){
            DiskWriter<T>::WriteVectorsToDisk(pMatrix, &N, &P, configurations->GetLoggerPath(), apLocation1);
        }
        delete[] pMatrix;
#else
        pMatrix = (T *) CHAM_descriptorZ->mat;
        string path = aConfigurations.GetLoggerPath();
        DiskWriter<T>::WriteVectorsToDisk(*pMatrix, n, P, path, *apLocation1);
#endif
        VERBOSE(" Done.\n")
    }

    CHAMELEON_dlaset_Tile(ChamUpperLower, 0, 0, p_descriptor);
    delete[] n_rand;
    VERBOSE("Done Z Vector Generation Phase. (Chameleon Synchronous)")
}

template<typename T>
void
ChameleonImplementationDST<T>::CopyDescriptorZ(dataunits::DescriptorData<T> *apDescriptorData, void *apDescA,
                                               T *apDoubleVector) {

    // Check for initialize the Chameleon context.
    if (!this->mpContext) {
        throw std::runtime_error(
                "ExaGeoStat hardware is not initialized, please use 'ExaGeoStat<double/float>::ExaGeoStatInitializeHardware(configurations)'.");
    }

    RUNTIME_option_t options;
    RUNTIME_options_init(&options, (CHAM_context_t *)
                                 this->mpContext, (RUNTIME_sequence_t *) apDescriptorData->GetSequence(),
                         (RUNTIME_request_t *) apDescriptorData->GetRequest());

    int m, m0;
    int tempmm;
    auto A = (CHAM_desc_t *) apDescA;
    struct starpu_codelet *cl = &this->cl_dzcpy;

    for (m = 0; m < A->mt; m++) {
        tempmm = m == A->mt - 1 ? A->m - m * A->mb : A->mb;
        m0 = m * A->mb;

        starpu_insert_task(cl,
                           STARPU_VALUE, &tempmm, sizeof(int),
                           STARPU_VALUE, &m0, sizeof(int),
                           STARPU_VALUE, &apDoubleVector, sizeof(double),
                           STARPU_W, RUNTIME_data_getaddr(A, m, 0),
#if defined(CHAMELEON_CODELETS_HAVE_NAME)
                STARPU_NAME, "dzcpy",
#endif
                           0);
    }
    RUNTIME_options_ws_free(&options);

}

template<typename T>
T ChameleonImplementationDST<T>::ExaGeoStatMleTile(const hardware::ExaGeoStatHardware &aHardware,
                                                   ExaGeoStatData<T> &aData, Configurations &aConfigurations,
                                                   const double *theta, T *apMeasurementsMatrix) {
    throw std::runtime_error("unimplemented for now");
}

template<typename T>
T *
ChameleonImplementationDST<T>::ExaGeoStatMLEPredictTILE(exageostat::dataunits::ExaGeoStatData<T> &aData, T *apTheta,
                                                        int aZMissNumber, int aZObsNumber,
                                                        T *apZObs, T *apZActual, T *apZMiss,
                                                        const hardware::ExaGeoStatHardware &aHardware,
                                                        configurations::Configurations &aConfiguration,
                                                        exageostat::dataunits::Locations<T> &aMissLocations,
                                                        exageostat::dataunits::Locations<T> &aObsLocations) {
    throw std::runtime_error("unimplemented for now");
}

template<typename T>
int ChameleonImplementationDST<T>::ExaGeoStatLapackCopyTile(exageostat::common::UpperLower aUpperLower, void *apA,
                                                            void *apB) {
    throw std::runtime_error("unimplemented for now");

}

template<typename T>
int
ChameleonImplementationDST<T>::ExaGeoStatLapackToDescriptor(exageostat::common::UpperLower aUpperLower, void *apAf77,
                                                            int aLda, void *apA) {
    throw std::runtime_error("unimplemented for now");
}

template<typename T>
int ChameleonImplementationDST<T>::ExaGeoStatSequenceWait(void *apSequence) {
    throw std::runtime_error("unimplemented for now");
}

template<typename T>
int ChameleonImplementationDST<T>::ExaGeoStatPotrfTile(exageostat::common::UpperLower aUpperLower, void *apA) {
    throw std::runtime_error("unimplemented for now");
}

template<typename T>
int ChameleonImplementationDST<T>::ExaGeoStatTrsmTile(common::Side aSide, common::UpperLower aUpperLower,
                                                      common::Trans aTrans, common::Diag aDiag, T aAlpha, void *apA,
                                                      void *apB) {
    throw std::runtime_error("unimplemented for now");
}

template<typename T>
int ChameleonImplementationDST<T>::ExaGeoStatGemmTile(common::Trans aTransA, common::Trans aTransB, T aAlpha, void *apA,
                                                      void *apB, T aBeta, void *apC) {
    throw std::runtime_error("unimplemented for now");
}

template<typename T>
int ChameleonImplementationDST<T>::ExaGeoStaStrideVectorTileAsync(void *apDescA, void *apDescB, void *apDescC,
                                                                  void *apSequence, void *apRequest) {
    throw std::runtime_error("unimplemented for now");
}

template<typename T>
int ChameleonImplementationDST<T>::ExaGeoStatMeasureDetTileAsync(void *apDescA, void *apSequence, void *apRequest,
                                                                 void *apDescDet) {
    throw std::runtime_error("unimplemented for now");
}

template<typename T>
int ChameleonImplementationDST<T>::ExaGeoStatMleMseTileAsync(void *apDescZPredict, void *apDescZMiss,
                                                             void *apDescError, void *apSequence,
                                                             void *apRequest) {
    throw std::runtime_error("unimplemented for now");
}

template<typename T>
int
ChameleonImplementationDST<T>::ExaGeoStatPosvTile(common::UpperLower aUpperLower, void *apA, void *apB) {
    throw std::runtime_error("unimplemented for now");
}

template<typename T>
void ChameleonImplementationDST<T>::ExaGeoStatLap2Desc(T *apA, int aLDA, void *apDescA,
                                                       common::UpperLower aUpperLower) {
    throw std::runtime_error("unimplemented for now");
}

template<typename T>
void ChameleonImplementationDST<T>::ExaGeoStatDesc2Lap(T *apA, int aLDA, void *apDescA,
                                                       common::UpperLower aUpperLower) {
    throw std::runtime_error("unimplemented for now");
}

template<typename T>
void ChameleonImplementationDST<T>::GetZObs(T *apZ, int aSize, exageostat::dataunits::DescriptorData<T> &aDescData) {
    throw std::runtime_error("unimplemented for now");
}