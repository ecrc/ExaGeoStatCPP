
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

using namespace std;

using namespace exageostat::linearAlgebra::diagonalSuperTile;
using namespace exageostat::common;
using namespace exageostat::dataunits;
using namespace exageostat::helpers;
using namespace exageostat::configurations;
using namespace exageostat::hardware;

template<typename T>
void
ChameleonImplementationDST<T>::InitiateDescriptors(Configurations &aConfigurations, DescriptorData<T> &aDescriptorData,
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
    ExaGeoStatCreateSequence(&pSequence);

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
void ChameleonImplementationDST<T>::InitiatePredictionDescriptors(Configurations &aConfigurations,
                                                                  ExaGeoStatData<T> &aData) {
    throw std::runtime_error("unimplemented for now");
}

template<typename T>
void
ChameleonImplementationDST<T>::InitiateMloeMmomDescriptors(Configurations &aConfigurations, ExaGeoStatData<T> &aData) {
    throw std::runtime_error("unimplemented for now");

}

template<typename T>
void ChameleonImplementationDST<T>::CovarianceMatrixCodelet(DescriptorData<T> &aDescriptorData, void *apDescriptor,
                                                            const int &aTriangularPart, Locations<T> *apLocation1,
                                                            Locations<T> *apLocation2, Locations<T> *apLocation3,
                                                            T *aLocalTheta, const int &aDistanceMetric,
                                                            const std::string &aKernelName) {

    // Check for initialize the Chameleon context.
    if (!this->mpContext) {
        throw std::runtime_error(
                "ExaGeoStat hardware is not initialized, please use 'ExaGeoStat<double/float>::ExaGeoStatInitializeHardware(configurations)'.");
    }

    RUNTIME_option_t options;
    ExaGeoStatOptionsInit(&options, this->mpContext, aDescriptorData.GetSequence(), aDescriptorData.GetRequest());

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
    ExaGeoStatOptionsFree(&options);
    ExaGeoStatOptionsFinalize(&options, this->mpContext);
    ExaGeoStatSequenceWait((RUNTIME_sequence_t *) aDescriptorData.GetSequence());
    delete pKernel;
}

template<typename T>
void ChameleonImplementationDST<T>::ExaGeoStatGaussianToNonTileAsync(
        DescriptorData<T> &aDescriptorData, void *apDesc, T *apTheta) {
    throw std::runtime_error("unimplemented for now");
}


template<typename T>
void ChameleonImplementationDST<T>::GenerateObservationsVector(Configurations &aConfigurations,
                                                               DescriptorData<T> &aDescriptorData,
                                                               const BaseDescriptor &aDescriptor,
                                                               Locations<T> *apLocation1, Locations<T> *apLocation2,
                                                               Locations<T> *apLocation3, const int &aDistanceMetric) {

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
    this->CovarianceMatrixCodelet(aDescriptorData, p_descriptor, lower_upper, apLocation1, apLocation2, apLocation3,
                                  theta, aDistanceMetric, aConfigurations.GetKernelName());
    delete[] theta;
    VERBOSE("Done.")

    //Copy n_rand to Z
    VERBOSE("Generate Normal Random Distribution Vector Z (Synthetic Dataset Generation Phase) .....")
    auto *CHAM_descriptorZ = aDescriptorData.GetDescriptor(CHAMELEON_DESCRIPTOR, DESCRIPTOR_Z).chameleon_desc;
    CopyDescriptorZ(aDescriptorData, CHAM_descriptorZ, n_rand);
    VERBOSE("Done.")

    //Cholesky factorization for the Co-variance matrix C
    VERBOSE("Cholesky factorization of Sigma (Synthetic Dataset Generation Phase) .....")
    int potential_failure = ExaGeoStatPotrfTile(EXAGEOSTAT_LOWER, p_descriptor);
    FAILURE_LOGGER(potential_failure, "Factorization cannot be performed..\nThe matrix is not positive definite")
    VERBOSE("Done.")

    //Triangular matrix-matrix multiplication
    VERBOSE("Triangular matrix-matrix multiplication Z=L.e (Synthetic Dataset Generation Phase) .....")
    ExaGeoStatTrmmTile(EXAGEOSTAT_LEFT, EXAGEOSTAT_LOWER, EXAGEOSTAT_NO_TRANS, EXAGEOSTAT_NON_UNIT, 1, p_descriptor, CHAM_descriptorZ);

    VERBOSE("Done.")

    const int P = aConfigurations.GetP();
    if (aConfigurations.GetLogger()) {
        T *pMatrix;
        VERBOSE("Writing generated data to the disk (Synthetic Dataset Generation Phase) .....")
#ifdef CHAMELEON_USE_MPI
        pMatrix = new T[n];
        ExaGeoStatDesc2Lap( *CHAM_descriptorZ, pMatrix, N, EXAGEOSTAT_UPPER_LOWER);
        if ( CHAMELEON_My_Mpi_Rank() == 0 ){
            DiskWriter<T>::WriteVectorsToDisk(pMatrix, &N, &P, configurations->GetLoggerPath(), apLocation1);
        }
        delete[] pMatrix;
#else
        pMatrix = (T *) CHAM_descriptorZ->mat;
        string path = aConfigurations.GetLoggerPath();
        DiskWriter<T>::WriteVectorsToDisk(*pMatrix, n, P, path, *apLocation1);
#endif
        VERBOSE("Done.")
    }

    ExaGeoStatLaSetTile(EXAGEOSTAT_UPPER_LOWER, 0, 0, p_descriptor);
    delete[] n_rand;
    VERBOSE("Done Z Vector Generation Phase. (Chameleon Synchronous)")
}

template<typename T>
void
ChameleonImplementationDST<T>::CopyDescriptorZ(DescriptorData<T> &aDescriptorData, void *apDescA,
                                               T *apDoubleVector) {

    // Check for initialize the Chameleon context.
    if (!this->mpContext) {
        throw std::runtime_error(
                "ExaGeoStat hardware is not initialized, please use 'ExaGeoStat<double/float>::ExaGeoStatInitializeHardware(configurations)'.");
    }

    RUNTIME_option_t options;
    ExaGeoStatOptionsInit(&options, this->mpContext, aDescriptorData.GetSequence(), aDescriptorData.GetRequest());

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
    ExaGeoStatOptionsFree(&options);

}

template<typename T>
T ChameleonImplementationDST<T>::ExaGeoStatMLETile(const ExaGeoStatHardware &aHardware, ExaGeoStatData<T> &aData,
                                                   Configurations &aConfigurations, const double *theta,
                                                   T *apMeasurementsMatrix) {
    throw std::runtime_error("unimplemented for now");
}

template<typename T>
T *
ChameleonImplementationDST<T>::ExaGeoStatMLEPredictTile(ExaGeoStatData<T> &aData, T *apTheta, const int &aZMissNumber,
                                                        const int &aZObsNumber, T *apZObs, T *apZActual, T *apZMiss,
                                                        const hardware::ExaGeoStatHardware &aHardware,
                                                        Configurations &aConfiguration, Locations<T> &aMissLocations,
                                                        Locations<T> &aObsLocations) {
    throw std::runtime_error("unimplemented for now");
}

template<typename T>
int ChameleonImplementationDST<T>::ExaGeoStatLapackCopyTile(const UpperLower &aUpperLower, void *apA, void *apB) {
    throw std::runtime_error("unimplemented for now");

}

template<typename T>
int
ChameleonImplementationDST<T>::ExaGeoStatLapackToDescriptor(const UpperLower &aUpperLower, void *apAf77,
                                                            const int &aLda, void *apA) {
    throw std::runtime_error("unimplemented for now");
}

template<typename T>
void
ChameleonImplementationDST<T>::ExaGeoStatOptionsInit(void *apOptoins, void *apContext, void *apSequence,
                                                       void *apRequest) {
    RUNTIME_options_init((RUNTIME_option_t *) &apOptoins, (CHAM_context_t *) apContext,
                         (RUNTIME_sequence_t *) apSequence,
                         (RUNTIME_request_t *) apRequest);
}

template<typename T>
void ChameleonImplementationDST<T>::ExaGeoStatOptionsFree(void *apOptions) {
    RUNTIME_options_ws_free((RUNTIME_option_t *) apOptions);
}

template<typename T>
void ChameleonImplementationDST<T>::ExaGeoStatOptionsFinalize(void *apOptions, void *apContext) {
    RUNTIME_options_finalize((RUNTIME_option_t *)apOptions, (CHAM_context_t *) apContext);
}

template<typename T>
int ChameleonImplementationDST<T>::ExaGeoStatSequenceWait(void *apSequence) {
    throw std::runtime_error("unimplemented for now");
}


template<typename T>
int ChameleonImplementationDST<T>::ExaGeoStatCreateSequence(void *apSequence) {
    int status = CHAMELEON_Sequence_Create(((RUNTIME_sequence_t ** )apSequence));
    if(status != CHAMELEON_SUCCESS){
        throw std::runtime_error("CHAMELEON_Sequence_Create Failed!");
    }
    return status;
}

template<typename T>
int ChameleonImplementationDST<T>::ExaGeoStatPotrfTile(const UpperLower &aUpperLower, void *apA) {
    throw std::runtime_error("unimplemented for now");
}

template<typename T>
int
ChameleonImplementationDST<T>::ExaGeoStatTrsmTile(const Side &aSide, const UpperLower &aUpperLower, const Trans &aTrans,
                                                  const Diag &aDiag, const T &aAlpha, void *apA, void *apB) {
    throw std::runtime_error("unimplemented for now");
}

template<typename T>
int
ChameleonImplementationDST<T>::ExaGeoStatGemmTile(const Trans &aTransA, const Trans &aTransB, const T &aAlpha,
                                                  void *apA, void *apB, const T &aBeta, void *apC) {
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
int ChameleonImplementationDST<T>::ExaGeoStatMLEMseTileAsync(void *apDescZPredict, void *apDescZMiss, void *apDescError,
                                                             void *apSequence, void *apRequest) {
    throw std::runtime_error("unimplemented for now");
}

template<typename T>
int
ChameleonImplementationDST<T>::ExaGeoStatPosvTile(const UpperLower &aUpperLower, void *apA, void *apB) {
    throw std::runtime_error("unimplemented for now");
}

template<typename T>
void ChameleonImplementationDST<T>::ExaGeoStatLap2Desc(T *apA, const int &aLDA, void *apDescA,
                                                       const UpperLower &aUpperLower) {
    throw std::runtime_error("unimplemented for now");
}

template<typename T>
void ChameleonImplementationDST<T>::ExaGeoStatDesc2Lap(T *apA, const int &aLDA, void *apDescA,
                                                       const UpperLower &aUpperLower) {
    throw std::runtime_error("unimplemented for now");
}

template<typename T>
int ChameleonImplementationDST<T>::ExaGeoStatLaSetTile(const common::UpperLower &aUpperLower, T alpha, T beta,
                                                         void *apDescriptor) {
    return CHAMELEON_dlaset_Tile((cham_uplo_t) aUpperLower, alpha, beta, (CHAM_desc_t *)apDescriptor);
}

template<typename T>
void
ChameleonImplementationDST<T>::ExaGeoStatGetZObs(exageostat::configurations::Configurations &aConfigurations, T *apZ,
                                                 const int &aSize,
                                                 exageostat::dataunits::DescriptorData<T> &aDescData,
                                                 T *apMeasurementsMatrix) {
    throw std::runtime_error("unimplemented for now");
}

template<typename T>
void
ChameleonImplementationDST<T>::ExaGeoStatMLEMloeMmomTile(Configurations &aConfigurations, ExaGeoStatData<T> &aData,
                                                         const ExaGeoStatHardware &aHardware, T *apTruthTheta,
                                                         T *apEstimatedTheta, Locations<T> &aMissLocations,
                                                         Locations<T> &aObsLocations) {
    throw std::runtime_error("unimplemented for now");
}

template<typename T>
int
ChameleonImplementationDST<T>::ExaGeoStatMLEMloeMmomTileAsync(void *apDescExpr2, void *apDescExpr3, void *apDescExpr4,
                                                              void *apDescMloe, void *apDescMmom, void *apSequence,
                                                              void *apRequest) {
    throw std::runtime_error("unimplemented for now");
}

template<typename T>
int
ChameleonImplementationDST<T>::ExaGeoStatGeaddTile(const Trans &aTrans, const T &aAlpha, void *apDescA, const T &aBeta,
                                                   void *apDescB) {
    throw std::runtime_error("unimplemented for now");
}

template<typename T>
void
ChameleonImplementationDST<T>::ExaGeoStatTrmmTile(const Side &aSide, const UpperLower &aUpperLower, const Trans &aTrans,
                                                  const Diag &aDiag, const T &alpha, void *apDescA, void *apDescB) {
    throw std::runtime_error("unimplemented for now");
}