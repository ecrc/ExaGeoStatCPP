
// Copyright (c) 2017-2023 King Abdullah University of Science and Technology,
// All rights reserved.
// ExaGeoStat is a software package, provided by King Abdullah University of Science and Technology (KAUST).

/**
 * @file ChameleonAllocateDescriptors.cpp
 * @brief Sets up the Chameleon descriptors needed for the diagonal super tile computations in ExaGeoStat.
 * @version 1.0.0
 * @author Sameh Abdulah
 * @date 2023-03-20
**/

#include <lapacke.h>

extern "C" {
#include <chameleon/struct.h>
#include <chameleon.h>
#include <control/context.h>
}

#include <linear-algebra-solvers/concrete/diagonal-super-tile/ChameleonImplementationDST.hpp>
#include "data-units/ExaGeoStatData.hpp"

using namespace std;

using namespace exageostat::linearAlgebra::diagonalSuperTile;
using namespace exageostat::common;
using namespace exageostat::kernels;
using namespace exageostat::dataunits;
using namespace exageostat::helpers;
using namespace exageostat::configurations;

template<typename T>
void ChameleonImplementationDST<T>::InitiateDescriptors(configurations::Configurations *apConfigurations, dataunits::DescriptorData<T> *apDescriptorData) {

//    
//    // Check for Initialise the Chameleon context.
//    if (!this->apContext) {
//        throw std::runtime_error(
//                "ExaGeoStat hardware is not initialized, please use 'ExaGeoStat<double/float>::ExaGeoStatInitializeHardware(configurations)'.");
//    }
//    vector<void *> &pDescriptorC = configurations->GetDescriptorC();
////    vector<void *> &pDescriptorZ = configurations->GetDescriptorZ();
////    auto pChameleonDescriptorZcpy = &configurations->GetDescriptorZcpy();
////    vector<void *> &pDescriptorProduct = configurations->GetDescriptorProduct();
////    auto pChameleonDescriptorDeterminant = &configurations->GetDescriptorDeterminant();
////
//    auto temp= &pDescriptorC;
//    delete temp;
//    pDescriptorC.push_back(nullptr);

//    auto **pChameleonDescriptorC = &pDescriptorC[0];
//
//    pDescriptorZ.push_back(nullptr);
//    auto **pChameleonDescriptorZ = &pDescriptorZ[0];
//
//    // Get the problem size and other configuration parameters
//    int N = configurations->GetProblemSize() * configurations->GetP();
//    int dts = configurations->GetDenseTileSize();
//    int p_grid = configurations->GetPGrid();
//    int q_grid = configurations->GetQGrid();
//    bool is_OOC = configurations->GetIsOOC();
//    int vector_size = 1;
//
//    // For distributed system and should be removed
//    T *Zcpy = (T *) malloc(N * sizeof(T));
//    T dot_product_value;
//
//    // Create a Chameleon sequence
//    RUNTIME_sequence_t *pSequence;
//    RUNTIME_request_t request[2] = {CHAMELEON_SUCCESS, CHAMELEON_SUCCESS};
//    CHAMELEON_Sequence_Create(&pSequence);
//
//    FloatPoint float_point;
//    if (sizeof(T) == SIZE_OF_FLOAT) {
//        float_point = EXAGEOSTAT_REAL_FLOAT;
//    } else {
//        float_point = EXAGEOSTAT_REAL_DOUBLE;
//        vector_size = 3;
//    }
//
//    //Identifies a set of routines sharing common exception handling.
//    CHAMELEON_Sequence_Create(&pSequence);
//
//    ExaGeoStatAllocateMatrixTile(pChameleonDescriptorC, is_OOC, nullptr, (cham_flttype_t) float_point, dts, dts,
//                                 dts * dts, N, N, 0, 0, N, N, p_grid, q_grid);
//    ExaGeoStatAllocateMatrixTile(pChameleonDescriptorZ, is_OOC, nullptr, (cham_flttype_t) float_point, dts, dts,
//                                 dts * dts, N, 1, 0, 0, N, 1, p_grid, q_grid);
//    ExaGeoStatAllocateMatrixTile(pChameleonDescriptorZcpy, is_OOC, Zcpy, (cham_flttype_t) float_point, dts,
//                                 dts, dts * dts, N, 1, 0, 0, N, 1, p_grid, q_grid);
//
//    for (int idx = 0; idx < vector_size; idx++) {
//        pDescriptorProduct.push_back(nullptr);
//        auto **pChameleonDescriptorProduct = &pDescriptorProduct[idx];
//        ExaGeoStatAllocateMatrixTile(pChameleonDescriptorProduct, is_OOC, &dot_product_value,
//                                     (cham_flttype_t) float_point, dts, dts, dts * dts, 1, 1, 0, 0, 1, 1, p_grid,
//                                     q_grid);
//
//    }
//    ExaGeoStatAllocateMatrixTile(pChameleonDescriptorDeterminant, is_OOC, &dot_product_value,
//                                 (cham_flttype_t) float_point, dts, dts, dts * dts, 1, 1, 0, 0, 1, 1, p_grid,
//                                 q_grid);
//
//    configurations->SetSequence(pSequence);
//    configurations->SetRequest(request);
//    //stop gsl error handler
//    gsl_set_error_handler_off();
//    free(Zcpy);
}

template<typename T>
void ChameleonImplementationDST<T>::ExaGeoStatInitContext(int aCoresNumber, int aGPUsNumbers) {

//    if (!this->apContext) {
//        CHAMELEON_user_tag_size(31, 26);
//        CHAMELEON_Init(Configurations::GetConfigurations()->GetCoresNumber(), Configurations::GetConfigurations()->GetGPUsNumbers())
//        this->apContext = chameleon_context_self();
//    }
}

template<typename T>
void ChameleonImplementationDST<T>::ExaGeoStatFinalizeContext() {

    if (!this->apContext) {
        cout
                << "No initialised context of Chameleon, Please use 'ExaGeoStat<double/or/float>::ExaGeoStatInitializeHardware(configurations);'"
                << endl;
    } else {
        CHAMELEON_Finalize()
        this->apContext = nullptr;
    }
}

template<typename T>
void ChameleonImplementationDST<T>::CovarianceMatrixCodelet(dataunits::DescriptorData<T> *apDescriptorData, void *apDescriptor, int &aTriangularPart,
                                                            dataunits::Locations<T> *apLocation1,
                                                            dataunits::Locations<T> *apLocation2,
                                                            dataunits::Locations<T> *apLocation3,
                                                            T *aLocalTheta, int aDistanceMetric,
                                                            kernels::Kernel<T> *apKernel) {

//    // Check for Initialise the Chameleon context.
//    if (!this->apContext) {
//        throw std::runtime_error(
//                "ExaGeoStat hardware is not initialized, please use 'ExaGeoStat<double/float>::ExaGeoStatInitializeHardware(configurations)'.");
//    }
//
//    RUNTIME_option_t options;
//    RUNTIME_options_init(&options, (CHAM_context_t *) this->apContext,
//                         (RUNTIME_sequence_t *) Configurations::GetConfigurations()->GetSequence(),
//                         (RUNTIME_request_t *) Configurations::GetConfigurations()->GetRequest());
//
//    int tempmm, tempnn;
//    auto *CHAM_apDescriptor = (CHAM_desc_t *) apDescriptor;
//    CHAM_desc_t A = *CHAM_apDescriptor;
//
//    struct starpu_codelet *cl = &this->cl_dcmg;
//    int m, n, m0 = 0, n0 = 0;
//
//    for (n = 0; n < A.nt; n++) {
//        tempnn = n == A.nt - 1 ? A.n - n * A.nb : A.nb;
//        if (aTriangularPart == ChamUpperLower) {
//            m = 0;
//        } else {
//            m = A.m == A.n ? n : 0;
//        }
//        for (; m < A.mt; m++) {
//
//            tempmm = m == A.mt - 1 ? A.m - m * A.mb : A.mb;
//            m0 = m * A.mb;
//            n0 = n * A.nb;
//
//            // Register the data with StarPU
//            starpu_insert_task(cl,
//                               STARPU_VALUE, &tempmm, sizeof(int),
//                               STARPU_VALUE, &tempnn, sizeof(int),
//                               STARPU_VALUE, &m0, sizeof(int),
//                               STARPU_VALUE, &n0, sizeof(int),
//                               STARPU_W, (starpu_data_handle_t) RUNTIME_data_getaddr(CHAM_apDescriptor, m, n),
//                               STARPU_VALUE, &apLocation1, sizeof(dataunits::Locations<T> *),
//                               STARPU_VALUE, &apLocation2, sizeof(dataunits::Locations<T> *),
//                               STARPU_VALUE, &apLocation3, sizeof(dataunits::Locations<T> *),
//                               STARPU_VALUE, &aLocalTheta, sizeof(double *),
//                               STARPU_VALUE, &aDistanceMetric, sizeof(int),
//                               STARPU_VALUE, &apKernel, sizeof(kernels::Kernel<T> *),
//                               0);
//
//            auto handle = (starpu_data_handle_t) RUNTIME_data_getaddr(CHAM_apDescriptor, m, n);
//            this->apMatrix = (double *) starpu_variable_get_local_ptr(handle);
//        }
//    }
//    RUNTIME_options_ws_free(&options);
//    RUNTIME_options_finalize(&options, (CHAM_context_t *) this->apContext);
//
//    CHAMELEON_Sequence_Wait((RUNTIME_sequence_t *) Configurations::GetConfigurations()->GetSequence());

}

template<typename T>
void ChameleonImplementationDST<T>::GenerateObservationsVector(configurations::Configurations *apConfigurations, dataunits::DescriptorData<T> *apDescriptorData, BaseDescriptor aDescriptor, Locations<T> *apLocation1,
                                                               Locations<T> *apLocation2, Locations<T> *apLocation3,
                                                               int aDistanceMetric, Kernel<T> *apKernel, T* Nrand) {
//
//    auto configurations = Configurations::GetConfigurations();
//    // Check for Initialise the Chameleon context.
//    if (!this->apContext) {
//        throw std::runtime_error(
//                "ExaGeoStat hardware is not initialized, please use 'ExaGeoStat<double/float>::ExaGeoStatInitializeHardware(configurations)'.");
//    }
//    int N = configurations->GetProblemSize();
//    int seed = configurations->GetSeed();
//    int iseed[4] = {seed, seed, seed, 1};
//
//    //nomral random generation of e -- ei~N(0, 1) to generate Z
//    auto *Nrand = (double *) malloc(N * sizeof(double));
//    LAPACKE_dlarnv(3, iseed, N, Nrand);
//
//    //Generate the co-variance matrix C
//    auto *theta = (double *) malloc(aLocalTheta.size() * sizeof(double));
//    for (int i = 0; i < aLocalTheta.size(); i++) {
//        theta[i] = aLocalTheta[i];
//    }
//
//    VERBOSE("Initializing Covariance Matrix (Synthetic Dataset Generation Phase).....")
//    int lower_upper = EXAGEOSTAT_LOWER;
//    this->CovarianceMatrixCodelet(apDescriptor, lower_upper, apLocation1, apLocation2, apLocation3, theta,
//                                  aDistanceMetric, apKernel);
//
//    free(theta);
//    VERBOSE("Done.\n")
//
//    //Copy Nrand to Z
//    VERBOSE("Generate Normal Random Distribution Vector Z (Synthetic Dataset Generation Phase) .....")
//    auto **CHAM_descriptorZ = (CHAM_desc_t **) &configurations->GetDescriptorZ()[0];
//    CopyDescriptorZ(*CHAM_descriptorZ, Nrand);
//    VERBOSE("Done.\n")
//
//    //Cholesky factorization for the Co-variance matrix C
//    VERBOSE("Cholesky factorization of Sigma (Synthetic Dataset Generation Phase) .....")
//    int potential_failure = CHAMELEON_dpotrf_Tile(ChamLower, (CHAM_desc_t *) apDescriptor);
//    FAILURE_LOGGER(potential_failure, "Factorization cannot be performed..\nThe matrix is not positive definite")
//    VERBOSE("Done.\n")
//
//    //Triangular matrix-matrix multiplication
//    VERBOSE("Triangular matrix-matrix multiplication Z=L.e (Synthetic Dataset Generation Phase) .....")
//    CHAMELEON_dtrmm_Tile(ChamLeft, ChamLower, ChamNoTrans, ChamNonUnit, 1, (CHAM_desc_t *) apDescriptor,
//                         *CHAM_descriptorZ);
//    VERBOSE("Done.\n")
//
//    const int P = configurations->GetP();
//    if (configurations->GetLogger()) {
//        T *pMatrix;
//        VERBOSE("Writing generated data to the disk (Synthetic Dataset Generation Phase) .....")
//#ifdef CHAMELEON_USE_MPI
//        pMatrix = (T*) malloc(N * sizeof(T));
//        CHAMELEON_Tile_to_Lapack( *CHAM_descriptorZ, pMatrix, N);
//        if ( CHAMELEON_My_Mpi_Rank() == 0 ){
//            DiskWriter<T>::WriteVectorsToDisk(pMatrix, &N, &P, configurations->GetLoggerPath(), apLocation1);
//        }
//        free(pMatrix);
//#else
//        pMatrix = (T *) (*CHAM_descriptorZ)->mat;
//        DiskWriter<T>::WriteVectorsToDisk(pMatrix, &N, &P, configurations->GetLoggerPath(), apLocation1);
//        free(pMatrix);
//#endif
//        VERBOSE(" Done.\n")
//    }
//
//    CHAMELEON_dlaset_Tile(ChamUpperLower, 0, 0, (CHAM_desc_t *) apDescriptor);
//    free(Nrand);
//    VERBOSE("Done Z Vector Generation Phase. (Chameleon Synchronous)")
}

template<typename T>
void
ChameleonImplementationDST<T>::CopyDescriptorZ(dataunits::DescriptorData<T> *apDescriptorData, void *apDescA, T *apDoubleVector) {
//
//    // Check for Initialise the Chameleon context.
//    if (!this->apContext) {
//        throw std::runtime_error(
//                "ExaGeoStat hardware is not initialized, please use 'ExaGeoStat<double/float>::ExaGeoStatInitializeHardware(configurations)'.");
//    }
//
//    RUNTIME_option_t options;
//    RUNTIME_options_init(&options, (CHAM_context_t *) this->apContext,
//                         (RUNTIME_sequence_t *) Configurations::GetConfigurations()->GetSequence(),
//                         (RUNTIME_request_t *) Configurations::GetConfigurations()->GetRequest());
//
//    int m, m0;
//    int tempmm;
//    auto A = (CHAM_desc_t *) apDescA;
//    struct starpu_codelet *cl = &this->cl_dzcpy;
//
//    for (m = 0; m < A->mt; m++) {
//        tempmm = m == A->mt - 1 ? A->m - m * A->mb : A->mb;
//        m0 = m * A->mb;
//
//        starpu_insert_task(cl,
//                           STARPU_VALUE, &tempmm, sizeof(int),
//                           STARPU_VALUE, &m0, sizeof(int),
//                           STARPU_VALUE, &apDoubleVector, sizeof(double),
//                           STARPU_W, RUNTIME_data_getaddr(A, m, 0),
//#if defined(CHAMELEON_CODELETS_HAVE_NAME)
//                STARPU_NAME, "dzcpy",
//#endif
//                           0);
//    }
//    RUNTIME_options_ws_free(&options);

}

template<typename T>
void ChameleonImplementationDST<T>::DestroyDescriptors(DescriptorData<T> *apDescriptorData) {
//
//    auto configurations = Configurations::GetConfigurations();
//    vector<void *> &pDescriptorC = configurations->GetDescriptorC();
//    vector<void *> &pDescriptorZ = configurations->GetDescriptorZ();
//    auto pChameleonDescriptorZcpy = (CHAM_desc_t **) &configurations->GetDescriptorZcpy();
//    vector<void *> &pDescriptorProduct = configurations->GetDescriptorProduct();
//    auto pChameleonDescriptorDeterminant = (CHAM_desc_t **) &configurations->GetDescriptorDeterminant();
//
//    if (!pDescriptorC.empty() && pDescriptorC[0]) {
//        CHAMELEON_Desc_Destroy((CHAM_desc_t **) &pDescriptorC[0]);
//    }
//    if (!pDescriptorZ.empty() && pDescriptorZ[0]) {
//        CHAMELEON_Desc_Destroy((CHAM_desc_t **) &pDescriptorZ[0]);
//    }
//    if (!pDescriptorProduct.empty() && pDescriptorProduct[0]) {
//        CHAMELEON_Desc_Destroy((CHAM_desc_t **) &pDescriptorProduct[0]);
//    }
//    if (*pChameleonDescriptorZcpy) {
//        CHAMELEON_Desc_Destroy(pChameleonDescriptorZcpy);
//    }
//    if (*pChameleonDescriptorDeterminant) {
//        CHAMELEON_Desc_Destroy(pChameleonDescriptorDeterminant);
//    }
//
//    if ((RUNTIME_sequence_t *) configurations->GetSequence()) {
//        CHAMELEON_Sequence_Destroy((RUNTIME_sequence_t *) configurations->GetSequence());
//    }
}

template<typename T>
void
ChameleonImplementationDST<T>::ExaGeoStatAllocateMatrixTile(void **apDescriptor, bool ais_OOC, T *apMemSpace,
                                                            int aType2,
                                                            int aMB,
                                                            int aNB, int aMBxNB, int aLda, int aN, int aSMB, int aSNB,
                                                            int aM, int aN2, int aP, int aQ) {
    if (ais_OOC && apMemSpace == nullptr && aMB != 1 && aNB != 1) {
        CHAMELEON_Desc_Create_OOC((CHAM_desc_t **) apDescriptor, (cham_flttype_t) aType2, aMB, aNB, aMBxNB, aLda, aN,
                                  aSMB, aSNB, aM, aN2, aP, aQ);
    } else {
        CHAMELEON_Desc_Create((CHAM_desc_t **) apDescriptor, apMemSpace, (cham_flttype_t) aType2, aMB, aNB, aMBxNB,
                              aLda, aN, aSMB, aSNB, aM, aN2, aP, aQ);
    }
}

template<typename T>
T ChameleonImplementationDST<T>::ExaGeoStatMleTile(ExaGeoStatData<T> *apData, configurations::Configurations *apConfigurations, const double* theta){
    throw std::domain_error("unimplemented for now");
}

template<typename T>
int ChameleonImplementationDST<T>::ExaGeoStatLapackCopyTile(exageostat::common::UpperLower aUpperLower, void *apA, void *apB){
    throw std::domain_error("unimplemented for now");

}

template<typename T>
int ChameleonImplementationDST<T>::ExaGeoStatLapackToDescriptor(exageostat::common::UpperLower aUpperLower, void *apAf77, int aLda, void * apA){
    throw std::domain_error("unimplemented for now");
}

template<typename T>
int ChameleonImplementationDST<T>::ExaGeoStatSequenceWait(void * apSequence) {
    throw std::domain_error("unimplemented for now");
}

template<typename T>
int ChameleonImplementationDST<T>::ExaGeoStatPotrfTile(exageostat::common::UpperLower aUpperLower, void *apA){
    throw std::domain_error("unimplemented for now");
}

template<typename T>
int ChameleonImplementationDST<T>::ExaGeoStatTrsmTile(common::Side aSide, common::UpperLower aUpperLower, common::Trans aTrans, common::Diag aDiag, T aAlpha, void *apA, void *apB) {
    throw std::domain_error("unimplemented for now");
}

template<typename T>
int ChameleonImplementationDST<T>::ExaGeoStatGemmTile(common::Trans aTransA, common::Trans aTransB, T aAlpha,
                                                       void *apA, void *apB, T aBeta, void *apC) {
    throw std::domain_error("unimplemented for now");
}

template<typename T>
int ChameleonImplementationDST<T>::ExaGeoStaStrideVectorTileAsync(void *apDescA, void *apDescB, void *apDescC,
                                                                 void * apSequence, void *apRequest){
    throw std::domain_error("unimplemented for now");
}

template<typename T>
int ChameleonImplementationDST<T>::ExaGeoStatMeasureDetTileAsync(void *apDescA, void * apSequence, void *apRequest, void *apDescDet){
    throw std::domain_error("unimplemented for now");
}

namespace exageostat::linearAlgebra::diagonalSuperTile {
    template<typename T> void *ChameleonImplementationDST<T>::apContext = nullptr;
}
