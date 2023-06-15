
// Copyright (c) 2017-2023 King Abdullah University of Science and Technology,
// All rights reserved.
// ExaGeoStat is a software package, provided by King Abdullah University of Science and Technology (KAUST).

/**
 * @file HicmaImplementation.cpp
 * @brief Sets up the HiCMA descriptors needed for the tile low rank computations in ExaGeoStat.
 * @version 1.0.0
 * @author Sameh Abdulah
 * @date 2023-03-26
**/

#include <lapacke.h>
#include <linear-algebra-solvers/concrete/tile-low-rank/HicmaImplementation.hpp>

extern "C" {
#include <hicma.h>
#include <control/hicma_context.h>
}

using namespace exageostat::linearAlgebra::tileLowRank;
using namespace exageostat::common;
using namespace exageostat::dataunits;
using namespace exageostat::kernels;

using namespace std;

template<typename T>
void HicmaImplementation<T>::InitiateDescriptors() {

    // Check for Initialise the Hicma context.
    if (!this->apContext) {
        throw std::runtime_error(
                "ExaGeoStat hardware is not initialized, please use 'ExaGeoStat<double/float>::ExaGeoStatInitializeHardware(configurations)'.");
    }
    vector<void *> &pDescriptorC = this->mpConfigurations->GetDescriptorC();
    vector<void *> &pDescriptorZ = this->mpConfigurations->GetDescriptorZ();
    auto pDescriptorZcpy = &this->mpConfigurations->GetDescriptorZcpy();
    auto pDescriptorDeterminant = &this->mpConfigurations->GetDescriptorDeterminant();

    // Create a Hicma sequence
    HICMA_sequence_t *pSequence;
    HICMA_request_t request[2] = {HICMA_SUCCESS, HICMA_SUCCESS};
    HICMA_Sequence_Create(&pSequence);

    int N = this->mpConfigurations->GetProblemSize();
    int lts = this->mpConfigurations->GetLowTileSize();
    int pGrid = this->mpConfigurations->GetPGrid();
    int qGrid = this->mpConfigurations->GetQGrid();
    bool isOOC = this->mpConfigurations->GetIsOOC();
    int maxRank = this->mpConfigurations->GetMaxRank();
    int nZmiss = this->mpConfigurations->GetUnknownObservationsNb();
    T meanSquareError = this->mpConfigurations->GetMeanSquareError();
    int approximationMode = this->mpConfigurations->GetApproximationMode();
    string actualObservationsFilePath = this->mpConfigurations->GetActualObservationsFilePath();
    T determinantValue = this->mpConfigurations->GetDeterminantValue();

    int nZobsValue;
    if (actualObservationsFilePath.empty()) {
        nZobsValue = N - nZmiss;
    } else {
        nZobsValue = N;
    }

    this->mpConfigurations->SetKnownObservationsValues(nZobsValue);
    int nZobs = this->mpConfigurations->GetKnownObservationsValues();

    // For distributed system and should be removed
    T *Zcpy = (T *) malloc(N * sizeof(T));

    pDescriptorC.push_back(nullptr);
    auto **pHicmaDescriptorC = &pDescriptorC[0];

    pDescriptorZ.push_back(nullptr);
    auto **pHicmaDescriptorZ = &pDescriptorZ[0];

    auto **pDescriptorCD = &this->mpConfigurations->GetDescriptorCD()[0];
    auto **pDescriptorC12D = &this->mpConfigurations->GetDescriptorCD()[1];
    auto **pDescriptorC22D = &this->mpConfigurations->GetDescriptorCD()[2];

    auto **pDescriptorCUV = &this->mpConfigurations->GetDescriptorCUV()[0];
    auto **pDescriptorC12UV = &this->mpConfigurations->GetDescriptorCUV()[1];
    auto **pDescriptorC22UV = &this->mpConfigurations->GetDescriptorCUV()[2];

    auto **pDescriptorCrk = &this->mpConfigurations->GetDescriptorCrk()[0];
    auto **pDescriptorC12rk = &this->mpConfigurations->GetDescriptorCrk()[1];
    auto **pDescriptorC22rk = &this->mpConfigurations->GetDescriptorCrk()[2];

    auto **pDescriptorZObservations = &this->mpConfigurations->GetDescriptorZObservations();
    auto **pDescriptorZactual = &this->mpConfigurations->GetDescriptorZActual();
    auto **pDescriptorMSE = &this->mpConfigurations->GetDescriptorMSE();

    int MBC, NBC, MC, NC;
    int MBD, NBD, MD, ND;
    int MBUV, NBUV, MUV, NUV;
    int MBrk, NBrk, Mrk, Nrk;

    FloatPoint floatPoint;
    if (sizeof(T) == SIZE_OF_FLOAT) {
        floatPoint = EXAGEOSTAT_REAL_FLOAT;
    } else {
        floatPoint = EXAGEOSTAT_REAL_DOUBLE;
    }

    //CDense Descriptor
    if (approximationMode == 1) {
        MBC = lts;
        NBC = lts;
        MC = N;
        NC = N;
    } else {
        MBC = 1;
        NBC = 1;
        MC = lts;
        NC = lts;
    }

    ExageostatAllocateMatrixTile(pHicmaDescriptorC, isOOC, nullptr, (HICMA_enum) floatPoint, MBC, NBC,
                                 MBC * NBC, MC, NC, 0, 0, MC, NC, pGrid, qGrid);

    //CAD Descriptor
    MBD = lts;
    NBD = lts;
    MD = N;
    ND = MBD;
    ExageostatAllocateMatrixTile(pDescriptorCD, isOOC, nullptr, (HICMA_enum) floatPoint, MBD, NBD, MBD * NBD,
                                 MD, ND, 0, 0, MD, ND, pGrid, qGrid);

    //CUV Descriptor
    MBUV = lts;
    NBUV = 2 * maxRank;
    int N_over_lts_times_lts = N / lts * lts;
    if (N_over_lts_times_lts < N) {
        MUV = N_over_lts_times_lts + lts;
    } else if (N_over_lts_times_lts == N) {
        MUV = N_over_lts_times_lts;
    } else {
        throw range_error("Invalid value. This case should not happen, Please make sure of N and lts values.");
    }

    T expr = (T) MUV / (T) lts;
    NUV = 2 * expr * maxRank;
    ExageostatAllocateMatrixTile(pDescriptorCUV, isOOC, nullptr, (HICMA_enum) floatPoint, MBUV, NBUV,
                                 MBUV * NBUV, MUV, NUV, 0, 0, MUV, NUV, pGrid, qGrid);

    //Crk Descriptor
    MBrk = 1;
    NBrk = 1;
    auto **desc_cuv = (HICMA_desc_t **) pDescriptorCUV;
    Mrk = (*desc_cuv)->mt;
    Nrk = (*desc_cuv)->mt;
    ExageostatAllocateMatrixTile(pDescriptorCrk, isOOC, nullptr, (HICMA_enum) floatPoint, MBrk, NBrk,
                                 MBrk * NBrk, Mrk, Nrk, 0, 0, Mrk, Nrk, pGrid, qGrid);

    HICMA_Sequence_Create(&pSequence);
    ExageostatAllocateMatrixTile(pHicmaDescriptorZ, isOOC, nullptr, (HICMA_enum) floatPoint, lts, lts,
                                 lts * lts, N, 1, 0, 0, N, 1, pGrid, qGrid);

    ExageostatAllocateMatrixTile(pDescriptorZcpy, isOOC, Zcpy, (HICMA_enum) floatPoint, lts, lts, lts * lts,
                                 N, 1, 0, 0, N, 1, pGrid, qGrid);
    ExageostatAllocateMatrixTile(pDescriptorDeterminant, isOOC, &determinantValue, (HICMA_enum) floatPoint,
                                 lts, lts, lts * lts, 1, 1, 0, 0, 1, 1, pGrid, qGrid);

    if (nZmiss != 0) {
        if (actualObservationsFilePath.empty()) {
            ExageostatAllocateMatrixTile(pDescriptorZObservations, isOOC, &Zcpy[nZmiss],
                                         (HICMA_enum) floatPoint, lts, lts, lts * lts, nZobs, 1, 0, 0, nZobs,
                                         1, pGrid, qGrid);
            ExageostatAllocateMatrixTile(pDescriptorZactual, isOOC, Zcpy, (HICMA_enum) floatPoint,
                                         lts, lts, lts * lts, nZmiss, 1, 0, 0, nZmiss, 1, pGrid, qGrid);
        } else {

            ExageostatAllocateMatrixTile(pDescriptorZObservations, isOOC, nullptr, (HICMA_enum) floatPoint,
                                         lts,
                                         lts, lts * lts, nZmiss, 1, 0, 0, nZmiss, 1, pGrid, qGrid);
            ExageostatAllocateMatrixTile(pDescriptorZactual, isOOC, nullptr, (HICMA_enum) floatPoint, lts,
                                         lts, lts * lts, nZmiss, 1, 0, 0, nZmiss, 1, pGrid, qGrid);
        }
        //C12AD Descriptor
        MBD = lts;
        NBD = lts;
        MD = nZmiss;
        ND = MBD;
        ExageostatAllocateMatrixTile(pDescriptorC12D, isOOC, nullptr, (HICMA_enum) floatPoint, MBD, NBD,
                                     MBD * NBD, MD, ND, 0, 0, MD, ND, pGrid, qGrid);

        //C12UV Descriptor
        MBUV = lts;
        NBUV = 2 * maxRank;
        ExageostatAllocateMatrixTile(pDescriptorC12UV, isOOC, nullptr, (HICMA_enum) floatPoint, MBUV, NBUV,
                                     MBUV * NBUV, MBUV, NBUV, 0, 0, MBUV, NBUV, pGrid, qGrid);

        //C12Ark Descriptor
        MBrk = 1;
        NBrk = 1;
        auto **desc_c12uv = (HICMA_desc_t **) pDescriptorC12UV;
        Mrk = (*desc_c12uv)->mt;
        Nrk = (*desc_c12uv)->mt;
        ExageostatAllocateMatrixTile(pDescriptorC12rk, isOOC, nullptr, (HICMA_enum) floatPoint, MBrk, NBrk,
                                     MBrk * NBrk, Mrk, Nrk, 0, 0, Mrk, Nrk, pGrid, qGrid);

        //C22D Descriptor
        MBD = lts;
        NBD = lts;
        MD = nZobs;
        ND = MBD;
        ExageostatAllocateMatrixTile(pDescriptorC22D, isOOC, nullptr, (HICMA_enum) floatPoint, MBD, NBD,
                                     MBD * NBD, MD, ND, 0, 0, MD, ND, pGrid, qGrid);

        //C22UV Descriptor
        MBUV = lts;
        NBUV = 2 * maxRank;
        ExageostatAllocateMatrixTile(pDescriptorC22UV, isOOC, nullptr, (HICMA_enum) floatPoint, MBUV, NBUV,
                                     MBUV * NBUV, MBUV, NBUV, 0, 0, MBUV, NBUV, pGrid, qGrid);

        //C22Ark Descriptor
        MBrk = 1;
        NBrk = 1;
        auto **desc_c22uv = (HICMA_desc_t **) pDescriptorC22UV;
        Mrk = (*desc_c22uv)->mt;
        Nrk = (*desc_c22uv)->mt;
        ExageostatAllocateMatrixTile(pDescriptorC22rk, isOOC, nullptr, (HICMA_enum) floatPoint, MBrk, NBrk,
                                     MBrk * NBrk, Mrk, Nrk, 0, 0, Mrk, Nrk, pGrid, qGrid);

        //Other descriptors
        ExageostatAllocateMatrixTile(pDescriptorMSE, isOOC, &meanSquareError, (HICMA_enum) floatPoint, lts,
                                     lts, lts * lts, 1, 1, 0, 0, 1, 1, pGrid, qGrid);
    }

    this->mpConfigurations->SetSequence(pSequence);
    this->mpConfigurations->SetRequest(request);

    //stop gsl error handler
    gsl_set_error_handler_off();
    free(Zcpy);
}

template<typename T>
void HicmaImplementation<T>::ExaGeoStatInitContext(const int &apCoresNumber, const int &apGPUs) {

    if (!this->apContext) {
        HICMA_user_tag_size(31, 26);
        HICMA_Init(apCoresNumber, apGPUs);
        this->apContext = hicma_context_self();
    }
}

template<typename T>
void HicmaImplementation<T>::ExaGeoStatFinalizeContext() {

    if (!this->apContext) {
        cout
                << "No initialised context of HiCMA, Please use 'ExaGeoStat<double/or/float>::ExaGeoStatInitializeHardware(configurations);'"
                << endl;
    } else {
        HICMA_Finalize();
        this->apContext = nullptr;
    }
}

template<typename T>
void
HicmaImplementation<T>::CovarianceMatrixCodelet(void *descA, int &uplo, dataunits::Locations *apLocation1,
                                                dataunits::Locations *apLocation2,
                                                dataunits::Locations *apLocation3, double *apLocalTheta,
                                                int aDistanceMetric,
                                                exageostat::kernels::Kernel *apKernel) {

    // Check for Initialise the Hicma context.
    if (!this->apContext) {
        throw std::runtime_error(
                "ExaGeoStat hardware is not initialized, please use 'ExaGeoStat<double/float>::ExaGeoStatInitializeHardware(configurations)'.");
    }

    HICMA_option_t options;
    HICMA_RUNTIME_options_init(&options, (HICMA_context_t * )
    this->apContext,
            (HICMA_sequence_t * )
    this->mpConfigurations->GetSequence(),
            (HICMA_request_t * )
    this->mpConfigurations->GetRequest());

    int tempmm, tempnn;

    auto *HICMA_descA = (HICMA_desc_t *) descA;
    HICMA_desc_t A = *HICMA_descA;
    struct starpu_codelet *cl = &this->cl_dcmg;
    int m, n, m0 = 0, n0 = 0;

    for (n = 0; n < A.nt; n++) {
        tempnn = n == A.nt - 1 ? A.n - n * A.nb : A.nb;
        if (uplo == HicmaUpperLower) {
            m = 0;
        } else {
            m = A.m == A.n ? n : 0;
        }
        for (; m < A.mt; m++) {

            tempmm = m == A.mt - 1 ? A.m - m * A.mb : A.mb;
            m0 = m * A.mb;
            n0 = n * A.nb;

            // Register the data with StarPU
            starpu_insert_task(starpu_mpi_codelet(cl),
                               STARPU_VALUE, &tempmm, sizeof(int),
                               STARPU_VALUE, &tempnn, sizeof(int),
                               STARPU_VALUE, &m0, sizeof(int),
                               STARPU_VALUE, &n0, sizeof(int),
                               STARPU_W, (starpu_data_handle_t) HICMA_RUNTIME_data_getaddr(HICMA_descA, m, n),
                               STARPU_VALUE, &apLocation1, sizeof(dataunits::Locations *),
                               STARPU_VALUE, &apLocation2, sizeof(dataunits::Locations *),
                               STARPU_VALUE, &apLocation3, sizeof(dataunits::Locations *),
                               STARPU_VALUE, &apLocalTheta, sizeof(double *),
                               STARPU_VALUE, &aDistanceMetric, sizeof(int),
                               STARPU_VALUE, &apKernel, sizeof(exageostat::kernels::Kernel *),
                               0);

            auto handle = (starpu_data_handle_t) HICMA_RUNTIME_data_getaddr(HICMA_descA, m, n);
            this->apMatrix = (double *) starpu_variable_get_local_ptr(handle);
        }
    }
    HICMA_RUNTIME_options_ws_free(&options);
    HICMA_RUNTIME_options_finalize(&options, (HICMA_context_t * )
    this->apContext);

    HICMA_Sequence_Wait((HICMA_sequence_t * )
    this->mpConfigurations->GetSequence());

}

template<typename T>
void HicmaImplementation<T>::GenerateObservationsVector(void *descA, Locations *apLocation1,
                                                        Locations *apLocation2, Locations *apLocation3,
                                                        vector<double> aLocalTheta, int aDistanceMetric,
                                                        Kernel *apKernel) {

    // Check for Initialise the Hicma context.
    if (!this->apContext) {
        throw std::runtime_error(
                "ExaGeoStat hardware is not initialized, please use 'ExaGeoStat<double/float>::ExaGeoStatInitializeHardware(configurations)'.");
    }
    int N = this->mpConfigurations->GetProblemSize();

    int seed = this->mpConfigurations->GetSeed();
    int iseed[4] = {seed, seed, seed, 1};

    //nomral random generation of e -- ei~N(0, 1) to generate Z
    auto *Nrand = (double *) malloc(N * sizeof(double));
    LAPACKE_dlarnv(3, iseed, N, Nrand);

    //Generate the co-variance matrix C
    auto *theta = (double *) malloc(aLocalTheta.size() * sizeof(double));
    for (int i = 0; i < aLocalTheta.size(); i++) {
        theta[i] = aLocalTheta[i];
    }
    VERBOSE("Initializing Covariance Matrix (Synthetic Dataset Generation Phase).....")
    int upper_lower = EXAGEOSTAT_LOWER;
    this->CovarianceMatrixCodelet(descA, upper_lower, apLocation1, apLocation2, apLocation3, theta,
                                  aDistanceMetric, apKernel);
    free(theta);
    VERBOSE("Done.")

    //Copy Nrand to Z
    VERBOSE("Generate Normal Random Distribution Vector Z (Synthetic Dataset Generation Phase) .....")
    auto **HICMA_descriptorZ = (HICMA_desc_t * *) & this->mpConfigurations->GetDescriptorZ()[0];
    CopyDescriptorZ(*HICMA_descriptorZ, Nrand);
    VERBOSE("Done.")
    free(Nrand);
    //// RESET OF THE IMPLEMENTATION WILL BE ADDED AFTER FINALIZING ALL MODULES WITH EXACT.
}

template<typename T>
void HicmaImplementation<T>::DestoryDescriptors() {

    vector<void *> &pDescriptorC = this->mpConfigurations->GetDescriptorC();
    vector<void *> &pDescriptorZ = this->mpConfigurations->GetDescriptorZ();
    auto pHicmaDescriptorZcpy = (HICMA_desc_t * *) & this->mpConfigurations->GetDescriptorZcpy();
    auto pHicmaDescriptorDeterminant = (HICMA_desc_t * *) & this->mpConfigurations->GetDescriptorDeterminant();

    auto **pDescriptorCD = (HICMA_desc_t * *) & this->mpConfigurations->GetDescriptorCD()[0];
    auto **pDescriptorCUV = (HICMA_desc_t * *) & this->mpConfigurations->GetDescriptorCUV()[0];
    auto **pDescriptorCrk = (HICMA_desc_t * *) & this->mpConfigurations->GetDescriptorCrk()[0];
    auto **pDescriptorZObservations = (HICMA_desc_t * *) & this->mpConfigurations->GetDescriptorZObservations();
    auto **pDescriptorZactual = (HICMA_desc_t * *) & this->mpConfigurations->GetDescriptorZActual();
    auto **pDescriptorMSE = (HICMA_desc_t * *) & this->mpConfigurations->GetDescriptorMSE();


    if (!pDescriptorC.empty() && pDescriptorC[0]) {
        HICMA_Desc_Destroy((HICMA_desc_t * *) & pDescriptorC[0]);
    }
    if (!pDescriptorZ.empty() && pDescriptorZ[0]) {
        HICMA_Desc_Destroy((HICMA_desc_t * *) & pDescriptorZ[0]);
    }

    if (*pHicmaDescriptorZcpy) {
        HICMA_Desc_Destroy(pHicmaDescriptorZcpy);
    }
    if (*pHicmaDescriptorDeterminant) {
        HICMA_Desc_Destroy(pHicmaDescriptorDeterminant);
    }

    if (*pDescriptorCD) {
        HICMA_Desc_Destroy(pDescriptorCD);
    }
    if (*pDescriptorCUV) {
        HICMA_Desc_Destroy(pDescriptorCUV);
    }
    if (*pDescriptorCrk) {
        HICMA_Desc_Destroy(pDescriptorCrk);
    }
    if (*pDescriptorZObservations) {
        HICMA_Desc_Destroy(pDescriptorZObservations);
    }
    if (*pDescriptorZactual) {
        HICMA_Desc_Destroy(pDescriptorZactual);
    }
    if (*pDescriptorMSE) {
        HICMA_Desc_Destroy(pDescriptorMSE);
    }

    if ((HICMA_sequence_t * ) this->mpConfigurations->GetSequence()) {
        HICMA_Sequence_Destroy((HICMA_sequence_t * )
        this->mpConfigurations->GetSequence());
    }
}

namespace exageostat::linearAlgebra::tileLowRank {
    template<typename T> void *HicmaImplementation<T>::apContext = nullptr;
}

template<typename T>
void
HicmaImplementation<T>::CopyDescriptorZ(void *apDescA, double *apDoubleVector) {

}

template<typename T>
void HicmaImplementation<T>::ExageostatAllocateMatrixTile(void **apDescriptor, bool aIsOOC, T *apMemSpace, int aType2,
                                                          int aMB,
                                                          int aNB, int aMBxNB, int aLda, int aN, int aSMB, int aSNB,
                                                          int aM, int aN2, int aP, int aQ) {

    if (aIsOOC && apMemSpace == nullptr && aMB != 1 && aNB != 1) {
        HICMA_Desc_Create_OOC((HICMA_desc_t **) apDescriptor, (HICMA_enum) aType2, aMB, aNB, aMBxNB, aLda, aN, aSMB,
                              aSNB, aM, aN2, aP, aQ);
    } else {
        HICMA_Desc_Create((HICMA_desc_t **) apDescriptor, apMemSpace, (HICMA_enum) aType2, aMB, aNB, aMBxNB, aLda, aN,
                          aSMB, aSNB, aM, aN2, aP, aQ);
    }
}