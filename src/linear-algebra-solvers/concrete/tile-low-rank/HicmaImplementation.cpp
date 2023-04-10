
// Copyright (c) 2017-2023 King Abdullah University of Science and Technology,
// Copyright (C) 2023 by Brightskies inc,
// All rights reserved.
// ExaGeoStat is a software package, provided by King Abdullah University of Science and Technology (KAUST).

/**
 * @file HicmaImplementation.cpp
 *
 * @version 1.0.0
 * @author Sameh Abdulah
 * @date 2023-03-26
**/

#include <linear-algebra-solvers/concrete/tile-low-rank/HicmaImplementation.hpp>

extern "C" {
#include <hicma.h>
#include <control/hicma_context.h>
}

using namespace exageostat::linearAlgebra::tileLowRank;
using namespace exageostat::common;
using namespace std;

template<typename T>
void HicmaImplementation<T>::InitiateDescriptors() {

    this->ExaGeoStatInitContext(this->mpConfigurations->GetCoresNumber(), this->mpConfigurations->GetGPUsNumber());

    vector<void *> &pDescriptorC = this->mpConfigurations->GetDescriptorC();
    vector<void *> &pDescriptorZ = this->mpConfigurations->GetDescriptorZ();
    auto pDescriptorZcpy = (HICMA_desc_t **) &this->mpConfigurations->GetDescriptorZcpy();
    auto pDescriptorDeterminant = (HICMA_desc_t **) &this->mpConfigurations->GetDescriptorDeterminant();

    HICMA_sequence_t *pSequence;

    int N = this->mpConfigurations->GetProblemSize() * this->mpConfigurations->GetP();
    int lts = this->mpConfigurations->GetLowTileSize();
    int pGrid = this->mpConfigurations->GetPGrid();
    int qGrid = this->mpConfigurations->GetQGrid();
    bool isOOC = this->mpConfigurations->GetIsOOC();
    int maxRank = this->mpConfigurations->GetMaxRank();
    int nZmiss = this->mpConfigurations->GetUnknownObservationsNb();
    double meanSquareError = this->mpConfigurations->GetMeanSquareError();
    int approximationMode = this->mpConfigurations->GetApproximationMode();
    string actualObservationsFilePath = this->mpConfigurations->GetActualObservationsFilePath();
    double determinantValue = this->mpConfigurations->GetDeterminantValue();

    int nZobsValue;
    if (actualObservationsFilePath.empty()){
        nZobsValue = N - nZmiss;
    }
    else{
        nZobsValue = N;
    }

    this->mpConfigurations->SetKnownObservationsValues(nZobsValue);
    int nZobs = this->mpConfigurations->GetKnownObservationsValues();

    // For distributed system and should be removed
    T *Zcpy = (T *) malloc(N * sizeof(T));

    pDescriptorC.push_back(nullptr);
    auto **pHicmaDescriptorC = (HICMA_desc_t **) &pDescriptorC[0];

    pDescriptorZ.push_back(nullptr);
    auto **pHicmaDescriptorZ = (HICMA_desc_t **) &pDescriptorZ[0];

    auto **pDescriptorCD = (HICMA_desc_t **) &this->mpConfigurations->GetDescriptorCD()[0];
    auto **pDescriptorC12D = (HICMA_desc_t **) &this->mpConfigurations->GetDescriptorCD()[1];
    auto **pDescriptorC22D = (HICMA_desc_t **) &this->mpConfigurations->GetDescriptorCD()[2];

    auto **pDescriptorCUV = (HICMA_desc_t **) &this->mpConfigurations->GetDescriptorCUV()[0];
    auto **pDescriptorC12UV = (HICMA_desc_t **) &this->mpConfigurations->GetDescriptorCUV()[1];
    auto **pDescriptorC22UV = (HICMA_desc_t **) &this->mpConfigurations->GetDescriptorCUV()[2];

    auto **pDescriptorCrk = (HICMA_desc_t **) &this->mpConfigurations->GetDescriptorCrk()[0];
    auto **pDescriptorC12rk = (HICMA_desc_t **) &this->mpConfigurations->GetDescriptorCrk()[1];
    auto **pDescriptorC22rk = (HICMA_desc_t **) &this->mpConfigurations->GetDescriptorCrk()[2];

    auto **pDescriptorZObservations = (HICMA_desc_t **) &this->mpConfigurations->GetDescriptorZObservations();
    auto **pDescriptorZactual = (HICMA_desc_t **) &this->mpConfigurations->GetDescriptorZActual();
    auto **pDescriptorMSE = (HICMA_desc_t **) &this->mpConfigurations->GetDescriptorMSE();

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

    EXAGEOSTAT_ALLOCATE_APPROX_MATRIX_TILE(pHicmaDescriptorC, isOOC, nullptr, (HICMA_enum) floatPoint, MBC, NBC,
                                           MBC * NBC, MC, NC, 0, 0, MC, NC, pGrid, qGrid);
    //CAD Descriptor
    MBD = lts;
    NBD = lts;
    MD = N;
    ND = MBD;
    EXAGEOSTAT_ALLOCATE_APPROX_MATRIX_TILE(pDescriptorCD, isOOC, nullptr, (HICMA_enum) floatPoint, MBD, NBD, MBD * NBD,
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
    EXAGEOSTAT_ALLOCATE_APPROX_MATRIX_TILE(pDescriptorCUV, isOOC, nullptr, (HICMA_enum) floatPoint, MBUV, NBUV,
                                           MBUV * NBUV, MUV, NUV, 0, 0, MUV, NUV, pGrid, qGrid);

    //Crk Descriptor
    MBrk = 1;
    NBrk = 1;
    Mrk = (*pDescriptorCUV)->mt;
    Nrk = (*pDescriptorCUV)->mt;
    EXAGEOSTAT_ALLOCATE_APPROX_MATRIX_TILE(pDescriptorCrk, isOOC, nullptr, (HICMA_enum) floatPoint, MBrk, NBrk,
                                           MBrk * NBrk, Mrk, Nrk, 0, 0, Mrk, Nrk, pGrid, qGrid);

    HICMA_Sequence_Create(&pSequence);
    EXAGEOSTAT_ALLOCATE_APPROX_MATRIX_TILE(pHicmaDescriptorZ, isOOC, nullptr, (HICMA_enum) floatPoint, lts, lts,
                                           lts * lts, N, 1, 0, 0, N, 1, pGrid, qGrid);

    EXAGEOSTAT_ALLOCATE_APPROX_MATRIX_TILE(pDescriptorZcpy, isOOC, Zcpy, (HICMA_enum) floatPoint, lts, lts, lts * lts,
                                           N, 1, 0, 0, N, 1, pGrid, qGrid);
    EXAGEOSTAT_ALLOCATE_APPROX_MATRIX_TILE(pDescriptorDeterminant, isOOC, &determinantValue, (HICMA_enum) floatPoint,
                                           lts, lts, lts * lts, 1, 1, 0, 0, 1, 1, pGrid, qGrid);

    if (nZmiss != 0) {
        if (actualObservationsFilePath.empty()) {
            EXAGEOSTAT_ALLOCATE_APPROX_MATRIX_TILE(pDescriptorZObservations, isOOC, &pDescriptorZcpy[nZmiss],
                                                   (HICMA_enum) floatPoint, lts, lts, lts * lts, nZobs, 1, 0, 0, nZobs,
                                                   1, pGrid, qGrid);
            EXAGEOSTAT_ALLOCATE_APPROX_MATRIX_TILE(pDescriptorZactual, isOOC, pDescriptorZcpy, (HICMA_enum) floatPoint,
                                                   lts, lts, lts * lts, nZmiss, 1, 0, 0, nZmiss, 1, pGrid, qGrid);
        } else {

            EXAGEOSTAT_ALLOCATE_APPROX_MATRIX_TILE(pDescriptorZObservations, isOOC, nullptr, (HICMA_enum) floatPoint, lts,
                                                   lts, lts * lts, nZmiss, 1, 0, 0, nZmiss, 1, pGrid, qGrid);
            EXAGEOSTAT_ALLOCATE_APPROX_MATRIX_TILE(pDescriptorZactual, isOOC, nullptr, (HICMA_enum) floatPoint, lts,
                                                   lts, lts * lts, nZmiss, 1, 0, 0, nZmiss, 1, pGrid, qGrid);
        }
        //C12AD Descriptor
        MBD = lts;
        NBD = lts;
        MD = nZmiss;
        ND = MBD;
        EXAGEOSTAT_ALLOCATE_APPROX_MATRIX_TILE(pDescriptorC12D, isOOC, nullptr, (HICMA_enum) floatPoint, MBD, NBD,
                                               MBD * NBD, MD, ND, 0, 0, MD, ND, pGrid, qGrid);

        //C12UV Descriptor
        MBUV = lts;
        NBUV = 2 * maxRank;
        MUV = nZmiss;
        NUV = 2 * MUV / lts * maxRank;
        EXAGEOSTAT_ALLOCATE_APPROX_MATRIX_TILE(pDescriptorC12UV, isOOC, nullptr, (HICMA_enum) floatPoint, MBUV, NBUV,
                                               MBUV * NBUV, MBUV, NBUV, 0, 0, MBUV, NBUV, pGrid, qGrid);

        //C12Ark Descriptor
        MBrk = 1;
        NBrk = 1;
        Mrk = (*pDescriptorC12UV)->mt;
        Nrk = (*pDescriptorC12UV)->mt;
        EXAGEOSTAT_ALLOCATE_APPROX_MATRIX_TILE(pDescriptorC12rk, isOOC, nullptr, (HICMA_enum) floatPoint, MBrk, NBrk,
                                               MBrk * NBrk, Mrk, Nrk, 0, 0, Mrk, Nrk, pGrid, qGrid);

        //C11AD Descriptor
        MBD = lts;
        NBD = lts;
        MD = nZobs;
        ND = MBD;
        EXAGEOSTAT_ALLOCATE_APPROX_MATRIX_TILE(pDescriptorC22D, isOOC, nullptr, (HICMA_enum) floatPoint, MBD, NBD,
                                               MBD * NBD, MD, ND, 0, 0, MD, ND, pGrid, qGrid);

        //C12UV Descriptor
        MBUV = lts;
        NBUV = 2 * maxRank;
        MUV = nZobs;
        NUV = 2 * MUV / lts * maxRank;
        EXAGEOSTAT_ALLOCATE_APPROX_MATRIX_TILE(pDescriptorC22UV, isOOC, nullptr, (HICMA_enum) floatPoint, MBUV, NBUV,
                                               MBUV * NBUV, MBUV, NBUV, 0, 0, MBUV, NBUV, pGrid, qGrid);

        //C12Ark Descriptor
        MBrk = 1;
        NBrk = 1;
        Mrk = (*pDescriptorC22UV)->mt;
        Nrk = (*pDescriptorC22UV)->mt;
        EXAGEOSTAT_ALLOCATE_APPROX_MATRIX_TILE(pDescriptorC22rk, isOOC, nullptr, (HICMA_enum) floatPoint, MBrk, NBrk,
                                               MBrk * NBrk, Mrk, Nrk, 0, 0, Mrk, Nrk, pGrid, qGrid);

        //Other descriptors
        EXAGEOSTAT_ALLOCATE_APPROX_MATRIX_TILE(pDescriptorMSE, isOOC, &meanSquareError, (HICMA_enum) floatPoint, lts,
                                               lts, lts * lts, 1, 1, 0, 0, 1, 1, pGrid, qGrid);
    }
    this->ExaGeoStatFinalizeContext();
    //stop gsl error handler
    gsl_set_error_handler_off();
}

template<typename T>
void HicmaImplementation<T>::ExaGeoStatInitContext(const int &apCoresNumber, const int &apGPUs) {

    HICMA_context_t *hicmaContext;
    hicmaContext = hicma_context_self();
    if (hicmaContext != nullptr) {
        printf("Another instance of HiCMA is already running...!");
    } else {
        HICMA_user_tag_size(31, 26);
        HICMA_Init(apCoresNumber, apGPUs);
    }
}

template<typename T>
void HicmaImplementation<T>::ExaGeoStatFinalizeContext() {

    HICMA_context_t *hicmaContext;
    hicmaContext = hicma_context_self();
    if (hicmaContext == nullptr) {
        printf("No active instance of HICMA...please use ExaGeoStatInitContext() function to initiate a new instance!\n");
    } else
        HICMA_Finalize();
}
