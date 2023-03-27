
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
#include <hicma_struct.h>
#include <hicma.h>
#include <vector>
#include <gsl/gsl_errno.h>

using namespace exageostat::linearAlgebra::tileLowRank;
using namespace exageostat::common;
using namespace std;

template<typename T>
void HicmaImplementation<T>::InitiateDescriptors() {

    vector<void *> pDescriptorC =  this->mpConfigurations->GetDescriptorC();
    vector<void *> pDescriptorZ = this->mpConfigurations->GetDescriptorZ();
    auto* pDescriptorZcpy = (HICMA_desc_t*) this->mpConfigurations->GetDescriptorZcpy();
    vector<void *> pDescriptorProduct = this->mpConfigurations->GetDescriptorProduct();
    auto* pDescriptorDeterminant = (HICMA_desc_t*) this->mpConfigurations->GetDescriptorDeterminant();


    HICMA_sequence_t *pSequence;
    HICMA_request_t request[2] = {HICMA_SUCCESS, HICMA_SUCCESS};

    int N = this->mpConfigurations->GetProblemSize() * this->mpConfigurations->GetP();
    int lts = this->mpConfigurations->GetLowTileSize();
    int pGrid = this->mpConfigurations->GetPGrid();
    int qGrid = this->mpConfigurations->GetQGrid();
    bool isOOC = this->mpConfigurations->GetIsOOC();
    int maxRank = this->mpConfigurations->GetMaxRank();
    int nZmiss = this->mpConfigurations->GetUnknownObservationsNb();
    int nZobs = this->mpConfigurations->GetKnownObservationsValues();

    // For distributed system and should be removed
    T *Zcpy = (T *) malloc(N * sizeof(T));


//    CHAM_desc_t *descZ = NULL;
//    CHAM_desc_t *cham_descZcpy = NULL;

//    CHAM_desc_t *CHAMELEON_descZmiss = NULL;
//    CHAM_desc_t *CHAM_descmse = NULL;
//    CHAM_desc_t *CHAMELEON_descZactual = NULL;
//    CHAM_desc_t *CHAMELEON_descZobs = NULL;

    pDescriptorC.push_back(nullptr);
    auto* pHicmaDescriptorC = (HICMA_desc_t*) pDescriptorC[0];

    pDescriptorZ.push_back(nullptr);
    auto* pHicmaDescriptorZ = (HICMA_desc_t*) pDescriptorZ[0];

    auto* pDescriptorCD = (HICMA_desc_t*) this->mpConfigurations->GetDescriptorCD()[0];
    auto* pDescriptorC12D = (HICMA_desc_t*) this->mpConfigurations->GetDescriptorCD()[1];
    auto* pDescriptorC22D = (HICMA_desc_t*) this->mpConfigurations->GetDescriptorCD()[2];

    auto* pDescriptorCUV = (HICMA_desc_t*) this->mpConfigurations->GetDescriptorCUV()[0];
    auto* pDescriptorC12UV = (HICMA_desc_t*) this->mpConfigurations->GetDescriptorCUV()[1];
    auto* pDescriptorC22UV = (HICMA_desc_t*) this->mpConfigurations->GetDescriptorCUV()[2];

    auto* pDescriptorCrk = (HICMA_desc_t*) this->mpConfigurations->GetDescriptorCrk()[0];
    auto* pDescriptorC12rk = (HICMA_desc_t*) this->mpConfigurations->GetDescriptorCrk()[1];
    auto* pDescriptorC22rk = (HICMA_desc_t*) this->mpConfigurations->GetDescriptorCrk()[2];

    auto* pDescriptorZObservations = (HICMA_desc_t*) this->mpConfigurations->GetDescriptorZObservations();

    int MBC, NBC, MC, NC;
    int MBD, NBD, MD, ND;
    int MBUV, NBUV, MUV, NUV;
    int MBrk, NBrk, Mrk, Nrk;

    FloatPoint floatPoint = EXAGEOSTAT_REAL_FLOAT;
    if (sizeof(T) == SIZE_OF_FLOAT) {
        floatPoint = EXAGEOSTAT_REAL_FLOAT;
    } else {
        floatPoint = EXAGEOSTAT_REAL_DOUBLE;
    }

    //CDense Descriptor
    if (data->check == 1) {
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
    EXAGEOSTAT_ALLOCATE_APPROX_MATRIX_TILE(&pHicmaDescriptorC, isOOC, nullptr, (HICMA_enum) floatPoint, MBC, NBC, MBC * NBC, MC, NC, 0, 0, MC, NC, pGrid, qGrid);

    //CAD Descriptor
    MBD = lts;
    NBD = lts;
    MD = N;
    ND = MBD;
    EXAGEOSTAT_ALLOCATE_APPROX_MATRIX_TILE(&pDescriptorCD, isOOC, nullptr, (HICMA_enum) floatPoint, MBD, NBD, MBD * NBD, MD, ND, 0, 0, MD, ND, pGrid, qGrid);

    //CAD Descriptor
    MBUV = lts;
    NBUV = 2 * maxRank;
    MUV = -1;
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
    EXAGEOSTAT_ALLOCATE_APPROX_MATRIX_TILE(&pDescriptorCUV, isOOC, nullptr, (HICMA_enum) floatPoint, MBUV, NBUV, MBUV * NBUV, MUV, NUV, 0, 0, MUV, NUV, pGrid, qGrid);

    //CUV Descriptor
    MBrk = 1;
    NBrk = 1;
    Mrk = pDescriptorCUV->mt;
    Nrk = pDescriptorCUV->mt;
    EXAGEOSTAT_ALLOCATE_APPROX_MATRIX_TILE(&pDescriptorCrk, isOOC, nullptr, (HICMA_enum) floatPoint, MBrk, NBrk, MBrk * NBrk, Mrk, Nrk, 0, 0, Mrk, Nrk, pGrid, qGrid);

    HICMA_Sequence_Create(&pSequence);
    EXAGEOSTAT_ALLOCATE_APPROX_MATRIX_TILE(&pHicmaDescriptorZ, isOOC, nullptr, (HICMA_enum) floatPoint, lts, lts, lts * lts, N, 1, 0, 0, N, 1, pGrid, qGrid);
    //// I replaced this with Hicma desc create instead of Chameleon
//    EXAGEOSTAT_ALLOCATE_DENSE_MATRIX_TILE(&descZ, NULL, ChamRealDouble, lts, lts, lts * lts, N, 1, 0, 0, N, 1, p_grid, q_grid);

    EXAGEOSTAT_ALLOCATE_APPROX_MATRIX_TILE(&pDescriptorZcpy, isOOC, nullptr, (HICMA_enum) floatPoint, lts, lts, lts * lts, N, 1, 0, 0, N, 1, pGrid, qGrid);
    EXAGEOSTAT_ALLOCATE_APPROX_MATRIX_TILE(&pDescriptorDeterminant, isOOC, nullptr, (HICMA_enum) floatPoint, lts, lts, lts * lts, 1, 1, 0, 0, 1, 1, pGrid, qGrid);
    //// I replaced this with Hicma desc create instead of Chameleon
//    EXAGEOSTAT_ALLOCATE_DENSE_MATRIX_TILE(&cham_descZcpy, Zcpy, ChamRealDouble, lts, lts, lts * lts, N, 1, 0, 0, N, 1, p_grid, q_grid);

    if (nZmiss != 0) {
//        if (strcmp(data->actualZFPath, "") == 0) {
//            EXAGEOSTAT_ALLOCATE_DENSE_MATRIX_TILE(&CHAMELEON_descZobs, &Zcpy[nZmiss], ChamRealDouble, lts, lts, lts * lts, nZobs, 1, 0, 0, nZobs, 1, p_grid, q_grid);
//            EXAGEOSTAT_ALLOCATE_DENSE_MATRIX_TILE(&CHAMELEON_descZactual, Zcpy, ChamRealDouble, lts, lts, lts * lts, nZmiss, 1, 0, 0, nZmiss, 1, p_grid, q_grid);
//        } else {
//            CHAMELEON_descZobs = cham_descZcpy;
//            EXAGEOSTAT_ALLOCATE_DENSE_MATRIX_TILE(&CHAMELEON_descZactual, NULL, ChamRealDouble, lts, lts, lts * lts, nZmiss, 1, 0, 0, nZmiss, 1, p_grid, q_grid);
//        }

        //C12AD Descriptor
        MBD = lts;
        NBD = lts;
        MD = nZmiss;
        ND = MBD;
        EXAGEOSTAT_ALLOCATE_APPROX_MATRIX_TILE(&pDescriptorC12D, isOOC, nullptr, (HICMA_enum) floatPoint, MBD, NBD, MBD * NBD, MD, ND, 0, 0, MD,ND, pGrid, qGrid);

        //C12UV Descriptor
        MBUV = lts;
        NBUV = 2 * maxRank;
        MUV = nZmiss;
        NUV = 2 * MUV / lts * maxRank;
        EXAGEOSTAT_ALLOCATE_APPROX_MATRIX_TILE(&pDescriptorC12UV, isOOC, nullptr, (HICMA_enum) floatPoint, MBUV, NBUV, MBUV * NBUV, MBUV, NBUV, 0,0, MBUV, NBUV, pGrid, qGrid);

        //C12Ark Descriptor
        MBrk = 1;
        NBrk = 1;
        Mrk = pDescriptorC12UV->mt;
        Nrk = pDescriptorC12UV->mt;
        EXAGEOSTAT_ALLOCATE_APPROX_MATRIX_TILE(&pDescriptorC12rk, isOOC, nullptr, (HICMA_enum) floatPoint, MBrk, NBrk, MBrk * NBrk, Mrk, Nrk, 0, 0,Mrk, Nrk, pGrid, qGrid);

        //C11AD Descriptor
        MBD = lts;
        NBD = lts;
        MD = nZobs;
        ND = MBD;
        EXAGEOSTAT_ALLOCATE_APPROX_MATRIX_TILE(&pDescriptorC22D, isOOC, nullptr, (HICMA_enum) floatPoint, MBD, NBD, MBD * NBD, MD, ND, 0, 0, MD, ND, pGrid, qGrid);

        //C12UV Descriptor
        MBUV = lts;
        NBUV = 2 * maxRank;
        MUV = nZobs;
        NUV = 2 * MUV / lts * maxRank;
        EXAGEOSTAT_ALLOCATE_APPROX_MATRIX_TILE(&pDescriptorC22UV, isOOC, nullptr, (HICMA_enum) floatPoint, MBUV, NBUV, MBUV * NBUV, MBUV, NBUV, 0,0, MBUV, NBUV, pGrid, qGrid);

        //C12Ark Descriptor
        MBrk = 1;
        NBrk = 1;
        Mrk = pDescriptorC22UV->mt;
        Nrk = pDescriptorC22UV->mt;
        EXAGEOSTAT_ALLOCATE_APPROX_MATRIX_TILE(&pDescriptorC22rk, isOOC, nullptr, (HICMA_enum) floatPoint, MBrk, NBrk, MBrk * NBrk, Mrk, Nrk, 0, 0,Mrk, Nrk, pGrid, qGrid);

        //Other descriptors
        // I changed this chameleon to Hicma
        EXAGEOSTAT_ALLOCATE_APPROX_MATRIX_TILE(&pDescriptorZObservations, isOOC, nullptr, (HICMA_enum) floatPoint, lts, lts, lts * lts, nZmiss, 1, 0, 0, nZmiss, 1, pGrid, qGrid);
//        EXAGEOSTAT_ALLOCATE_DENSE_MATRIX_TILE(&CHAM_descmse, &data->mserror, ChamRealDouble, lts, lts, lts * lts, 1, 1, 0, 0, 1, 1, p_grid, q_grid);
    }


//    data->descmse = CHAM_descmse;
//    data->descZactual = CHAMELEON_descZactual;
//    data->descZobs = CHAMELEON_descZobs;


    //stop gsl error handler
    gsl_set_error_handler_off();
}