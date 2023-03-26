
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
using namespace std;

template<typename T>
void HicmaImplementation<T>::InitiateDescriptors() {

    vector<void *> pDescriptorC =  this->mpConfigurations->GetDescriptorC();
    vector<void *> pDescriptorZ = this->mpConfigurations->GetDescriptorZ();
    auto* pDescriptorZcpy = (HICMA_desc_t*) this->mpConfigurations->GetDescriptorZcpy();
    vector<void *> pDescriptorProduct = this->mpConfigurations->GetDescriptorProduct();
    auto* pDescriptorDeterminant = (HICMA_desc_t*) this->mpConfigurations->GetDescriptorDeterminant();


    int vectorSize = 0;
    HICMA_sequence_t *pSequence;
    HICMA_request_t request[2] = {HICMA_SUCCESS, HICMA_SUCCESS};

    int N = this->mpConfigurations->GetProblemSize() * this->mpConfigurations->GetP();
    int dts = this->mpConfigurations->GetTileSize();
    int pGrid = this->mpConfigurations->GetPGrid();
    int qGrid = this->mpConfigurations->GetQGrid();

    // For distributed system and should be removed
    T *Zcpy = (T *) malloc(N * sizeof(T));

    //// TODO:: what is this?
    //CDense Descriptor
//    if (data->check == 1) {
//        MBC = lts;
//        NBC = lts;
//        MC = N;
//        NC = N;
//    } else {
//        MBC = 1;
//        NBC = 1;
//        MC = lts;
//        NC = lts;
//    }

    T dotProductValue;

    //Identifies a set of routines sharing common exception handling.
//    CHAMELEON_Sequence_Create(&pSequence);
//
//    FloatPoint floatPoint;
//    if (sizeof(T) == SIZE_OF_FLOAT) {
//        floatPoint = EXAGEOSTAT_REAL_FLOAT;
//        vectorSize = 1;
//    } else {
//        floatPoint = EXAGEOSTAT_REAL_DOUBLE;
//        vectorSize = 3;
//    }
//    // Depending on the passed Precession, A for loop with value 1 or 3 will initialize the vectors
//    for (int idx = 0; idx < vectorSize; idx++){
//        pDescriptorC.push_back(nullptr);
//        pDescriptorZ.push_back(nullptr);
//        pDescriptorProduct.push_back(nullptr);
//    }
//
//    CHAM_desc_t* CHAM_descriptorC = (CHAM_desc_t*) pDescriptorC[0];
//    if(vectorSize > 1){
//        pDescriptorC.push_back(nullptr);
//    }
//
//    //// TODO: change tile size one for hicma and one for chameleon
//    EXAGEOSTAT_ALLOCATE_MATRIX_TILE(&CHAM_descriptorC, nullptr, (cham_flttype_t) floatPoint, dts, dts, dts * dts, N, N, 0, 0, N, N, pGrid, qGrid);
//    CHAM_desc_t* CHAM_descriptorZ = (CHAM_desc_t*) pDescriptorZ[0];
//    EXAGEOSTAT_ALLOCATE_MATRIX_TILE(&CHAM_descriptorZ, nullptr, (cham_flttype_t) floatPoint, dts, dts, dts * dts, N, 1, 0,0, N, 1, pGrid, qGrid);
//    EXAGEOSTAT_ALLOCATE_MATRIX_TILE(&pDescriptorZcpy, Zcpy, (cham_flttype_t) floatPoint, dts, dts, dts * dts, N, 1, 0, 0, N, 1, pGrid, qGrid);
//    EXAGEOSTAT_ALLOCATE_MATRIX_TILE(&pDescriptorDeterminant, &dotProductValue, (cham_flttype_t) floatPoint, dts, dts, dts * dts, 1, 1, 0, 0, 1, 1, pGrid, qGrid);
//
//    for (int idx = 1; idx < pDescriptorZ.size(); idx++) {
//        CHAM_desc_t* CHAM_descriptorZ_ = (CHAM_desc_t*) pDescriptorZ[idx];
//        EXAGEOSTAT_ALLOCATE_MATRIX_TILE(&CHAM_descriptorZ_, nullptr, (cham_flttype_t) floatPoint, dts, dts, dts * dts, N / 2, 1, 0, 0, N / 2, 1, pGrid, qGrid);
//    }
//
//    for (int idx = 0; idx < pDescriptorZ.size(); idx++) {
//        CHAM_desc_t* CHAM_descriptorProduct = (CHAM_desc_t*) pDescriptorProduct[idx];
//        EXAGEOSTAT_ALLOCATE_MATRIX_TILE(&CHAM_descriptorProduct, &dotProductValue, (cham_flttype_t) floatPoint, dts, dts, dts * dts, 1, 1, 0, 0, 1, 1, pGrid, qGrid)
//    }
//
    //stop gsl error handler
    gsl_set_error_handler_off();
}