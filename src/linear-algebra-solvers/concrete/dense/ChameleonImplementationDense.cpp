
// Copyright (c) 2017-2023 King Abdullah University of Science and Technology,
// Copyright (C) 2023 by Brightskies inc,
// All rights reserved.
// ExaGeoStat is a software package, provided by King Abdullah University of Science and Technology (KAUST).

/**
 * @file ChameleonAllocateDescriptors.cpp
 * @brief Sets up the Chameleon descriptors needed for the dense matrix computations in ExaGeoStat.
 * @version 1.0.0
 * @author Sameh Abdulah
 * @date 2023-03-20
**/

#include <linear-algebra-solvers/concrete/dense/ChameleonImplementationDense.hpp>

extern "C" {
#include <chameleon/struct.h>
#include <chameleon.h>
#include <control/descriptor.h>
#include <control/context.h>
}
using namespace exageostat::linearAlgebra::dense;
using namespace exageostat::common;
using namespace std;

template<typename T>
void ChameleonImplementationDense<T>::InitiateDescriptors() {

    this->ExaGeoStatInitContext(this->mpConfigurations->GetCoresNumber(), this->mpConfigurations->GetGPUsNumber());

    vector<void *> &pDescriptorC = this->mpConfigurations->GetDescriptorC();
    vector<void *> &pDescriptorZ = this->mpConfigurations->GetDescriptorZ();
    auto pDescriptorZcpy = (CHAM_desc_t **) &this->mpConfigurations->GetDescriptorZcpy();
    vector<void *> &pDescriptorProduct = this->mpConfigurations->GetDescriptorProduct();
    auto pDescriptorDeterminant = (CHAM_desc_t **) &this->mpConfigurations->GetDescriptorDeterminant();

    int vectorSize;
    RUNTIME_sequence_t *pSequence;

    int N = this->mpConfigurations->GetProblemSize() * this->mpConfigurations->GetP();
    int dts = this->mpConfigurations->GetDenseTileSize();
    int pGrid = this->mpConfigurations->GetPGrid();
    int qGrid = this->mpConfigurations->GetQGrid();
    bool isOOC = this->mpConfigurations->GetIsOOC();

    // For distributed system and should be removed
    T *Zcpy = (T *) malloc(N * sizeof(T));
    T dotProductValue;

    //Identifies a set of routines sharing common exception handling.
    CHAMELEON_Sequence_Create(&pSequence);

    FloatPoint floatPoint;
    if (sizeof(T) == SIZE_OF_FLOAT) {
        floatPoint = EXAGEOSTAT_REAL_FLOAT;
        vectorSize = 1;
    } else {
        floatPoint = EXAGEOSTAT_REAL_DOUBLE;
        vectorSize = 3;
    }

    // Depending on the passed Precession, the descriptor will resize for value 1 or 3.
    pDescriptorC.resize(vectorSize + 1, nullptr);
    pDescriptorZ.resize(vectorSize, nullptr);
    pDescriptorProduct.resize(vectorSize, nullptr);

    auto **CHAM_descriptorC = (CHAM_desc_t **) &pDescriptorC[0];
    EXAGEOSTAT_ALLOCATE_DENSE_MATRIX_TILE(CHAM_descriptorC, isOOC, nullptr, (cham_flttype_t) floatPoint, dts, dts,
                                          dts * dts, N, N, 0, 0, N, N, pGrid, qGrid)

    if (vectorSize > 1) {

        auto **CHAM_descsubC11 = (CHAM_desc_t **) &pDescriptorC[1];
        auto **CHAM_descsubC12 = (CHAM_desc_t **) &pDescriptorC[2];
        auto **CHAM_descsubC22 = (CHAM_desc_t **) &pDescriptorC[3];

        *CHAM_descsubC11 = chameleon_desc_submatrix(*CHAM_descriptorC, 0, 0, (*CHAM_descriptorC)->m / 2,
                                                    (*CHAM_descriptorC)->n / 2);
        *CHAM_descsubC12 = chameleon_desc_submatrix(*CHAM_descriptorC, (*CHAM_descriptorC)->m / 2, 0,
                                                    (*CHAM_descriptorC)->m / 2, (*CHAM_descriptorC)->n / 2);
        *CHAM_descsubC22 = chameleon_desc_submatrix(*CHAM_descriptorC, (*CHAM_descriptorC)->m / 2,
                                                    (*CHAM_descriptorC)->n / 2,
                                                    (*CHAM_descriptorC)->m / 2, (*CHAM_descriptorC)->n / 2);
    }
    EXAGEOSTAT_ALLOCATE_DENSE_MATRIX_TILE(((CHAM_desc_t **) &pDescriptorZ[0]), isOOC, nullptr,
                                          (cham_flttype_t) floatPoint, dts, dts, dts * dts, N, 1, 0, 0, N, 1, pGrid,
                                          qGrid)
    EXAGEOSTAT_ALLOCATE_DENSE_MATRIX_TILE(pDescriptorZcpy, isOOC, Zcpy, (cham_flttype_t) floatPoint, dts, dts,
                                          dts * dts, N, 1, 0, 0, N, 1, pGrid, qGrid)
    EXAGEOSTAT_ALLOCATE_DENSE_MATRIX_TILE(pDescriptorDeterminant, isOOC, &dotProductValue, (cham_flttype_t) floatPoint,
                                          dts, dts, dts * dts, 1, 1, 0, 0, 1, 1, pGrid, qGrid)

    for (int idx = 1; idx < pDescriptorZ.size(); idx++) {
        auto **CHAM_descriptorZ_ = (CHAM_desc_t **) &pDescriptorZ[idx];
        EXAGEOSTAT_ALLOCATE_DENSE_MATRIX_TILE(CHAM_descriptorZ_, isOOC, nullptr, (cham_flttype_t) floatPoint, dts, dts,
                                              dts * dts, N / 2, 1, 0, 0, N / 2, 1, pGrid, qGrid)
    }

    for (auto & idx : pDescriptorProduct) {
        auto **CHAM_descriptorProduct = (CHAM_desc_t **) &idx;
        EXAGEOSTAT_ALLOCATE_DENSE_MATRIX_TILE(CHAM_descriptorProduct, isOOC, &dotProductValue,
                                              (cham_flttype_t) floatPoint, dts, dts, dts * dts, 1, 1, 0, 0, 1, 1, pGrid,
                                              qGrid)
    }
    this->ExaGeoStatFinalizeContext();
    //stop gsl error handler
    gsl_set_error_handler_off();
}

template<typename T>
void ChameleonImplementationDense<T>::ExaGeoStatInitContext(const int &apCoresNumber, const int &apGPUs) {

    CHAM_context_t *chameleonContext;
    chameleonContext = chameleon_context_self();
    if (chameleonContext != nullptr) {
        printf("Another instance of Chameleon is already running...!");
    } else {
        CHAMELEON_user_tag_size(31, 26);
        CHAMELEON_Init(apCoresNumber, apGPUs)
    }
}

template<typename T>
void ChameleonImplementationDense<T>::ExaGeoStatFinalizeContext() {

    CHAM_context_t *chameleonContext;
    chameleonContext = chameleon_context_self();
    if (chameleonContext == nullptr) {
        printf("No active instance oh Chameleon...please use ExaGeoStatInitContext() function to initiate a new instance!\n");
    } else
        CHAMELEON_Finalize();
}