
// Copyright (c) 2017-2023 King Abdullah University of Science and Technology,
// Copyright (C) 2023 by Brightskies inc,
// All rights reserved.
// ExaGeoStat is a software package, provided by King Abdullah University of Science and Technology (KAUST).

/**
 * @file TestChameleonImplmentationDense.cpp
 *
 * @version 1.0.0
 * @author Sameh Abdulah
 * @date 2023-04-06
**/
#include <libraries/catch/catch.hpp>
#include <cmath>
#include <vector>
#include <linear-algebra-solvers/LinearAlgebraFactory.hpp>
#include <configurations/data-generation/concrete/SyntheticDataConfigurations.hpp>
#include <control/context.h>

using namespace exageostat::linearAlgebra::dense;
using namespace exageostat::linearAlgebra;
using namespace exageostat::common;
using namespace exageostat::configurations::data_configurations;
using namespace std;

void INIT_HARDWARE() {
    ChameleonImplementationDense<double> chameleonImpl;
    chameleonImpl.ExaGeoStatInitContext(4, 0);
    CHAM_context_t *chameleonContext = chameleon_context_self();
    REQUIRE(chameleonContext != nullptr);
}
void FINALIZE_HARDWARE() {
    ChameleonImplementationDense<double> chameleonImpl;
    chameleonImpl.ExaGeoStatInitContext(4, 0);
    chameleonImpl.ExaGeoStatFinalizeContext();
    REQUIRE(chameleon_context_self() == nullptr);
}

// Test that the function initializes all the required descriptors without errors.
void TEST_INITIALIZETION() {
    auto *syntheticDataConfigurations = new SyntheticDataConfigurations();

    SECTION("Single") {
        auto linearAlgebraSolver = LinearAlgebraFactory<float>::CreateLinearAlgebraSolver(EXACT_DENSE);

        syntheticDataConfigurations->SetProblemSize(1000);
        syntheticDataConfigurations->SetDenseTileSize(64);
        linearAlgebraSolver->SetConfigurations(syntheticDataConfigurations);

        REQUIRE(syntheticDataConfigurations->GetDescriptorZcpy() == nullptr);
        REQUIRE(syntheticDataConfigurations->GetDescriptorDeterminant() == nullptr);

        linearAlgebraSolver->InitiateDescriptors();

        REQUIRE(syntheticDataConfigurations->GetDescriptorC().size() == 2);
        REQUIRE(syntheticDataConfigurations->GetDescriptorZ().size() == 1);
        REQUIRE(syntheticDataConfigurations->GetDescriptorProduct().size() == 1);

        REQUIRE(syntheticDataConfigurations->GetDescriptorC()[0] != nullptr);
        REQUIRE(syntheticDataConfigurations->GetDescriptorZ()[0] != nullptr);
        REQUIRE(syntheticDataConfigurations->GetDescriptorProduct()[0] != nullptr);
        REQUIRE(syntheticDataConfigurations->GetDescriptorZcpy() != nullptr);
        REQUIRE(syntheticDataConfigurations->GetDescriptorDeterminant() != nullptr);
    }

    SECTION("Double") {
        auto linearAlgebraSolver = LinearAlgebraFactory<double>::CreateLinearAlgebraSolver(EXACT_DENSE);

        syntheticDataConfigurations->SetProblemSize(1024);
        syntheticDataConfigurations->SetDenseTileSize(64);
        linearAlgebraSolver->SetConfigurations(syntheticDataConfigurations);

        REQUIRE(syntheticDataConfigurations->GetDescriptorZcpy() == nullptr);
        REQUIRE(syntheticDataConfigurations->GetDescriptorDeterminant() == nullptr);

        linearAlgebraSolver->InitiateDescriptors();

        REQUIRE(syntheticDataConfigurations->GetDescriptorC().size() == 4);
        REQUIRE(syntheticDataConfigurations->GetDescriptorZ().size() == 3);
        REQUIRE(syntheticDataConfigurations->GetDescriptorProduct().size() == 3);

        for (auto &descriptorC: syntheticDataConfigurations->GetDescriptorC()) {
            REQUIRE(descriptorC != nullptr);
        }
        for (auto &i: syntheticDataConfigurations->GetDescriptorZ()) {
            REQUIRE(i != nullptr);
        }
        for (auto &i: syntheticDataConfigurations->GetDescriptorProduct()) {
            REQUIRE(i != nullptr);
        }
        REQUIRE(syntheticDataConfigurations->GetDescriptorZcpy() != nullptr);
        REQUIRE(syntheticDataConfigurations->GetDescriptorDeterminant() != nullptr);
    }
    free(syntheticDataConfigurations);
}

//Test that the function initializes the CHAM_descriptorC descriptor correctly.
void TEST_CHAMELEON_DESCRIPTORS_VALUES() {

    SECTION("SINGLE") {
        auto *syntheticDataConfigurations = new SyntheticDataConfigurations();
        auto linearAlgebraSolver = LinearAlgebraFactory<float>::CreateLinearAlgebraSolver(EXACT_DENSE);

        syntheticDataConfigurations->SetProblemSize(6656);
        syntheticDataConfigurations->SetDenseTileSize(512);
        linearAlgebraSolver->SetConfigurations(syntheticDataConfigurations);

        linearAlgebraSolver->InitiateDescriptors();
        auto *CHAM_descriptorC = (CHAM_desc_t *) syntheticDataConfigurations->GetDescriptorC()[0];
        auto *CHAM_descriptorZ = (CHAM_desc_t *) syntheticDataConfigurations->GetDescriptorZ()[0];
        auto *CHAM_descriptorZcpy = (CHAM_desc_t *) syntheticDataConfigurations->GetDescriptorZcpy();
        auto *CHAM_descriptorDeterminant = (CHAM_desc_t *) syntheticDataConfigurations->GetDescriptorDeterminant();
        auto *CHAM_descriptorProduct = (CHAM_desc_t *) syntheticDataConfigurations->GetDescriptorProduct()[0];

        int N = syntheticDataConfigurations->GetProblemSize() * syntheticDataConfigurations->GetP();
        int dts = syntheticDataConfigurations->GetDenseTileSize();
        int pGrid = syntheticDataConfigurations->GetPGrid();
        int qGrid = syntheticDataConfigurations->GetQGrid();

        REQUIRE(CHAM_descriptorC->m == N);
        REQUIRE(CHAM_descriptorZ->m == N);
        REQUIRE(CHAM_descriptorZcpy->m == N);
        REQUIRE(CHAM_descriptorDeterminant->m == 1);
        REQUIRE(CHAM_descriptorProduct->m == 1);

        REQUIRE(CHAM_descriptorC->n == N);
        REQUIRE(CHAM_descriptorZ->n == 1);
        REQUIRE(CHAM_descriptorZcpy->n == 1);
        REQUIRE(CHAM_descriptorDeterminant->n == 1);
        REQUIRE(CHAM_descriptorProduct->n == 1);

        REQUIRE(CHAM_descriptorC->mb == dts);
        REQUIRE(CHAM_descriptorZ->mb == dts);
        REQUIRE(CHAM_descriptorZcpy->mb == dts);
        REQUIRE(CHAM_descriptorDeterminant->mb == dts);
        REQUIRE(CHAM_descriptorProduct->mb == dts);

        REQUIRE(CHAM_descriptorC->nb == dts);
        REQUIRE(CHAM_descriptorZ->nb == dts);
        REQUIRE(CHAM_descriptorZcpy->nb == dts);
        REQUIRE(CHAM_descriptorDeterminant->nb == dts);
        REQUIRE(CHAM_descriptorProduct->nb == dts);

        REQUIRE(CHAM_descriptorC->bsiz == dts * dts);
        REQUIRE(CHAM_descriptorZ->bsiz == dts * dts);
        REQUIRE(CHAM_descriptorZcpy->bsiz == dts * dts);
        REQUIRE(CHAM_descriptorDeterminant->bsiz == dts * dts);
        REQUIRE(CHAM_descriptorProduct->bsiz == dts * dts);

        REQUIRE(CHAM_descriptorC->i == 0);
        REQUIRE(CHAM_descriptorZ->i == 0);
        REQUIRE(CHAM_descriptorZcpy->i == 0);
        REQUIRE(CHAM_descriptorDeterminant->i == 0);
        REQUIRE(CHAM_descriptorProduct->i == 0);

        REQUIRE(CHAM_descriptorC->j == 0);
        REQUIRE(CHAM_descriptorZ->j == 0);
        REQUIRE(CHAM_descriptorZcpy->j == 0);
        REQUIRE(CHAM_descriptorDeterminant->j == 0);
        REQUIRE(CHAM_descriptorProduct->j == 0);

        REQUIRE(CHAM_descriptorC->mt == ceil((N * 1.0) / (dts * 1.0)));
        REQUIRE(CHAM_descriptorZ->mt == ceil((N * 1.0) / (dts * 1.0)));
        REQUIRE(CHAM_descriptorZcpy->mt == ceil((N * 1.0) / (dts * 1.0)));
        REQUIRE(CHAM_descriptorDeterminant->mt == 1);
        REQUIRE(CHAM_descriptorProduct->mt == 1);

        REQUIRE(CHAM_descriptorC->nt == ceil((N * 1.0) / (dts * 1.0)));
        REQUIRE(CHAM_descriptorZ->nt == 1);
        REQUIRE(CHAM_descriptorZcpy->nt == 1);
        REQUIRE(CHAM_descriptorDeterminant->nt == 1);
        REQUIRE(CHAM_descriptorProduct->nt == 1);

        REQUIRE(CHAM_descriptorC->lm == N);
        REQUIRE(CHAM_descriptorZ->lm == N);
        REQUIRE(CHAM_descriptorZcpy->lm == N);
        REQUIRE(CHAM_descriptorDeterminant->lm == 1);
        REQUIRE(CHAM_descriptorProduct->lm == 1);

        REQUIRE(CHAM_descriptorC->ln == N);
        REQUIRE(CHAM_descriptorZ->ln == 1);
        REQUIRE(CHAM_descriptorZcpy->ln == 1);
        REQUIRE(CHAM_descriptorDeterminant->ln == 1);
        REQUIRE(CHAM_descriptorProduct->ln == 1);

        REQUIRE(CHAM_descriptorC->p == pGrid);
        REQUIRE(CHAM_descriptorZ->p == pGrid);
        REQUIRE(CHAM_descriptorZcpy->p == pGrid);
        REQUIRE(CHAM_descriptorDeterminant->p == pGrid);
        REQUIRE(CHAM_descriptorProduct->p == pGrid);

        REQUIRE(CHAM_descriptorC->q == qGrid);
        REQUIRE(CHAM_descriptorZ->q == qGrid);
        REQUIRE(CHAM_descriptorZcpy->q == qGrid);
        REQUIRE(CHAM_descriptorDeterminant->q == qGrid);
        REQUIRE(CHAM_descriptorProduct->q == qGrid);


        auto *mat = (float *) CHAM_descriptorC->mat;
        for (auto i = 0;
             i < (CHAM_descriptorC->mt - 1) * (CHAM_descriptorC->nt - 1) * (CHAM_descriptorC->bsiz - 1); i++) {
                REQUIRE(mat[i] == 0.0f);
        }

        mat = (float *) CHAM_descriptorZ->mat;
        auto *matZcpy = (float *) CHAM_descriptorZcpy->mat;
        for (auto i = 0;
             i < (CHAM_descriptorZ->mt - 1) * (CHAM_descriptorZ->nt - 1) * (CHAM_descriptorZ->bsiz - 1); i++) {
            REQUIRE(mat[i] == 0.0f);
            REQUIRE(matZcpy[i] == 0.0f);
        }

        mat = (float *) CHAM_descriptorDeterminant->mat;
        auto *matProduct = (float *) CHAM_descriptorProduct->mat;
        for (auto i = 0; i < (CHAM_descriptorDeterminant->mt - 1) * (CHAM_descriptorDeterminant->nt - 1) *
                             (CHAM_descriptorDeterminant->bsiz - 1); i++) {
            REQUIRE(mat[i] == 0.0f);
            REQUIRE(matProduct[i] == 0.0f);
        }
        free(syntheticDataConfigurations);
    }

    SECTION("DOUBLE") {
        auto *syntheticDataConfigurations = new SyntheticDataConfigurations();
        auto linearAlgebraSolver = LinearAlgebraFactory<double>::CreateLinearAlgebraSolver(EXACT_DENSE);

        syntheticDataConfigurations->SetProblemSize(6400);
        syntheticDataConfigurations->SetDenseTileSize(512);
        linearAlgebraSolver->SetConfigurations(syntheticDataConfigurations);

        linearAlgebraSolver->InitiateDescriptors();
        auto *CHAM_descriptorC = (CHAM_desc_t *) syntheticDataConfigurations->GetDescriptorC()[0];
        auto *CHAM_descsubC11 = (CHAM_desc_t *) syntheticDataConfigurations->GetDescriptorC()[1];
        auto *CHAM_descsubC12 = (CHAM_desc_t *) syntheticDataConfigurations->GetDescriptorC()[2];
        auto *CHAM_descsubC22 = (CHAM_desc_t *) syntheticDataConfigurations->GetDescriptorC()[3];
        auto *CHAM_descriptorZ = (CHAM_desc_t *) syntheticDataConfigurations->GetDescriptorZ()[0];
        auto *CHAM_descriptorZcpy = (CHAM_desc_t *) syntheticDataConfigurations->GetDescriptorZcpy();
        auto *CHAM_descriptorDeterminant = (CHAM_desc_t *) syntheticDataConfigurations->GetDescriptorDeterminant();
        vector<void *> &pDescriptorProduct = syntheticDataConfigurations->GetDescriptorProduct();
        vector<void *> &pDescriptorZ = syntheticDataConfigurations->GetDescriptorZ();


        int N = syntheticDataConfigurations->GetProblemSize() * syntheticDataConfigurations->GetP();
        int dts = syntheticDataConfigurations->GetDenseTileSize();
        int pGrid = syntheticDataConfigurations->GetPGrid();
        int qGrid = syntheticDataConfigurations->GetQGrid();

        REQUIRE(CHAM_descriptorC->m == N);
        REQUIRE(CHAM_descriptorZ->m == N);
        REQUIRE(CHAM_descriptorZcpy->m == N);
        REQUIRE(CHAM_descriptorDeterminant->m == 1);
        REQUIRE(CHAM_descsubC11->m == N / 2);
        REQUIRE(CHAM_descsubC12->m == N / 2);
        REQUIRE(CHAM_descsubC22->m == N / 2);
        for (auto & idx : pDescriptorProduct) {
            auto **CHAM_descriptorProduct = (CHAM_desc_t **) &idx;
            REQUIRE((*CHAM_descriptorProduct)->m == 1);
        }
        for (int idx = 1; idx < pDescriptorZ.size(); idx++) {
            auto **CHAM_descriptorZ_ = (CHAM_desc_t **) &pDescriptorZ[idx];
            REQUIRE((*CHAM_descriptorZ_)->m == N / 2);
        }

        REQUIRE(CHAM_descriptorC->n == N);
        REQUIRE(CHAM_descriptorZ->n == 1);
        REQUIRE(CHAM_descriptorZcpy->n == 1);
        REQUIRE(CHAM_descriptorDeterminant->n == 1);
        REQUIRE(CHAM_descsubC11->n == N / 2);
        REQUIRE(CHAM_descsubC12->n == N / 2);
        REQUIRE(CHAM_descsubC22->n == N / 2);
        for (auto & idx : pDescriptorProduct) {
            auto **CHAM_descriptorProduct = (CHAM_desc_t **) &idx;
            REQUIRE((*CHAM_descriptorProduct)->n == 1);
        }
        for (int idx = 1; idx < pDescriptorZ.size(); idx++) {
            auto **CHAM_descriptorZ_ = (CHAM_desc_t **) &pDescriptorZ[idx];
            REQUIRE((*CHAM_descriptorZ_)->n == 1);
        }

        REQUIRE(CHAM_descriptorC->mb == dts);
        REQUIRE(CHAM_descriptorZ->mb == dts);
        REQUIRE(CHAM_descriptorZcpy->mb == dts);
        REQUIRE(CHAM_descriptorDeterminant->mb == dts);
        REQUIRE(CHAM_descsubC11->mb == dts);
        REQUIRE(CHAM_descsubC12->mb == dts);
        REQUIRE(CHAM_descsubC22->mb == dts);
        for (auto & idx : pDescriptorProduct) {
            auto **CHAM_descriptorProduct = (CHAM_desc_t **) &idx;
            REQUIRE((*CHAM_descriptorProduct)->mb == dts);
        }
        for (int idx = 1; idx < pDescriptorZ.size(); idx++) {
            auto **CHAM_descriptorZ_ = (CHAM_desc_t **) &pDescriptorZ[idx];
            REQUIRE((*CHAM_descriptorZ_)->mb == dts);
        }

        REQUIRE(CHAM_descriptorC->nb == dts);
        REQUIRE(CHAM_descriptorZ->nb == dts);
        REQUIRE(CHAM_descriptorZcpy->nb == dts);
        REQUIRE(CHAM_descriptorDeterminant->nb == dts);
        REQUIRE(CHAM_descsubC11->nb == dts);
        REQUIRE(CHAM_descsubC12->nb == dts);
        REQUIRE(CHAM_descsubC22->nb == dts);
        for (auto & idx : pDescriptorProduct) {
            auto **CHAM_descriptorProduct = (CHAM_desc_t **) &idx;
            REQUIRE((*CHAM_descriptorProduct)->nb == dts);
        }
        for (int idx = 1; idx < pDescriptorZ.size(); idx++) {
            auto **CHAM_descriptorZ_ = (CHAM_desc_t **) &pDescriptorZ[idx];
            REQUIRE((*CHAM_descriptorZ_)->nb == dts);
        }

        REQUIRE(CHAM_descriptorC->bsiz == dts * dts);
        REQUIRE(CHAM_descriptorZ->bsiz == dts * dts);
        REQUIRE(CHAM_descriptorZcpy->bsiz == dts * dts);
        REQUIRE(CHAM_descriptorDeterminant->bsiz == dts * dts);
        REQUIRE(CHAM_descsubC11->bsiz == dts * dts);
        REQUIRE(CHAM_descsubC12->bsiz == dts * dts);
        REQUIRE(CHAM_descsubC22->bsiz == dts * dts);
        for (auto & idx : pDescriptorProduct) {
            auto **CHAM_descriptorProduct = (CHAM_desc_t **) &idx;
            REQUIRE((*CHAM_descriptorProduct)->bsiz == dts * dts);
        }
        for (int idx = 1; idx < pDescriptorZ.size(); idx++) {
            auto **CHAM_descriptorZ_ = (CHAM_desc_t **) &pDescriptorZ[idx];
            REQUIRE((*CHAM_descriptorZ_)->bsiz == dts * dts);
        }

        REQUIRE(CHAM_descriptorC->i == 0);
        REQUIRE(CHAM_descriptorZ->i == 0);
        REQUIRE(CHAM_descriptorZcpy->i == 0);
        REQUIRE(CHAM_descriptorDeterminant->i == 0);
        REQUIRE(CHAM_descsubC11->i == 0);
        REQUIRE(CHAM_descsubC12->i == N / 2);
        REQUIRE(CHAM_descsubC22->i == N / 2);
        for (auto & idx : pDescriptorProduct) {
            auto **CHAM_descriptorProduct = (CHAM_desc_t **) &idx;
            REQUIRE((*CHAM_descriptorProduct)->i == 0);
        }
        for (int idx = 1; idx < pDescriptorZ.size(); idx++) {
            auto **CHAM_descriptorZ_ = (CHAM_desc_t **) &pDescriptorZ[idx];
            REQUIRE((*CHAM_descriptorZ_)->i == 0);
        }

        REQUIRE(CHAM_descriptorC->j == 0);
        REQUIRE(CHAM_descriptorZ->j == 0);
        REQUIRE(CHAM_descriptorZcpy->j == 0);
        REQUIRE(CHAM_descriptorDeterminant->j == 0);
        REQUIRE(CHAM_descsubC11->j == 0);
        REQUIRE(CHAM_descsubC12->j == 0);
        REQUIRE(CHAM_descsubC22->j == N / 2);
        for (auto & idx : pDescriptorProduct) {
            auto **CHAM_descriptorProduct = (CHAM_desc_t **) &idx;
            REQUIRE((*CHAM_descriptorProduct)->j == 0);
        }
        for (int idx = 1; idx < pDescriptorZ.size(); idx++) {
            auto **CHAM_descriptorZ_ = (CHAM_desc_t **) &pDescriptorZ[idx];
            REQUIRE((*CHAM_descriptorZ_)->j == 0);
        }

        REQUIRE(CHAM_descriptorC->mt == ceil((N * 1.0) / (dts * 1.0)));
        REQUIRE(CHAM_descriptorZ->mt == ceil((N * 1.0) / (dts * 1.0)));
        REQUIRE(CHAM_descriptorZcpy->mt == ceil((N * 1.0) / (dts * 1.0)));
        REQUIRE(CHAM_descriptorDeterminant->mt == 1);
        REQUIRE(CHAM_descsubC11->mt == ceil((N / 2 * 1.0) / (dts * 1.0)));
        REQUIRE(CHAM_descsubC12->mt == ceil((N / 2 * 1.0) / (dts * 1.0)));
        REQUIRE(CHAM_descsubC22->mt == ceil((N / 2 * 1.0) / (dts * 1.0)));
        for (auto & idx : pDescriptorProduct) {
            auto **CHAM_descriptorProduct = (CHAM_desc_t **) &idx;
            REQUIRE((*CHAM_descriptorProduct)->mt == 1);
        }
        for (int idx = 1; idx < pDescriptorZ.size(); idx++) {
            auto **CHAM_descriptorZ_ = (CHAM_desc_t **) &pDescriptorZ[idx];
            REQUIRE((*CHAM_descriptorZ_)->mt == ceil((N / 2 * 1.0) / (dts * 1.0)));
        }

        REQUIRE(CHAM_descriptorC->nt == ceil((N * 1.0) / (dts * 1.0)));
        REQUIRE(CHAM_descriptorZ->nt == 1);
        REQUIRE(CHAM_descriptorZcpy->nt == 1);
        REQUIRE(CHAM_descriptorDeterminant->nt == 1);
        REQUIRE(CHAM_descsubC11->nt == ceil((N / 2 * 1.0) / (dts * 1.0)));
        REQUIRE(CHAM_descsubC12->nt == ceil((N / 2 * 1.0) / (dts * 1.0)));
        REQUIRE(CHAM_descsubC22->nt == ceil((N / 2 * 1.0) / (dts * 1.0)));
        for (auto & idx : pDescriptorProduct) {
            auto **CHAM_descriptorProduct = (CHAM_desc_t **) &idx;
            REQUIRE((*CHAM_descriptorProduct)->nt == 1);
        }
        for (int idx = 1; idx < pDescriptorZ.size(); idx++) {
            auto **CHAM_descriptorZ_ = (CHAM_desc_t **) &pDescriptorZ[idx];
            REQUIRE((*CHAM_descriptorZ_)->nt == 1);
        }

        REQUIRE(CHAM_descriptorC->lm == N);
        REQUIRE(CHAM_descriptorZ->lm == N);
        REQUIRE(CHAM_descriptorZcpy->lm == N);
        REQUIRE(CHAM_descriptorDeterminant->lm == 1);
        REQUIRE(CHAM_descsubC11->lm == N);
        REQUIRE(CHAM_descsubC12->lm == N);
        REQUIRE(CHAM_descsubC22->lm == N);
        for (auto & idx : pDescriptorProduct) {
            auto **CHAM_descriptorProduct = (CHAM_desc_t **) &idx;
            REQUIRE((*CHAM_descriptorProduct)->lm == 1);
        }
        for (int idx = 1; idx < pDescriptorZ.size(); idx++) {
            auto **CHAM_descriptorZ_ = (CHAM_desc_t **) &pDescriptorZ[idx];
            REQUIRE((*CHAM_descriptorZ_)->lm == N / 2);
        }

        REQUIRE(CHAM_descriptorC->ln == N);
        REQUIRE(CHAM_descriptorZ->ln == 1);
        REQUIRE(CHAM_descriptorZcpy->ln == 1);
        REQUIRE(CHAM_descriptorDeterminant->ln == 1);
        REQUIRE(CHAM_descsubC11->ln == N);
        REQUIRE(CHAM_descsubC12->ln == N);
        REQUIRE(CHAM_descsubC22->ln == N);
        for (auto & idx : pDescriptorProduct) {
            auto **CHAM_descriptorProduct = (CHAM_desc_t **) &idx;
            REQUIRE((*CHAM_descriptorProduct)->ln == 1);
        }
        for (int idx = 1; idx < pDescriptorZ.size(); idx++) {
            auto **CHAM_descriptorZ_ = (CHAM_desc_t **) &pDescriptorZ[idx];
            REQUIRE((*CHAM_descriptorZ_)->ln == 1);
        }

        REQUIRE(CHAM_descriptorC->p == pGrid);
        REQUIRE(CHAM_descriptorZ->p == pGrid);
        REQUIRE(CHAM_descriptorZcpy->p == pGrid);
        REQUIRE(CHAM_descriptorDeterminant->p == pGrid);
        REQUIRE(CHAM_descsubC11->p == pGrid);
        REQUIRE(CHAM_descsubC12->p == pGrid);
        REQUIRE(CHAM_descsubC22->p == pGrid);
        for (auto & idx : pDescriptorProduct) {
            auto **CHAM_descriptorProduct = (CHAM_desc_t **) &idx;
            REQUIRE((*CHAM_descriptorProduct)->p == pGrid);
        }
        for (int idx = 1; idx < pDescriptorZ.size(); idx++) {
            auto **CHAM_descriptorZ_ = (CHAM_desc_t **) &pDescriptorZ[idx];
            REQUIRE((*CHAM_descriptorZ_)->p == pGrid);
        }

        REQUIRE(CHAM_descriptorC->q == qGrid);
        REQUIRE(CHAM_descriptorZ->q == qGrid);
        REQUIRE(CHAM_descriptorZcpy->q == qGrid);
        REQUIRE(CHAM_descriptorDeterminant->q == qGrid);
        REQUIRE(CHAM_descsubC11->q == qGrid);
        REQUIRE(CHAM_descsubC12->q == qGrid);
        REQUIRE(CHAM_descsubC22->q == qGrid);
        for (auto & idx : pDescriptorProduct) {
            auto **CHAM_descriptorProduct = (CHAM_desc_t **) &idx;
            REQUIRE((*CHAM_descriptorProduct)->q == qGrid);
        }
        for (int idx = 1; idx < pDescriptorZ.size(); idx++) {
            auto **CHAM_descriptorZ_ = (CHAM_desc_t **) &pDescriptorZ[idx];
            REQUIRE((*CHAM_descriptorZ_)->q == qGrid);
        }


        auto *mat = (double *) CHAM_descriptorC->mat;
        for (auto i = 0;
             i < (CHAM_descriptorC->mt - 1) * (CHAM_descriptorC->nt - 1) * (CHAM_descriptorC->bsiz - 1); i++) {
            REQUIRE(mat[i] == 0.0f);
        }
        auto *matC11 = (double *) CHAM_descsubC11->mat;
        auto *matC12 = (double *) CHAM_descsubC12->mat;
        auto *matC22 = (double *) CHAM_descsubC22->mat;
        for (auto i = 0; i < (CHAM_descsubC11->mt - 1) * (CHAM_descsubC11->nt - 1) * (CHAM_descsubC11->bsiz - 1); i++) {
            REQUIRE(matC11[i] == 0.0f);
            REQUIRE(matC12[i] == 0.0f);
            REQUIRE(matC22[i] == 0.0f);
        }

        mat = (double *) CHAM_descriptorZ->mat;
        auto *matZcpy = (double *) CHAM_descriptorZcpy->mat;
        for (auto i = 0;
             i < (CHAM_descriptorZ->mt - 1) * (CHAM_descriptorZ->nt - 1) * (CHAM_descriptorZ->bsiz - 1); i++) {
            REQUIRE(mat[i] == 0.0f);
            REQUIRE(matZcpy[i] == 0.0f);
        }
        for (auto & idx : pDescriptorZ) {
            auto **CHAM_descriptorZ_ = (CHAM_desc_t **) &idx;
            mat = (double *) (*CHAM_descriptorZ_)->mat;
            for (auto i = 0; i < ((*CHAM_descriptorZ_)->mt - 1) * ((*CHAM_descriptorZ_)->nt - 1) *
                                 ((*CHAM_descriptorZ_)->bsiz - 1); i++) {
                REQUIRE(mat[i] == 0.0f);
            }
        }

        mat = (double *) CHAM_descriptorDeterminant->mat;
        for (auto i = 0; i < (CHAM_descriptorDeterminant->mt - 1) * (CHAM_descriptorDeterminant->nt - 1) *
                             (CHAM_descriptorDeterminant->bsiz - 1); i++) {
            REQUIRE(mat[i] == 0.0f);
            for (auto & idx : pDescriptorProduct) {
                auto **CHAM_descriptorProduct = (CHAM_desc_t **) &idx;
                auto *matProduct = (double *) (*CHAM_descriptorProduct)->mat;
                REQUIRE(matProduct[i] == 0.0f);
            }
        }
        free(syntheticDataConfigurations);
    }

}
void blabla(){
    auto *syntheticDataConfigurations = new SyntheticDataConfigurations();
    auto linearAlgebraSolver = LinearAlgebraFactory<double>::CreateLinearAlgebraSolver(EXACT_DENSE);

    syntheticDataConfigurations->SetProblemSize(6400);
    syntheticDataConfigurations->SetDenseTileSize(512);
    linearAlgebraSolver->SetConfigurations(syntheticDataConfigurations);

    linearAlgebraSolver->InitiateDescriptors();
    linearAlgebraSolver->testKernelfornow();
}
TEST_CASE("Chameleon Implementation Dense") {
    INIT_HARDWARE();
    FINALIZE_HARDWARE();
    TEST_INITIALIZETION();
    TEST_CHAMELEON_DESCRIPTORS_VALUES();
//blabla();
}
