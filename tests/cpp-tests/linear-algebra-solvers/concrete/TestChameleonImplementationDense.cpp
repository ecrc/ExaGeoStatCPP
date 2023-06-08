
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
#include <api/ExaGeoStat.hpp>
#include <data-generators/DataGenerator.hpp>

using namespace exageostat::linearAlgebra::dense;
using namespace exageostat::linearAlgebra;
using namespace exageostat::common;
using namespace exageostat::configurations::data_configurations;
using namespace std;

void INIT_FINALIZE_HARDWARE() {

    ChameleonImplementationDense<float> chameleonImpl;
    // Initialise the context
    chameleonImpl.ExaGeoStatInitContext(4, 0);
    CHAM_context_t *chameleonContext1 = chameleon_context_self();
    REQUIRE(chameleonContext1 != nullptr);

    // Finalize the context.
    chameleonImpl.ExaGeoStatFinalizeContext();
    REQUIRE(chameleon_context_self() == nullptr);

    // Test using operations without initialise Hardware.
    REQUIRE_THROWS_WITH(chameleonImpl.InitiateDescriptors(),
                        "ExaGeoStat hardware is not initialized, please use 'ExaGeoStat<double/float>::ExaGeoStatInitializeHardware(configurations)'.");
}

// Test that the function initializes all the required descriptors without errors.
void TEST_DESCRIPTORS_INITIALIZATION() {
    SyntheticDataConfigurations synthetic_data_configurations;

    SECTION("Single") {

        // Initialise Hardware.
        exageostat::api::ExaGeoStat<float>::ExaGeoStatInitializeHardware(&synthetic_data_configurations);

        auto linearAlgebraSolver = LinearAlgebraFactory<float>::CreateLinearAlgebraSolver(EXACT_DENSE);

        synthetic_data_configurations.SetProblemSize(32);
        synthetic_data_configurations.SetDenseTileSize(32);
        linearAlgebraSolver->SetConfigurations(&synthetic_data_configurations);

        linearAlgebraSolver->InitiateDescriptors();

        REQUIRE(synthetic_data_configurations.GetDescriptorC().size() == 2);
        REQUIRE(synthetic_data_configurations.GetDescriptorZ().size() == 1);
        REQUIRE(synthetic_data_configurations.GetDescriptorProduct().size() == 1);

        REQUIRE(synthetic_data_configurations.GetDescriptorC()[0] != nullptr);
        REQUIRE(synthetic_data_configurations.GetDescriptorZ()[0] != nullptr);
        REQUIRE(synthetic_data_configurations.GetDescriptorProduct()[0] != nullptr);
        REQUIRE(synthetic_data_configurations.GetDescriptorZcpy() != nullptr);
        REQUIRE(synthetic_data_configurations.GetDescriptorDeterminant() != nullptr);

        linearAlgebraSolver->DestoryDescriptors();
        // Finalise Hardware.
        exageostat::api::ExaGeoStat<float>::ExaGeoStatFinalizeHardware(&synthetic_data_configurations);

    }

    SECTION("Double") {
        // Initialise Hardware.
        exageostat::api::ExaGeoStat<double>::ExaGeoStatInitializeHardware(&synthetic_data_configurations);

        auto linearAlgebraSolver = LinearAlgebraFactory<double>::CreateLinearAlgebraSolver(EXACT_DENSE);

        synthetic_data_configurations.SetProblemSize(16);
        synthetic_data_configurations.SetDenseTileSize(8);
        linearAlgebraSolver->SetConfigurations(&synthetic_data_configurations);

        linearAlgebraSolver->InitiateDescriptors();

        REQUIRE(synthetic_data_configurations.GetDescriptorC().size() == 4);
        REQUIRE(synthetic_data_configurations.GetDescriptorZ().size() == 3);
        REQUIRE(synthetic_data_configurations.GetDescriptorProduct().size() == 3);

        for (auto &descriptorC: synthetic_data_configurations.GetDescriptorC()) {
            REQUIRE(descriptorC != nullptr);
        }
        for (auto &i: synthetic_data_configurations.GetDescriptorZ()) {
            REQUIRE(i != nullptr);
        }
        for (auto &i: synthetic_data_configurations.GetDescriptorProduct()) {
            REQUIRE(i != nullptr);
        }
        REQUIRE(synthetic_data_configurations.GetDescriptorZcpy() != nullptr);
        REQUIRE(synthetic_data_configurations.GetDescriptorDeterminant() != nullptr);

        linearAlgebraSolver->DestoryDescriptors();
        // Finalise Hardware.
        exageostat::api::ExaGeoStat<double>::ExaGeoStatFinalizeHardware(&synthetic_data_configurations);
    }
}

//Test that the function initializes the CHAM_descriptorC descriptor correctly.
void TEST_CHAMELEON_DESCRIPTORS_VALUES() {
    SyntheticDataConfigurations synthetic_data_configurations;

    SECTION("SINGLE") {

        // Initialise Hardware.
        exageostat::api::ExaGeoStat<float>::ExaGeoStatInitializeHardware(&synthetic_data_configurations);

        auto linearAlgebraSolver = LinearAlgebraFactory<float>::CreateLinearAlgebraSolver(EXACT_DENSE);

        synthetic_data_configurations.SetProblemSize(4);
        synthetic_data_configurations.SetDenseTileSize(1);
        linearAlgebraSolver->SetConfigurations(&synthetic_data_configurations);

        linearAlgebraSolver->InitiateDescriptors();
        auto *CHAM_descriptorC = (CHAM_desc_t *) synthetic_data_configurations.GetDescriptorC()[0];
        auto *CHAM_descriptorZ = (CHAM_desc_t *) synthetic_data_configurations.GetDescriptorZ()[0];
        auto *CHAM_descriptorZcpy = (CHAM_desc_t *) synthetic_data_configurations.GetDescriptorZcpy();
        auto *CHAM_descriptorDeterminant = (CHAM_desc_t *) synthetic_data_configurations.GetDescriptorDeterminant();
        auto *CHAM_descriptorProduct = (CHAM_desc_t *) synthetic_data_configurations.GetDescriptorProduct()[0];

        int N = synthetic_data_configurations.GetProblemSize() * synthetic_data_configurations.GetP();
        int dts = synthetic_data_configurations.GetDenseTileSize();
        int pGrid = synthetic_data_configurations.GetPGrid();
        int qGrid = synthetic_data_configurations.GetQGrid();

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


        auto *mat = (float *) CHAM_descriptorZ->mat;
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
        linearAlgebraSolver->DestoryDescriptors();
        // Finalise Hardware.
        exageostat::api::ExaGeoStat<float>::ExaGeoStatFinalizeHardware(&synthetic_data_configurations);
    }

    SECTION("DOUBLE") {

        // Initialise Hardware.
        exageostat::api::ExaGeoStat<double>::ExaGeoStatInitializeHardware(&synthetic_data_configurations);

        auto linearAlgebraSolver = LinearAlgebraFactory<double>::CreateLinearAlgebraSolver(EXACT_DENSE);

        synthetic_data_configurations.SetProblemSize(64);
        synthetic_data_configurations.SetDenseTileSize(8);
        linearAlgebraSolver->SetConfigurations(&synthetic_data_configurations);

        linearAlgebraSolver->InitiateDescriptors();
        auto *CHAM_descriptorC = (CHAM_desc_t *) synthetic_data_configurations.GetDescriptorC()[0];
        auto *CHAM_descsubC11 = (CHAM_desc_t *) synthetic_data_configurations.GetDescriptorC()[1];
        auto *CHAM_descsubC12 = (CHAM_desc_t *) synthetic_data_configurations.GetDescriptorC()[2];
        auto *CHAM_descsubC22 = (CHAM_desc_t *) synthetic_data_configurations.GetDescriptorC()[3];
        auto *CHAM_descriptorZ = (CHAM_desc_t *) synthetic_data_configurations.GetDescriptorZ()[0];
        auto *CHAM_descriptorZcpy = (CHAM_desc_t *) synthetic_data_configurations.GetDescriptorZcpy();
        auto *CHAM_descriptorDeterminant = (CHAM_desc_t *) synthetic_data_configurations.GetDescriptorDeterminant();
        vector<void *> &pDescriptorProduct = synthetic_data_configurations.GetDescriptorProduct();
        vector<void *> &pDescriptorZ = synthetic_data_configurations.GetDescriptorZ();


        int N = synthetic_data_configurations.GetProblemSize() * synthetic_data_configurations.GetP();
        int dts = synthetic_data_configurations.GetDenseTileSize();
        int pGrid = synthetic_data_configurations.GetPGrid();
        int qGrid = synthetic_data_configurations.GetQGrid();

        REQUIRE(CHAM_descriptorC->m == N);
        REQUIRE(CHAM_descriptorZ->m == N);
        REQUIRE(CHAM_descriptorZcpy->m == N);
        REQUIRE(CHAM_descriptorDeterminant->m == 1);
        REQUIRE(CHAM_descsubC11->m == N / 2);
        REQUIRE(CHAM_descsubC12->m == N / 2);
        REQUIRE(CHAM_descsubC22->m == N / 2);
        for (auto &idx: pDescriptorProduct) {
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
        for (auto &idx: pDescriptorProduct) {
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
        for (auto &idx: pDescriptorProduct) {
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
        for (auto &idx: pDescriptorProduct) {
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
        for (auto &idx: pDescriptorProduct) {
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
        for (auto &idx: pDescriptorProduct) {
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
        for (auto &idx: pDescriptorProduct) {
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
        REQUIRE(CHAM_descsubC11->mt == ((N / 2) / dts));
        REQUIRE(CHAM_descsubC12->mt == ((N / 2) / dts));
        REQUIRE(CHAM_descsubC22->mt == ((N / 2) / dts));
        for (auto &idx: pDescriptorProduct) {
            auto **CHAM_descriptorProduct = (CHAM_desc_t **) &idx;
            REQUIRE((*CHAM_descriptorProduct)->mt == 1);
        }
        for (int idx = 1; idx < pDescriptorZ.size(); idx++) {
            auto **CHAM_descriptorZ_ = (CHAM_desc_t **) &pDescriptorZ[idx];
            REQUIRE((*CHAM_descriptorZ_)->mt == ((N / 2) / dts));
        }

        REQUIRE(CHAM_descriptorC->nt == ceil((N * 1.0) / (dts * 1.0)));
        REQUIRE(CHAM_descriptorZ->nt == 1);
        REQUIRE(CHAM_descriptorZcpy->nt == 1);
        REQUIRE(CHAM_descriptorDeterminant->nt == 1);
        REQUIRE(CHAM_descsubC11->nt == ((N / 2) / dts));
        REQUIRE(CHAM_descsubC12->nt == ((N / 2) / dts));
        REQUIRE(CHAM_descsubC22->nt == ((N / 2) / dts));
        for (auto &idx: pDescriptorProduct) {
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
        for (auto &idx: pDescriptorProduct) {
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
        for (auto &idx: pDescriptorProduct) {
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
        for (auto &idx: pDescriptorProduct) {
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
        for (auto &idx: pDescriptorProduct) {
            auto **CHAM_descriptorProduct = (CHAM_desc_t **) &idx;
            REQUIRE((*CHAM_descriptorProduct)->q == qGrid);
        }
        for (int idx = 1; idx < pDescriptorZ.size(); idx++) {
            auto **CHAM_descriptorZ_ = (CHAM_desc_t **) &pDescriptorZ[idx];
            REQUIRE((*CHAM_descriptorZ_)->q == qGrid);
        }

        auto *mat = (double *) CHAM_descriptorZ->mat;
        auto *matZcpy = (double *) CHAM_descriptorZcpy->mat;
        for (auto i = 0;
             i < (CHAM_descriptorZ->mt - 1) * (CHAM_descriptorZ->nt - 1) * (CHAM_descriptorZ->bsiz - 1); i++) {
            REQUIRE(mat[i] == 0.0f);
            REQUIRE(matZcpy[i] == 0.0f);
        }
        for (auto &idx: pDescriptorZ) {
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
            for (auto &idx: pDescriptorProduct) {
                auto **CHAM_descriptorProduct = (CHAM_desc_t **) &idx;
                auto *matProduct = (double *) (*CHAM_descriptorProduct)->mat;
                REQUIRE(matProduct[i] == 0.0f);
            }
        }
        linearAlgebraSolver->DestoryDescriptors();
        // Finalise Hardware.
        exageostat::api::ExaGeoStat<double>::ExaGeoStatFinalizeHardware(&synthetic_data_configurations);
    }
}

void TEST_CHAMELEON_GENERATE_OBSERVATIONS() {
    SECTION("Data generation - Observations") {
        // Create a new synthetic_data_configurations object with the provided command line arguments
        SyntheticDataConfigurations synthetic_data_configurations;
        synthetic_data_configurations.SetProblemSize(9);
        synthetic_data_configurations.SetKernel("UnivariateMaternStationary");
#ifdef EXAGEOSTAT_USE_CHAMELEON
        synthetic_data_configurations.SetDenseTileSize(5);
        synthetic_data_configurations.SetComputation(EXACT_DENSE);
#endif
#ifdef EXAGEOSTAT_USE_HiCMA
        synthetic_data_configurations.SetLowTileSize(5);
        synthetic_data_configurations.SetComputation(TILE_LOW_RANK);
#endif
        synthetic_data_configurations.SetDimension(Dimension2D);
        synthetic_data_configurations.SetIsSynthetic(true);
        synthetic_data_configurations.SetPrecision(DOUBLE);

        // Create a unique pointer to a DataGenerator object
        std::unique_ptr<exageostat::generators::DataGenerator<double>> synthetic_generator;

        vector<double> lb{0.1, 0.1, 0.1};
        synthetic_data_configurations.SetLowerBounds(lb);

        vector<double> ub{5, 5, 5};
        synthetic_data_configurations.SetUpperBounds(ub);

        vector<double> initial_theta{1, 0.1, 0.5};
        synthetic_data_configurations.SetInitialTheta(initial_theta);

        // Initialise ExaGeoStat Hardware.
        exageostat::api::ExaGeoStat<double>::ExaGeoStatInitializeHardware(&synthetic_data_configurations);

        // Create the DataGenerator object
        synthetic_generator = synthetic_generator->CreateGenerator(&synthetic_data_configurations);

        // Initialize the seed manually with zero, to get the first generated seeded numbers.
        srand(0);
        // Generated locations data
        synthetic_generator->GenerateLocations();
        synthetic_generator->GenerateDescriptors();

        auto descriptorC = synthetic_data_configurations.GetDescriptorC()[0];
        exageostat::dataunits::Locations *l1 = synthetic_generator->GetLocations();

        synthetic_generator->GetLinearAlgberaSolver()->GenerateObservationsVector(descriptorC, l1, l1, nullptr,
                                                                                  synthetic_data_configurations.GetInitialTheta(),
                                                                                  0, synthetic_generator->GetKernel());

//        auto *A = synthetic_generator->GetLinearAlgberaSolver()->GetMatrix();

        synthetic_generator->DestoryDescriptors();

        // Finalize ExaGeoStat Hardware.
        exageostat::api::ExaGeoStat<double>::ExaGeoStatFinalizeHardware(&synthetic_data_configurations);

    }
}

TEST_CASE("Chameleon Implementation Dense") {
    INIT_FINALIZE_HARDWARE();
    TEST_DESCRIPTORS_INITIALIZATION();
    TEST_CHAMELEON_DESCRIPTORS_VALUES();
    TEST_CHAMELEON_GENERATE_OBSERVATIONS();
}
