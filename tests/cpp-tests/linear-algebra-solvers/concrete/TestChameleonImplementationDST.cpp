
// Copyright (c) 2017-2023 King Abdullah University of Science and Technology,
// Copyright (C) 2023 by Brightskies inc,
// All rights reserved.
// ExaGeoStat is a software package, provided by King Abdullah University of Science and Technology (KAUST).

/**
 * @file TestChameleonImplementationDST.cpp
 *
 * @version 1.0.0
 * @author Sameh Abdulah
 * @date 2023-04-09
**/

#include <libraries/catch/catch.hpp>
#include <cmath>
#include <vector>
#include <linear-algebra-solvers/LinearAlgebraFactory.hpp>
#include <configurations/data-generation/concrete/SyntheticDataConfigurations.hpp>
#include <control/context.h>

using namespace exageostat::linearAlgebra::diagonalSuperTile;
using namespace exageostat::linearAlgebra;
using namespace exageostat::common;
using namespace exageostat::configurations::data_configurations;
using namespace std;

void INIT_HARDWARE_DST() {
    ChameleonImplementationDST<double> chameleonImpl;
    chameleonImpl.ExaGeoStatInitContext(4, 0);
    CHAM_context_t *chameleonContext = chameleon_context_self();
    REQUIRE(chameleonContext != nullptr);
}
void FINALIZE_HARDWARE_DST() {
    ChameleonImplementationDST<double> chameleonImpl;
    chameleonImpl.ExaGeoStatInitContext(4, 0);
    chameleonImpl.ExaGeoStatFinalizeContext();
    REQUIRE(chameleon_context_self() == nullptr);
}

// Test that the function initializes all the required descriptors without errors.
void TEST_INITIALIZETION_DST() {
    auto *syntheticDataConfigurations = new SyntheticDataConfigurations();

    SECTION("Single") {
        auto linearAlgebraSolver = LinearAlgebraFactory<float>::CreateLinearAlgebraSolver(DIAGONAL_APPROX);

        syntheticDataConfigurations->SetProblemSize(1000);
        syntheticDataConfigurations->SetDenseTileSize(64);
        linearAlgebraSolver->SetConfigurations(syntheticDataConfigurations);

        REQUIRE(syntheticDataConfigurations->GetDescriptorZcpy() == nullptr);
        REQUIRE(syntheticDataConfigurations->GetDescriptorDeterminant() == nullptr);

        linearAlgebraSolver->InitiateDescriptors();

        REQUIRE(syntheticDataConfigurations->GetDescriptorC().size() == 1);
        REQUIRE(syntheticDataConfigurations->GetDescriptorZ().size() == 1);
        REQUIRE(syntheticDataConfigurations->GetDescriptorProduct().size() == 1);

        REQUIRE(syntheticDataConfigurations->GetDescriptorC()[0] != nullptr);
        REQUIRE(syntheticDataConfigurations->GetDescriptorZ()[0] != nullptr);
        REQUIRE(syntheticDataConfigurations->GetDescriptorProduct()[0] != nullptr);
        REQUIRE(syntheticDataConfigurations->GetDescriptorZcpy() != nullptr);
        REQUIRE(syntheticDataConfigurations->GetDescriptorDeterminant() != nullptr);
    }

    SECTION("Double") {
        auto linearAlgebraSolver = LinearAlgebraFactory<double>::CreateLinearAlgebraSolver(DIAGONAL_APPROX);

        syntheticDataConfigurations->SetProblemSize(1024);
        syntheticDataConfigurations->SetDenseTileSize(64);
        linearAlgebraSolver->SetConfigurations(syntheticDataConfigurations);

        REQUIRE(syntheticDataConfigurations->GetDescriptorZcpy() == nullptr);
        REQUIRE(syntheticDataConfigurations->GetDescriptorDeterminant() == nullptr);

        linearAlgebraSolver->InitiateDescriptors();

        REQUIRE(syntheticDataConfigurations->GetDescriptorC().size() == 1);
        REQUIRE(syntheticDataConfigurations->GetDescriptorZ().size() == 1);
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
}

//Test that the function initializes the CHAM_descriptorC descriptor correctly.
void TEST_CHAMELEON_DESCRIPTORS_VALUES_DST() {
    auto *syntheticDataConfigurations = new SyntheticDataConfigurations();

    SECTION("SINGLE") {
        auto linearAlgebraSolver = LinearAlgebraFactory<float>::CreateLinearAlgebraSolver(DIAGONAL_APPROX);

        syntheticDataConfigurations->SetProblemSize(6400);
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
    }

    SECTION("DOUBLE") {
        auto linearAlgebraSolver = LinearAlgebraFactory<double>::CreateLinearAlgebraSolver(DIAGONAL_APPROX);

        syntheticDataConfigurations->SetProblemSize(6400);
        syntheticDataConfigurations->SetDenseTileSize(512);
        linearAlgebraSolver->SetConfigurations(syntheticDataConfigurations);

        linearAlgebraSolver->InitiateDescriptors();
        auto *CHAM_descriptorC = (CHAM_desc_t *) syntheticDataConfigurations->GetDescriptorC()[0];
        auto *CHAM_descriptorZ = (CHAM_desc_t *) syntheticDataConfigurations->GetDescriptorZ()[0];
        auto *CHAM_descriptorZcpy = (CHAM_desc_t *) syntheticDataConfigurations->GetDescriptorZcpy();
        auto *CHAM_descriptorDeterminant = (CHAM_desc_t *) syntheticDataConfigurations->GetDescriptorDeterminant();
        vector<void *> &pDescriptorProduct = syntheticDataConfigurations->GetDescriptorProduct();

        int N = syntheticDataConfigurations->GetProblemSize() * syntheticDataConfigurations->GetP();
        int dts = syntheticDataConfigurations->GetDenseTileSize();
        int pGrid = syntheticDataConfigurations->GetPGrid();
        int qGrid = syntheticDataConfigurations->GetQGrid();

        REQUIRE(CHAM_descriptorC->m == N);
        REQUIRE(CHAM_descriptorZ->m == N);
        REQUIRE(CHAM_descriptorZcpy->m == N);
        REQUIRE(CHAM_descriptorDeterminant->m == 1);
        for (auto & idx : pDescriptorProduct) {
            auto **CHAM_descriptorProduct = (CHAM_desc_t **) &idx;
            REQUIRE((*CHAM_descriptorProduct)->m == 1);
        }

        REQUIRE(CHAM_descriptorC->n == N);
        REQUIRE(CHAM_descriptorZ->n == 1);
        REQUIRE(CHAM_descriptorZcpy->n == 1);
        REQUIRE(CHAM_descriptorDeterminant->n == 1);
        for (auto & idx : pDescriptorProduct) {
            auto **CHAM_descriptorProduct = (CHAM_desc_t **) &idx;
            REQUIRE((*CHAM_descriptorProduct)->n == 1);
        }

        REQUIRE(CHAM_descriptorC->mb == dts);
        REQUIRE(CHAM_descriptorZ->mb == dts);
        REQUIRE(CHAM_descriptorZcpy->mb == dts);
        REQUIRE(CHAM_descriptorDeterminant->mb == dts);
        for (auto & idx : pDescriptorProduct) {
            auto **CHAM_descriptorProduct = (CHAM_desc_t **) &idx;
            REQUIRE((*CHAM_descriptorProduct)->mb == dts);
        }

        REQUIRE(CHAM_descriptorC->nb == dts);
        REQUIRE(CHAM_descriptorZ->nb == dts);
        REQUIRE(CHAM_descriptorZcpy->nb == dts);
        REQUIRE(CHAM_descriptorDeterminant->nb == dts);
        for (auto & idx : pDescriptorProduct) {
            auto **CHAM_descriptorProduct = (CHAM_desc_t **) &idx;
            REQUIRE((*CHAM_descriptorProduct)->nb == dts);
        }

        REQUIRE(CHAM_descriptorC->bsiz == dts * dts);
        REQUIRE(CHAM_descriptorZ->bsiz == dts * dts);
        REQUIRE(CHAM_descriptorZcpy->bsiz == dts * dts);
        REQUIRE(CHAM_descriptorDeterminant->bsiz == dts * dts);
        for (auto & idx : pDescriptorProduct) {
            auto **CHAM_descriptorProduct = (CHAM_desc_t **) &idx;
            REQUIRE((*CHAM_descriptorProduct)->bsiz == dts * dts);
        }

        REQUIRE(CHAM_descriptorC->i == 0);
        REQUIRE(CHAM_descriptorZ->i == 0);
        REQUIRE(CHAM_descriptorZcpy->i == 0);
        REQUIRE(CHAM_descriptorDeterminant->i == 0);
        for (auto & idx : pDescriptorProduct) {
            auto **CHAM_descriptorProduct = (CHAM_desc_t **) &idx;
            REQUIRE((*CHAM_descriptorProduct)->i == 0);
        }

        REQUIRE(CHAM_descriptorC->j == 0);
        REQUIRE(CHAM_descriptorZ->j == 0);
        REQUIRE(CHAM_descriptorZcpy->j == 0);
        REQUIRE(CHAM_descriptorDeterminant->j == 0);
        for (auto & idx : pDescriptorProduct) {
            auto **CHAM_descriptorProduct = (CHAM_desc_t **) &idx;
            REQUIRE((*CHAM_descriptorProduct)->j == 0);
        }

        REQUIRE(CHAM_descriptorC->mt == ceil((N * 1.0) / (dts * 1.0)));
        REQUIRE(CHAM_descriptorZ->mt == ceil((N * 1.0) / (dts * 1.0)));
        REQUIRE(CHAM_descriptorZcpy->mt == ceil((N * 1.0) / (dts * 1.0)));
        REQUIRE(CHAM_descriptorDeterminant->mt == 1);
        for (auto & idx : pDescriptorProduct) {
            auto **CHAM_descriptorProduct = (CHAM_desc_t **) &idx;
            REQUIRE((*CHAM_descriptorProduct)->mt == 1);
        }

        REQUIRE(CHAM_descriptorC->nt == ceil((N * 1.0) / (dts * 1.0)));
        REQUIRE(CHAM_descriptorZ->nt == 1);
        REQUIRE(CHAM_descriptorZcpy->nt == 1);
        REQUIRE(CHAM_descriptorDeterminant->nt == 1);
        for (auto & idx : pDescriptorProduct) {
            auto **CHAM_descriptorProduct = (CHAM_desc_t **) &idx;
            REQUIRE((*CHAM_descriptorProduct)->nt == 1);
        }

        REQUIRE(CHAM_descriptorC->lm == N);
        REQUIRE(CHAM_descriptorZ->lm == N);
        REQUIRE(CHAM_descriptorZcpy->lm == N);
        REQUIRE(CHAM_descriptorDeterminant->lm == 1);
        for (auto & idx : pDescriptorProduct) {
            auto **CHAM_descriptorProduct = (CHAM_desc_t **) &idx;
            REQUIRE((*CHAM_descriptorProduct)->lm == 1);
        }

        REQUIRE(CHAM_descriptorC->ln == N);
        REQUIRE(CHAM_descriptorZ->ln == 1);
        REQUIRE(CHAM_descriptorZcpy->ln == 1);
        REQUIRE(CHAM_descriptorDeterminant->ln == 1);
        for (auto & idx : pDescriptorProduct) {
            auto **CHAM_descriptorProduct = (CHAM_desc_t **) &idx;
            REQUIRE((*CHAM_descriptorProduct)->ln == 1);
        }

        REQUIRE(CHAM_descriptorC->p == pGrid);
        REQUIRE(CHAM_descriptorZ->p == pGrid);
        REQUIRE(CHAM_descriptorZcpy->p == pGrid);
        REQUIRE(CHAM_descriptorDeterminant->p == pGrid);
        for (auto & idx : pDescriptorProduct) {
            auto **CHAM_descriptorProduct = (CHAM_desc_t **) &idx;
            REQUIRE((*CHAM_descriptorProduct)->p == pGrid);
        }

        REQUIRE(CHAM_descriptorC->q == qGrid);
        REQUIRE(CHAM_descriptorZ->q == qGrid);
        REQUIRE(CHAM_descriptorZcpy->q == qGrid);
        REQUIRE(CHAM_descriptorDeterminant->q == qGrid);
        for (auto & idx : pDescriptorProduct) {
            auto **CHAM_descriptorProduct = (CHAM_desc_t **) &idx;
            REQUIRE((*CHAM_descriptorProduct)->q == qGrid);
        }
        auto *mat = (double *) CHAM_descriptorC->mat;
        for (auto i = 0;
             i < (CHAM_descriptorC->mt - 1) * (CHAM_descriptorC->nt - 1) * (CHAM_descriptorC->bsiz - 1); i++) {
            REQUIRE(mat[i] == 0.0f);
        }
        mat = (double *) CHAM_descriptorZ->mat;
        auto *matZcpy = (double *) CHAM_descriptorZcpy->mat;
        for (auto i = 0;
             i < (CHAM_descriptorZ->mt - 1) * (CHAM_descriptorZ->nt - 1) * (CHAM_descriptorZ->bsiz - 1); i++) {
            REQUIRE(mat[i] == 0.0f);
            REQUIRE(matZcpy[i] == 0.0f);
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
    }
}

TEST_CASE("Chameleon Implementation DST") {
    INIT_HARDWARE_DST();
    FINALIZE_HARDWARE_DST();
    TEST_INITIALIZETION_DST();
    TEST_CHAMELEON_DESCRIPTORS_VALUES_DST();
}