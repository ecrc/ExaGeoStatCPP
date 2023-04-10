
// Copyright (c) 2017-2023 King Abdullah University of Science and Technology,
// Copyright (C) 2023 by Brightskies inc,
// All rights reserved.
// ExaGeoStat is a software package, provided by King Abdullah University of Science and Technology (KAUST).

/**
 * @file TestHiCMAImplementationTLR.cpp
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
#include <control/hicma_context.h>

using namespace exageostat::linearAlgebra::tileLowRank;
using namespace exageostat::linearAlgebra;
using namespace exageostat::common;
using namespace exageostat::configurations::data_configurations;
using namespace std;

void INIT_HARDWARE_TLR() {
    HicmaImplementation<double> hicmaImpl;
    hicmaImpl.ExaGeoStatInitContext(4, 0);
    HICMA_context_t *hicmaContext = hicma_context_self();
    REQUIRE(hicmaContext != nullptr);
}

void FINALIZE_HARDWARE_TLR() {
    HicmaImplementation<double> hicmaImpl;
    hicmaImpl.ExaGeoStatInitContext(4, 0);
    hicmaImpl.ExaGeoStatFinalizeContext();
    REQUIRE(hicma_context_self() == nullptr);
}

// Test that the function initializes all the required descriptors without errors.
// ONLY DOUBLE IS AVAILABLE FOR NOW.
void TEST_INITIALIZETION_TLR() {
    auto *syntheticDataConfigurations = new SyntheticDataConfigurations();

    SECTION("Double without NZmiss") {
        auto linearAlgebraSolver = LinearAlgebraFactory<double>::CreateLinearAlgebraSolver(TILE_LOW_RANK);

        syntheticDataConfigurations->SetProblemSize(6400);
        syntheticDataConfigurations->SetLowTileSize(512);
        linearAlgebraSolver->SetConfigurations(syntheticDataConfigurations);
        int nZmiss = syntheticDataConfigurations->GetUnknownObservationsNb();

        REQUIRE(syntheticDataConfigurations->GetDescriptorZcpy() == nullptr);
        REQUIRE(syntheticDataConfigurations->GetDescriptorDeterminant() == nullptr);
        REQUIRE(syntheticDataConfigurations->GetDescriptorZObservations() == nullptr);
        REQUIRE(syntheticDataConfigurations->GetDescriptorZActual() == nullptr);
        REQUIRE(syntheticDataConfigurations->GetDescriptorMSE() == nullptr);

        linearAlgebraSolver->InitiateDescriptors();

        REQUIRE(syntheticDataConfigurations->GetDescriptorC().size() == 1);
        REQUIRE(syntheticDataConfigurations->GetDescriptorZ().size() == 1);

        for (auto &descriptorC: syntheticDataConfigurations->GetDescriptorC()) {
            REQUIRE(descriptorC != nullptr);
        }
        for (auto &i: syntheticDataConfigurations->GetDescriptorZ()) {
            REQUIRE(i != nullptr);
        }
        if(nZmiss != 0){
            for (auto &i: syntheticDataConfigurations->GetDescriptorCD()) {
                REQUIRE(i != nullptr);
            }
            for (auto &i: syntheticDataConfigurations->GetDescriptorCUV()) {
                REQUIRE(i != nullptr);
            }
            for (auto &i: syntheticDataConfigurations->GetDescriptorCrk()) {
                REQUIRE(i != nullptr);
            }
            REQUIRE(syntheticDataConfigurations->GetDescriptorZObservations() != nullptr);
            REQUIRE(syntheticDataConfigurations->GetDescriptorZActual() != nullptr);
            REQUIRE(syntheticDataConfigurations->GetDescriptorMSE() != nullptr);
        }
        // Since HiCMA doesn't need product descriptor
        for (auto &i: syntheticDataConfigurations->GetDescriptorProduct()) {
            REQUIRE(i == nullptr);
        }

        REQUIRE(syntheticDataConfigurations->GetDescriptorZcpy() != nullptr);
        REQUIRE(syntheticDataConfigurations->GetDescriptorDeterminant() != nullptr);
    }
    SECTION("Double WITH NZmiss") {
        auto linearAlgebraSolver = LinearAlgebraFactory<double>::CreateLinearAlgebraSolver(TILE_LOW_RANK);

        syntheticDataConfigurations->SetProblemSize(6400);
        syntheticDataConfigurations->SetLowTileSize(512);
        syntheticDataConfigurations->SetUnknownObservationsNb(10);
        linearAlgebraSolver->SetConfigurations(syntheticDataConfigurations);
        int nZmiss = syntheticDataConfigurations->GetUnknownObservationsNb();

        REQUIRE(syntheticDataConfigurations->GetDescriptorZcpy() == nullptr);
        REQUIRE(syntheticDataConfigurations->GetDescriptorDeterminant() == nullptr);
        REQUIRE(syntheticDataConfigurations->GetDescriptorZObservations() == nullptr);
        REQUIRE(syntheticDataConfigurations->GetDescriptorZActual() == nullptr);
        REQUIRE(syntheticDataConfigurations->GetDescriptorMSE() == nullptr);

        linearAlgebraSolver->InitiateDescriptors();

        REQUIRE(syntheticDataConfigurations->GetDescriptorC().size() == 1);
        REQUIRE(syntheticDataConfigurations->GetDescriptorZ().size() == 1);

        for (auto &descriptorC: syntheticDataConfigurations->GetDescriptorC()) {
            REQUIRE(descriptorC != nullptr);
        }
        for (auto &i: syntheticDataConfigurations->GetDescriptorZ()) {
            REQUIRE(i != nullptr);
        }
        if(nZmiss != 0){
            for (auto &i: syntheticDataConfigurations->GetDescriptorCD()) {
                REQUIRE(i != nullptr);
            }
            for (auto &i: syntheticDataConfigurations->GetDescriptorCUV()) {
                REQUIRE(i != nullptr);
            }
            for (auto &i: syntheticDataConfigurations->GetDescriptorCrk()) {
                REQUIRE(i != nullptr);
            }
            REQUIRE(syntheticDataConfigurations->GetDescriptorZObservations() != nullptr);
            REQUIRE(syntheticDataConfigurations->GetDescriptorZActual() != nullptr);
            REQUIRE(syntheticDataConfigurations->GetDescriptorMSE() != nullptr);
        }
        // Since HiCMA doesn't need product descriptor
        for (auto &i: syntheticDataConfigurations->GetDescriptorProduct()) {
            REQUIRE(i == nullptr);
        }

        REQUIRE(syntheticDataConfigurations->GetDescriptorZcpy() != nullptr);
        REQUIRE(syntheticDataConfigurations->GetDescriptorDeterminant() != nullptr);
    }
}

//Test that the function initializes the (*HICMA_descriptorC) descriptor correctly.
void TEST_HICMA_DESCRIPTORS_VALUES_TLR() {
    auto *syntheticDataConfigurations = new SyntheticDataConfigurations();

    SECTION("DOUBLE without NZmiss") {
        auto linearAlgebraSolver = LinearAlgebraFactory<double>::CreateLinearAlgebraSolver(TILE_LOW_RANK);

        syntheticDataConfigurations->SetProblemSize(6400);
        syntheticDataConfigurations->SetLowTileSize(512);
        syntheticDataConfigurations->SetApproximationMode(1);
        linearAlgebraSolver->SetConfigurations(syntheticDataConfigurations);

        linearAlgebraSolver->InitiateDescriptors();

        auto **HICMA_descriptorC = (HICMA_desc_t **) &syntheticDataConfigurations->GetDescriptorC()[0];
        auto **HICMA_descriptorZcpy = (HICMA_desc_t **) &syntheticDataConfigurations->GetDescriptorZcpy();
        auto **HICMA_descriptorDeterminant = (HICMA_desc_t **) &syntheticDataConfigurations->GetDescriptorDeterminant();
        auto **HICMA_descriptorCD = (HICMA_desc_t **) &syntheticDataConfigurations->GetDescriptorCD()[0];
        auto **HICMA_descriptorCUV = (HICMA_desc_t **) &syntheticDataConfigurations->GetDescriptorCUV()[0];
        auto **HICMA_descriptorCrk = (HICMA_desc_t **) &syntheticDataConfigurations->GetDescriptorCrk()[0];
        int approximationMode = syntheticDataConfigurations->GetApproximationMode();
        int N = syntheticDataConfigurations->GetProblemSize() * syntheticDataConfigurations->GetP();
        int lts = syntheticDataConfigurations->GetLowTileSize();
        int pGrid = syntheticDataConfigurations->GetPGrid();
        int qGrid = syntheticDataConfigurations->GetQGrid();
        int maxRank = syntheticDataConfigurations->GetMaxRank();
        int nZmiss = syntheticDataConfigurations->GetUnknownObservationsNb();
        double meanSquareError = syntheticDataConfigurations->GetMeanSquareError();
        string actualObservationsFilePath = syntheticDataConfigurations->GetActualObservationsFilePath();
        double determinantValue = syntheticDataConfigurations->GetDeterminantValue();
        int nZobs = syntheticDataConfigurations->GetKnownObservationsValues();

        if(approximationMode == 1){
            // Descriptor C.
            REQUIRE((*HICMA_descriptorC)->m == N);
            REQUIRE((*HICMA_descriptorC)->n == N);
            REQUIRE((*HICMA_descriptorC)->mb == lts);
            REQUIRE((*HICMA_descriptorC)->nb == lts);
            REQUIRE((*HICMA_descriptorC)->bsiz == lts * lts);
            REQUIRE((*HICMA_descriptorC)->i == 0);
            REQUIRE((*HICMA_descriptorC)->j == 0);
            REQUIRE((*HICMA_descriptorC)->mt == ceil((N * 1.0) / (lts * 1.0)));
            REQUIRE((*HICMA_descriptorC)->nt == ceil((N * 1.0) / (lts * 1.0)));
            REQUIRE((*HICMA_descriptorC)->lm == N);
            REQUIRE((*HICMA_descriptorC)->ln == N);
            REQUIRE((*HICMA_descriptorC)->p == pGrid);
            REQUIRE((*HICMA_descriptorC)->q == qGrid);
        }
        // Re-Run again but with approx mode OFF
        syntheticDataConfigurations->SetApproximationMode(0);
        linearAlgebraSolver->SetConfigurations(syntheticDataConfigurations);

        linearAlgebraSolver->InitiateDescriptors();
        approximationMode = syntheticDataConfigurations->GetApproximationMode();
        auto **HICMA_descriptorZ = (HICMA_desc_t **) &syntheticDataConfigurations->GetDescriptorZ()[0];


        if(approximationMode != 1){
            // Descriptor C.
            REQUIRE((*HICMA_descriptorC)->m == N);
            REQUIRE((*HICMA_descriptorC)->n == 1);
            REQUIRE((*HICMA_descriptorC)->mb == lts);
            REQUIRE((*HICMA_descriptorC)->nb == lts);
            REQUIRE((*HICMA_descriptorC)->bsiz == lts * lts);
            REQUIRE((*HICMA_descriptorC)->i == 0);
            REQUIRE((*HICMA_descriptorC)->j == 0);
            REQUIRE((*HICMA_descriptorC)->mt == ceil((N * 1.0) / (lts * 1.0)));
            REQUIRE((*HICMA_descriptorC)->nt == 1);
            REQUIRE((*HICMA_descriptorC)->lm == N);
            REQUIRE((*HICMA_descriptorC)->ln == 1);
            REQUIRE((*HICMA_descriptorC)->p == pGrid);
            REQUIRE((*HICMA_descriptorC)->q == qGrid);
        }

        // Descriptor CD.
        REQUIRE((*HICMA_descriptorCD)->m == N);
        REQUIRE((*HICMA_descriptorCD)->n == lts);
        REQUIRE((*HICMA_descriptorCD)->mb == lts);
        REQUIRE((*HICMA_descriptorCD)->nb == lts);
        REQUIRE((*HICMA_descriptorCD)->bsiz == lts * lts);
        REQUIRE((*HICMA_descriptorCD)->i == 0);
        REQUIRE((*HICMA_descriptorCD)->j == 0);
        REQUIRE((*HICMA_descriptorCD)->mt == ceil((N * 1.0) / (lts * 1.0)));
        REQUIRE((*HICMA_descriptorCD)->nt == 1);
        REQUIRE((*HICMA_descriptorCD)->lm == N);
        REQUIRE((*HICMA_descriptorCD)->ln == lts);
        REQUIRE((*HICMA_descriptorCD)->p == pGrid);
        REQUIRE((*HICMA_descriptorCD)->q == qGrid);

        // Descriptor CUV.
        int MUV = N / lts * lts + lts;
        int expr = MUV / lts;
        int NUV = 2 * expr * maxRank;
        int NBUV = 2 * maxRank;
        REQUIRE((*HICMA_descriptorCUV)->m == MUV);
        REQUIRE((*HICMA_descriptorCUV)->n == NUV);
        REQUIRE((*HICMA_descriptorCUV)->mb == lts);
        REQUIRE((*HICMA_descriptorCUV)->nb == NBUV);
        REQUIRE((*HICMA_descriptorCUV)->bsiz == lts * NBUV);
        REQUIRE((*HICMA_descriptorCUV)->i == 0);
        REQUIRE((*HICMA_descriptorCUV)->j == 0);
        REQUIRE((*HICMA_descriptorCUV)->mt == ceil((N * 1.0) / (lts * 1.0)));
        REQUIRE((*HICMA_descriptorCUV)->nt == ceil((N * 1.0) / (lts * 1.0)));
        REQUIRE((*HICMA_descriptorCUV)->lm == MUV);
        REQUIRE((*HICMA_descriptorCUV)->ln == NUV);
        REQUIRE((*HICMA_descriptorCUV)->p == pGrid);
        REQUIRE((*HICMA_descriptorCUV)->q == qGrid);

        // Descriptor Crk.
        REQUIRE((*HICMA_descriptorCrk)->m == (*HICMA_descriptorCD)->mt);
        REQUIRE((*HICMA_descriptorCrk)->n == (*HICMA_descriptorCD)->mt);
        REQUIRE((*HICMA_descriptorCrk)->mb == 1);
        REQUIRE((*HICMA_descriptorCrk)->nb == 1);
        REQUIRE((*HICMA_descriptorCrk)->bsiz == 1);
        REQUIRE((*HICMA_descriptorCrk)->i == 0);
        REQUIRE((*HICMA_descriptorCrk)->j == 0);
        REQUIRE((*HICMA_descriptorCrk)->mt == ceil((N * 1.0) / (lts * 1.0)));
        REQUIRE((*HICMA_descriptorCrk)->nt == ceil((N * 1.0) / (lts * 1.0)));
        REQUIRE((*HICMA_descriptorCrk)->lm == (*HICMA_descriptorCD)->mt);
        REQUIRE((*HICMA_descriptorCrk)->ln == (*HICMA_descriptorCD)->mt);
        REQUIRE((*HICMA_descriptorCrk)->p == pGrid);
        REQUIRE((*HICMA_descriptorCrk)->q == qGrid);

        // Descriptor Z.
        REQUIRE((*HICMA_descriptorZ)->m == N);
        REQUIRE((*HICMA_descriptorZ)->n == 1);
        REQUIRE((*HICMA_descriptorZ)->mb == lts);
        REQUIRE((*HICMA_descriptorZ)->nb == lts);
        REQUIRE((*HICMA_descriptorZ)->bsiz == lts * lts);
        REQUIRE((*HICMA_descriptorZ)->i == 0);
        REQUIRE((*HICMA_descriptorZ)->j == 0);
        REQUIRE((*HICMA_descriptorZ)->mt == ceil((N * 1.0) / (lts * 1.0)));
        REQUIRE((*HICMA_descriptorZ)->nt == 1);
        REQUIRE((*HICMA_descriptorZ)->lm == N);
        REQUIRE((*HICMA_descriptorZ)->ln == 1);
        REQUIRE((*HICMA_descriptorZ)->p == pGrid);
        REQUIRE((*HICMA_descriptorZ)->q == qGrid);

        // Descriptor Zcpy.
        REQUIRE((*HICMA_descriptorZcpy)->m == N);
        REQUIRE((*HICMA_descriptorZcpy)->n == 1);
        REQUIRE((*HICMA_descriptorZcpy)->mb == lts);
        REQUIRE((*HICMA_descriptorZcpy)->nb == lts);
        REQUIRE((*HICMA_descriptorZcpy)->bsiz == lts * lts);
        REQUIRE((*HICMA_descriptorZcpy)->i == 0);
        REQUIRE((*HICMA_descriptorZcpy)->j == 0);
        REQUIRE((*HICMA_descriptorZcpy)->mt == ceil((N * 1.0) / (lts * 1.0)));
        REQUIRE((*HICMA_descriptorZcpy)->nt == 1);
        REQUIRE((*HICMA_descriptorZcpy)->lm == N);
        REQUIRE((*HICMA_descriptorZcpy)->ln == 1);
        REQUIRE((*HICMA_descriptorZcpy)->p == pGrid);
        REQUIRE((*HICMA_descriptorZcpy)->q == qGrid);

        // Descriptor Determinant.
        REQUIRE((*HICMA_descriptorDeterminant)->m == 1);
        REQUIRE((*HICMA_descriptorDeterminant)->n == 1);
        REQUIRE((*HICMA_descriptorDeterminant)->mb == lts);
        REQUIRE((*HICMA_descriptorDeterminant)->nb == lts);
        REQUIRE((*HICMA_descriptorDeterminant)->bsiz == lts * lts);
        REQUIRE((*HICMA_descriptorDeterminant)->i == 0);
        REQUIRE((*HICMA_descriptorDeterminant)->j == 0);
        REQUIRE((*HICMA_descriptorDeterminant)->mt == 1);
        REQUIRE((*HICMA_descriptorDeterminant)->nt == 1);
        REQUIRE((*HICMA_descriptorDeterminant)->lm == 1);
        REQUIRE((*HICMA_descriptorDeterminant)->ln == 1);
        REQUIRE((*HICMA_descriptorDeterminant)->p == pGrid);
        REQUIRE((*HICMA_descriptorDeterminant)->q == qGrid);
    }

    SECTION("DOUBLE WITH NZmiss") {
        auto linearAlgebraSolver = LinearAlgebraFactory<double>::CreateLinearAlgebraSolver(TILE_LOW_RANK);

        syntheticDataConfigurations->SetProblemSize(6400);
        syntheticDataConfigurations->SetLowTileSize(512);
        syntheticDataConfigurations->SetUnknownObservationsNb(10);
        linearAlgebraSolver->SetConfigurations(syntheticDataConfigurations);

        linearAlgebraSolver->InitiateDescriptors();

        auto **HICMA_descriptorZObservations = (HICMA_desc_t **) &syntheticDataConfigurations->GetDescriptorZObservations();
        auto **HICMA_descriptorZactual = (HICMA_desc_t **) &syntheticDataConfigurations->GetDescriptorZActual();
        auto **HICMA_descriptorMSE = (HICMA_desc_t **) &syntheticDataConfigurations->GetDescriptorMSE();

        auto **HICMA_descriptorC12D = (HICMA_desc_t **) &syntheticDataConfigurations->GetDescriptorCD()[1];
        auto **HICMA_descriptorC22D = (HICMA_desc_t **) &syntheticDataConfigurations->GetDescriptorCD()[2];

        auto **HICMA_descriptorC12UV = (HICMA_desc_t **) &syntheticDataConfigurations->GetDescriptorCUV()[1];
        auto **HICMA_descriptorC22UV = (HICMA_desc_t **) &syntheticDataConfigurations->GetDescriptorCUV()[2];

        auto **HICMA_descriptorC12rk = (HICMA_desc_t **) &syntheticDataConfigurations->GetDescriptorCrk()[1];
        auto **HICMA_descriptorC22rk = (HICMA_desc_t **) &syntheticDataConfigurations->GetDescriptorCrk()[2];

        int N = syntheticDataConfigurations->GetProblemSize() * syntheticDataConfigurations->GetP();
        int lts = syntheticDataConfigurations->GetLowTileSize();
        int pGrid = syntheticDataConfigurations->GetPGrid();
        int qGrid = syntheticDataConfigurations->GetQGrid();
        int maxRank = syntheticDataConfigurations->GetMaxRank();
        int nZmiss = syntheticDataConfigurations->GetUnknownObservationsNb();
        string actualObservationsFilePath = syntheticDataConfigurations->GetActualObservationsFilePath();
        int nZobs = syntheticDataConfigurations->GetKnownObservationsValues();

        if (nZmiss != 0) {
            if (actualObservationsFilePath.empty()) {
                // Descriptor ZObservations.
                REQUIRE((*HICMA_descriptorZObservations)->m == nZobs);
                REQUIRE((*HICMA_descriptorZObservations)->n == 1);
                REQUIRE((*HICMA_descriptorZObservations)->mb == lts);
                REQUIRE((*HICMA_descriptorZObservations)->nb == lts);
                REQUIRE((*HICMA_descriptorZObservations)->bsiz == lts * lts);
                REQUIRE((*HICMA_descriptorZObservations)->i == 0);
                REQUIRE((*HICMA_descriptorZObservations)->j == 0);
                REQUIRE((*HICMA_descriptorZObservations)->mt == ceil((N * 1.0) / (lts * 1.0)));
                REQUIRE((*HICMA_descriptorZObservations)->nt == 1);
                REQUIRE((*HICMA_descriptorZObservations)->lm == nZobs);
                REQUIRE((*HICMA_descriptorZObservations)->ln == 1);
                REQUIRE((*HICMA_descriptorZObservations)->p == pGrid);
                REQUIRE((*HICMA_descriptorZObservations)->q == qGrid);

                // Descriptor Zactual.
                REQUIRE((*HICMA_descriptorZactual)->m == nZmiss);
                REQUIRE((*HICMA_descriptorZactual)->n == 1);
                REQUIRE((*HICMA_descriptorZactual)->mb == lts);
                REQUIRE((*HICMA_descriptorZactual)->nb == lts);
                REQUIRE((*HICMA_descriptorZactual)->bsiz == lts * lts);
                REQUIRE((*HICMA_descriptorZactual)->i == 0);
                REQUIRE((*HICMA_descriptorZactual)->j == 0);
                REQUIRE((*HICMA_descriptorZactual)->mt == 1);
                REQUIRE((*HICMA_descriptorZactual)->nt == 1);
                REQUIRE((*HICMA_descriptorZactual)->lm == nZmiss);
                REQUIRE((*HICMA_descriptorZactual)->ln == 1);
                REQUIRE((*HICMA_descriptorZactual)->p == pGrid);
                REQUIRE((*HICMA_descriptorZactual)->q == qGrid);
            }
            // Descriptor C12D.
            REQUIRE((*HICMA_descriptorC12D)->m == nZmiss);
            REQUIRE((*HICMA_descriptorC12D)->n == lts);
            REQUIRE((*HICMA_descriptorC12D)->mb == lts);
            REQUIRE((*HICMA_descriptorC12D)->nb == lts);
            REQUIRE((*HICMA_descriptorC12D)->bsiz == lts * lts);
            REQUIRE((*HICMA_descriptorC12D)->i == 0);
            REQUIRE((*HICMA_descriptorC12D)->j == 0);
            REQUIRE((*HICMA_descriptorC12D)->mt == 1);
            REQUIRE((*HICMA_descriptorC12D)->nt == 1);
            REQUIRE((*HICMA_descriptorC12D)->lm == nZmiss);
            REQUIRE((*HICMA_descriptorC12D)->ln == lts);
            REQUIRE((*HICMA_descriptorC12D)->p == pGrid);
            REQUIRE((*HICMA_descriptorC12D)->q == qGrid);

            // Descriptor C12UV.
            int NBUV = 2 * maxRank;
            REQUIRE((*HICMA_descriptorC12UV)->m == lts);
            REQUIRE((*HICMA_descriptorC12UV)->n == NBUV);
            REQUIRE((*HICMA_descriptorC12UV)->mb == lts);
            REQUIRE((*HICMA_descriptorC12UV)->nb == NBUV);
            REQUIRE((*HICMA_descriptorC12UV)->bsiz == lts * NBUV);
            REQUIRE((*HICMA_descriptorC12UV)->i == 0);
            REQUIRE((*HICMA_descriptorC12UV)->j == 0);
            REQUIRE((*HICMA_descriptorC12UV)->mt == 1);
            REQUIRE((*HICMA_descriptorC12UV)->nt == 1);
            REQUIRE((*HICMA_descriptorC12UV)->lm == lts);
            REQUIRE((*HICMA_descriptorC12UV)->ln == NBUV);
            REQUIRE((*HICMA_descriptorC12UV)->p == pGrid);
            REQUIRE((*HICMA_descriptorC12UV)->q == qGrid);

            // Descriptor C12rk.
            REQUIRE((*HICMA_descriptorC12rk)->m == 1);
            REQUIRE((*HICMA_descriptorC12rk)->n == 1);
            REQUIRE((*HICMA_descriptorC12rk)->mb == 1);
            REQUIRE((*HICMA_descriptorC12rk)->nb == 1);
            REQUIRE((*HICMA_descriptorC12rk)->bsiz == 1);
            REQUIRE((*HICMA_descriptorC12rk)->i == 0);
            REQUIRE((*HICMA_descriptorC12rk)->j == 0);
            REQUIRE((*HICMA_descriptorC12rk)->mt == 1);
            REQUIRE((*HICMA_descriptorC12rk)->nt == 1);
            REQUIRE((*HICMA_descriptorC12rk)->lm == 1);
            REQUIRE((*HICMA_descriptorC12rk)->ln == 1);
            REQUIRE((*HICMA_descriptorC12rk)->p == pGrid);
            REQUIRE((*HICMA_descriptorC12rk)->q == qGrid);

            // Descriptor C22D.
            REQUIRE((*HICMA_descriptorC22D)->m == nZobs);
            REQUIRE((*HICMA_descriptorC22D)->n == lts);
            REQUIRE((*HICMA_descriptorC22D)->mb == lts);
            REQUIRE((*HICMA_descriptorC22D)->nb == lts);
            REQUIRE((*HICMA_descriptorC22D)->bsiz == lts * lts);
            REQUIRE((*HICMA_descriptorC22D)->i == 0);
            REQUIRE((*HICMA_descriptorC22D)->j == 0);
            REQUIRE((*HICMA_descriptorC22D)->mt == ceil((N * 1.0) / (lts * 1.0)));
            REQUIRE((*HICMA_descriptorC22D)->nt == 1);
            REQUIRE((*HICMA_descriptorC22D)->lm == nZobs);
            REQUIRE((*HICMA_descriptorC22D)->ln == lts);
            REQUIRE((*HICMA_descriptorC22D)->p == pGrid);
            REQUIRE((*HICMA_descriptorC22D)->q == qGrid);

            // Descriptor C22UV.
            REQUIRE((*HICMA_descriptorC22UV)->m == lts);
            REQUIRE((*HICMA_descriptorC22UV)->n == NBUV);
            REQUIRE((*HICMA_descriptorC22UV)->mb == lts);
            REQUIRE((*HICMA_descriptorC22UV)->nb == NBUV);
            REQUIRE((*HICMA_descriptorC22UV)->bsiz == lts * NBUV);
            REQUIRE((*HICMA_descriptorC22UV)->i == 0);
            REQUIRE((*HICMA_descriptorC22UV)->j == 0);
            REQUIRE((*HICMA_descriptorC22UV)->mt == 1);
            REQUIRE((*HICMA_descriptorC22UV)->nt == 1);
            REQUIRE((*HICMA_descriptorC22UV)->lm == lts);
            REQUIRE((*HICMA_descriptorC22UV)->ln == NBUV);
            REQUIRE((*HICMA_descriptorC22UV)->p == pGrid);
            REQUIRE((*HICMA_descriptorC22UV)->q == qGrid);

            // Descriptor C22rk.
            REQUIRE((*HICMA_descriptorC22rk)->m == 1);
            REQUIRE((*HICMA_descriptorC22rk)->n == 1);
            REQUIRE((*HICMA_descriptorC22rk)->mb == 1);
            REQUIRE((*HICMA_descriptorC22rk)->nb == 1);
            REQUIRE((*HICMA_descriptorC22rk)->bsiz == 1);
            REQUIRE((*HICMA_descriptorC22rk)->i == 0);
            REQUIRE((*HICMA_descriptorC22rk)->j == 0);
            REQUIRE((*HICMA_descriptorC22rk)->mt == 1);
            REQUIRE((*HICMA_descriptorC22rk)->nt == 1);
            REQUIRE((*HICMA_descriptorC22rk)->lm == 1);
            REQUIRE((*HICMA_descriptorC22rk)->ln == 1);
            REQUIRE((*HICMA_descriptorC22rk)->p == pGrid);
            REQUIRE((*HICMA_descriptorC22rk)->q == qGrid);

            // Descriptor Determinant.
            REQUIRE((*HICMA_descriptorMSE)->m == 1);
            REQUIRE((*HICMA_descriptorMSE)->n == 1);
            REQUIRE((*HICMA_descriptorMSE)->mb == lts);
            REQUIRE((*HICMA_descriptorMSE)->nb == lts);
            REQUIRE((*HICMA_descriptorMSE)->bsiz == lts * lts);
            REQUIRE((*HICMA_descriptorMSE)->i == 0);
            REQUIRE((*HICMA_descriptorMSE)->j == 0);
            REQUIRE((*HICMA_descriptorMSE)->mt == 1);
            REQUIRE((*HICMA_descriptorMSE)->nt == 1);
            REQUIRE((*HICMA_descriptorMSE)->lm == 1);
            REQUIRE((*HICMA_descriptorMSE)->ln == 1);
            REQUIRE((*HICMA_descriptorMSE)->p == pGrid);
            REQUIRE((*HICMA_descriptorMSE)->q == qGrid);
        }
    }
}

TEST_CASE("HiCMA Implementation TLR") {
    INIT_HARDWARE_TLR();
//    FINALIZE_HARDWARE_TLR();
    TEST_HICMA_DESCRIPTORS_VALUES_TLR();
    TEST_INITIALIZETION_TLR();
}