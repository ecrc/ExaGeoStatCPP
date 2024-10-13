
// Copyright (c) 2017-2024 King Abdullah University of Science and Technology,
// All rights reserved.
// ExaGeoStat is a software package, provided by King Abdullah University of Science and Technology (KAUST).

/**
 * @file TestHiCMAImplementationTLR.cpp
 * @brief Unit tests for the Tile Low Rank computation in the ExaGeoStat software package.
 * @version 1.1.0
 * @author Mahmoud ElKarargy
 * @date 2023-04-09
**/

#include <catch2/catch_all.hpp>
#include <hardware/ExaGeoStatHardware.hpp>
#include <api/ExaGeoStat.hpp>

using namespace std;

using namespace exageostat::common;
using namespace exageostat::dataunits;
using namespace exageostat::configurations;

//Test that the function initializes the HICMA_descriptorC descriptor correctly.
void TEST_HICMA_DESCRIPTORS_VALUES_TLR() {

    SECTION("descriptors values") {

        Configurations synthetic_data_configurations;
        synthetic_data_configurations.SetComputation(exageostat::common::TILE_LOW_RANK);
        synthetic_data_configurations.SetMaxMleIterations(1);
        synthetic_data_configurations.SetTolerance(4);

        vector<double> lb{0.1, 0.1, 0.1};
        synthetic_data_configurations.SetLowerBounds(lb);
        synthetic_data_configurations.SetStartingTheta(lb);
        vector<double> ub{5, 5, 5};
        synthetic_data_configurations.SetUpperBounds(ub);
        vector<double> initial_theta{1, 0.1, 0.5};
        synthetic_data_configurations.SetInitialTheta(initial_theta);
        vector<double> estimated_theta{-1, -1, -1};
        synthetic_data_configurations.SetEstimatedTheta(estimated_theta);
        synthetic_data_configurations.SetKernelName("UnivariateMaternStationary");
        synthetic_data_configurations.SetMaxRank(500);
        synthetic_data_configurations.SetProblemSize(16);
        synthetic_data_configurations.SetLowTileSize(8);
        synthetic_data_configurations.SetDenseTileSize(8);

        // initialize Hardware.
        auto hardware = ExaGeoStatHardware(TILE_LOW_RANK, 1, 0);
        synthetic_data_configurations.SetApproximationMode(1);

        std::unique_ptr<ExaGeoStatData<double>> data;
        exageostat::api::ExaGeoStat<double>::ExaGeoStatLoadData(synthetic_data_configurations, data);
        exageostat::api::ExaGeoStat<double>::ExaGeoStatDataModeling(synthetic_data_configurations, data);

        auto *HICMA_descriptorC = data->GetDescriptorData()->GetDescriptor(HICMA_DESCRIPTOR, DESCRIPTOR_C).hicma_desc;
        int approximationMode = synthetic_data_configurations.GetApproximationMode();
        int N = synthetic_data_configurations.GetProblemSize();
        int lts = synthetic_data_configurations.GetLowTileSize();
        int pGrid = ExaGeoStatHardware::GetPGrid();
        int qGrid = ExaGeoStatHardware::GetQGrid();

        // Descriptor C.
        REQUIRE(HICMA_descriptorC->m == N);
        REQUIRE(HICMA_descriptorC->n == N);
        REQUIRE(HICMA_descriptorC->mb == lts);
        REQUIRE(HICMA_descriptorC->nb == lts);
        REQUIRE(HICMA_descriptorC->bsiz == lts * lts);
        REQUIRE(HICMA_descriptorC->i == 0);
        REQUIRE(HICMA_descriptorC->j == 0);
        REQUIRE(HICMA_descriptorC->mt == ceil((N * 1.0) / (lts * 1.0)));
        REQUIRE(HICMA_descriptorC->nt == ceil((N * 1.0) / (lts * 1.0)));
        REQUIRE(HICMA_descriptorC->lm == N);
        REQUIRE(HICMA_descriptorC->ln == N);
        REQUIRE(HICMA_descriptorC->p == pGrid);
        REQUIRE(HICMA_descriptorC->q == qGrid);

        int maxRank = synthetic_data_configurations.GetMaxRank();
        string actualObservationsFilePath = synthetic_data_configurations.GetActualObservationsFilePath();

        auto *HICMA_descriptorZ = data->GetDescriptorData()->GetDescriptor(HICMA_DESCRIPTOR, DESCRIPTOR_Z).hicma_desc;
        auto *HICMA_descriptorZcpy = data->GetDescriptorData()->GetDescriptor(HICMA_DESCRIPTOR,
                                                                              DESCRIPTOR_Z_COPY).hicma_desc;
        auto *HICMA_descriptorDeterminant = data->GetDescriptorData()->GetDescriptor(HICMA_DESCRIPTOR,
                                                                                     DESCRIPTOR_DETERMINANT).hicma_desc;
        auto *HICMA_descriptorCD = data->GetDescriptorData()->GetDescriptor(HICMA_DESCRIPTOR, DESCRIPTOR_CD).hicma_desc;
        auto *HICMA_descriptorCUV = data->GetDescriptorData()->GetDescriptor(HICMA_DESCRIPTOR,
                                                                             DESCRIPTOR_CUV).hicma_desc;
        auto *HICMA_descriptorCrk = data->GetDescriptorData()->GetDescriptor(HICMA_DESCRIPTOR,
                                                                             DESCRIPTOR_CRK).hicma_desc;

        // Descriptor CD.
        REQUIRE(HICMA_descriptorCD->m == N);
        REQUIRE(HICMA_descriptorCD->n == lts);
        REQUIRE(HICMA_descriptorCD->mb == lts);
        REQUIRE(HICMA_descriptorCD->nb == lts);
        REQUIRE(HICMA_descriptorCD->bsiz == lts * lts);
        REQUIRE(HICMA_descriptorCD->i == 0);
        REQUIRE(HICMA_descriptorCD->j == 0);
        REQUIRE(HICMA_descriptorCD->mt == ceil((N * 1.0) / (lts * 1.0)));
        REQUIRE(HICMA_descriptorCD->nt == 1);
        REQUIRE(HICMA_descriptorCD->lm == N);
        REQUIRE(HICMA_descriptorCD->ln == lts);
        REQUIRE(HICMA_descriptorCD->p == pGrid);
        REQUIRE(HICMA_descriptorCD->q == qGrid);

        // Descriptor CUV.
        int NBUV = 2 * maxRank;
        int MUV = -1;
        int N_over_lts_times_lts = N / lts * lts;
        if (N_over_lts_times_lts < N) {
            MUV = N_over_lts_times_lts + lts;
        } else if (N_over_lts_times_lts == N) {
            MUV = N_over_lts_times_lts;
        }
        int expr = MUV / lts;
        int NUV = 2 * expr * maxRank;
        REQUIRE(HICMA_descriptorCUV->m == N);
        REQUIRE(HICMA_descriptorCUV->n == NUV);
        REQUIRE(HICMA_descriptorCUV->mb == lts);
        REQUIRE(HICMA_descriptorCUV->nb == NBUV);
        REQUIRE(HICMA_descriptorCUV->bsiz == lts * NBUV);
        REQUIRE(HICMA_descriptorCUV->i == 0);
        REQUIRE(HICMA_descriptorCUV->j == 0);
        REQUIRE(HICMA_descriptorCUV->mt == ceil((N * 1.0) / (lts * 1.0)));
        REQUIRE(HICMA_descriptorCUV->nt == ceil((N * 1.0) / (lts * 1.0)));
        REQUIRE(HICMA_descriptorCUV->lm == N);
        REQUIRE(HICMA_descriptorCUV->ln == NUV);
        REQUIRE(HICMA_descriptorCUV->p == pGrid);
        REQUIRE(HICMA_descriptorCUV->q == qGrid);

        // Descriptor Crk.
        REQUIRE(HICMA_descriptorCrk->m == HICMA_descriptorCD->mt);
        REQUIRE(HICMA_descriptorCrk->n == HICMA_descriptorCD->mt);
        REQUIRE(HICMA_descriptorCrk->mb == 1);
        REQUIRE(HICMA_descriptorCrk->nb == 1);
        REQUIRE(HICMA_descriptorCrk->bsiz == 1);
        REQUIRE(HICMA_descriptorCrk->i == 0);
        REQUIRE(HICMA_descriptorCrk->j == 0);
        REQUIRE(HICMA_descriptorCrk->mt == ceil((N * 1.0) / (lts * 1.0)));
        REQUIRE(HICMA_descriptorCrk->nt == ceil((N * 1.0) / (lts * 1.0)));
        REQUIRE(HICMA_descriptorCrk->lm == HICMA_descriptorCD->mt);
        REQUIRE(HICMA_descriptorCrk->ln == HICMA_descriptorCD->mt);
        REQUIRE(HICMA_descriptorCrk->p == pGrid);
        REQUIRE(HICMA_descriptorCrk->q == qGrid);

        // Descriptor Z.
        REQUIRE(HICMA_descriptorZ->m == N);
        REQUIRE(HICMA_descriptorZ->n == 1);
        REQUIRE(HICMA_descriptorZ->mb == lts);
        REQUIRE(HICMA_descriptorZ->nb == lts);
        REQUIRE(HICMA_descriptorZ->bsiz == lts * lts);
        REQUIRE(HICMA_descriptorZ->i == 0);
        REQUIRE(HICMA_descriptorZ->j == 0);
        REQUIRE(HICMA_descriptorZ->mt == ceil((N * 1.0) / (lts * 1.0)));
        REQUIRE(HICMA_descriptorZ->nt == 1);
        REQUIRE(HICMA_descriptorZ->lm == N);
        REQUIRE(HICMA_descriptorZ->ln == 1);
        REQUIRE(HICMA_descriptorZ->p == pGrid);
        REQUIRE(HICMA_descriptorZ->q == qGrid);

        // Descriptor Zcpy.
        REQUIRE(HICMA_descriptorZcpy->m == N);
        REQUIRE(HICMA_descriptorZcpy->n == 1);
        REQUIRE(HICMA_descriptorZcpy->mb == lts);
        REQUIRE(HICMA_descriptorZcpy->nb == lts);
        REQUIRE(HICMA_descriptorZcpy->bsiz == lts * lts);
        REQUIRE(HICMA_descriptorZcpy->i == 0);
        REQUIRE(HICMA_descriptorZcpy->j == 0);
        REQUIRE(HICMA_descriptorZcpy->mt == ceil((N * 1.0) / (lts * 1.0)));
        REQUIRE(HICMA_descriptorZcpy->nt == 1);
        REQUIRE(HICMA_descriptorZcpy->lm == N);
        REQUIRE(HICMA_descriptorZcpy->ln == 1);
        REQUIRE(HICMA_descriptorZcpy->p == pGrid);
        REQUIRE(HICMA_descriptorZcpy->q == qGrid);

        // Descriptor Determinant.
        REQUIRE(HICMA_descriptorDeterminant->m == 1);
        REQUIRE(HICMA_descriptorDeterminant->n == 1);
        REQUIRE(HICMA_descriptorDeterminant->mb == lts);
        REQUIRE(HICMA_descriptorDeterminant->nb == lts);
        REQUIRE(HICMA_descriptorDeterminant->bsiz == lts * lts);
        REQUIRE(HICMA_descriptorDeterminant->i == 0);
        REQUIRE(HICMA_descriptorDeterminant->j == 0);
        REQUIRE(HICMA_descriptorDeterminant->mt == 1);
        REQUIRE(HICMA_descriptorDeterminant->nt == 1);
        REQUIRE(HICMA_descriptorDeterminant->lm == 1);
        REQUIRE(HICMA_descriptorDeterminant->ln == 1);
        REQUIRE(HICMA_descriptorDeterminant->p == pGrid);
        REQUIRE(HICMA_descriptorDeterminant->q == qGrid);
    }
}

TEST_CASE("HiCMA Implementation TLR") {
    TEST_HICMA_DESCRIPTORS_VALUES_TLR();
}
