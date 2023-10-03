
// Copyright (c) 2017-2023 King Abdullah University of Science and Technology,
// All rights reserved.
// ExaGeoStat is a software package, provided by King Abdullah University of Science and Technology (KAUST).

/**
 * @file TestBivariateSpacetimeMaternStationary.cpp
 * @brief Unit tests for the TestBivariateSpacetimeMaternStationary kernel in the ExaGeoStat software package.
 * @details This file contains Catch2 unit tests that validate the functionality of the TestBivariateSpacetimeMaternStationary kernel
 * in the ExaGeoStat software package. The tests cover the generation of data using this kernel with various configurations.
 * @version 1.0.0
 * @author Mahmoud ElKarargy
 * @author Sameh Abdulah
 * @date 2023-05-10
**/

#include <catch2/catch_all.hpp>
#include <api/ExaGeoStat.hpp>
#include <hardware/ExaGeoStatHardware.hpp>

using namespace std;

using namespace exageostat::configurations;
using namespace exageostat::api;
using namespace exageostat::common;
using namespace exageostat::hardware;

void TEST_KERNEL_GENERATION_BivariateSpacetimeMaternStationary() {

    SECTION("BivariateSpacetimeMaternStationary")
    {

        // Create a new synthetic_data_configurations object with the provided command line arguments
        Configurations synthetic_data_configurations;
        int N = 4;
        synthetic_data_configurations.SetProblemSize(N);
        synthetic_data_configurations.SetKernelName("BivariateSpacetimeMaternStationary");

        vector<double> initial_theta{1, 1, 1, 0.1, 0.5, 1, 1.5, 0.1, 0.1, 0.1};
        synthetic_data_configurations.SetInitialTheta(initial_theta);

        synthetic_data_configurations.SetDimension(DimensionST);
        synthetic_data_configurations.SetTimeSlot(5);
        int dts = 4;
        synthetic_data_configurations.SetDenseTileSize(dts);
        synthetic_data_configurations.SetComputation(EXACT_DENSE);

        // initialize ExaGeoStat Hardware.
        auto hardware = ExaGeoStatHardware(EXACT_DENSE, 3, 0);

        int seed = 0;
        srand(seed);
        exageostat::dataunits::ExaGeoStatData<double> data(synthetic_data_configurations.GetProblemSize(),
                                                           synthetic_data_configurations.GetDimension());
        exageostat::api::ExaGeoStat<double>::ExaGeoStatGenerateData(hardware, synthetic_data_configurations, data);
        auto *CHAM_descriptorZ = data.GetDescriptorData()->GetDescriptor(exageostat::common::CHAMELEON_DESCRIPTOR,
                                                                         exageostat::common::DESCRIPTOR_Z).chameleon_desc;
        auto *A = (double *) CHAM_descriptorZ->mat;
        // Define the expected output
        double expected_output_data[] = {
                -1.272336, -2.516013, -1.182511, -2.512958, -1.203093, -2.548053, -1.213645, -2.511269, -1.807922,
                -2.992990, -2.072972, -3.001448, -1.525724, -3.004085, -1.997874, -2.993439, -1.256771, -2.832258,
                -1.022281, -2.827732, -0.924332, -2.856855, -1.363071, -2.832165, -1.680104, -3.400946, -1.685995,
                -3.414854, -1.318928, -3.440336, -1.944915, -3.428945, -1.300296, -3.367150, -1.572847, -3.392844,
                -1.126261, -3.392112, -1.682665, -3.394862
        };

        for (size_t i = 0; i < N; i++) {
            double diff = A[i] - expected_output_data[i];
            REQUIRE(diff == Catch::Approx(0.0).margin(1e-6));
        }
    }
}

TEST_CASE("BivariateSpacetimeMaternStationary kernel test") {
    TEST_KERNEL_GENERATION_BivariateSpacetimeMaternStationary();

}