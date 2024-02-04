
// Copyright (c) 2017-2023 King Abdullah University of Science and Technology,
// All rights reserved.
// ExaGeoStat is a software package, provided by King Abdullah University of Science and Technology (KAUST).

/**
 * @file TestBivariateSpacetimeMaternStationary.cpp
 * @brief Unit tests for the TestBivariateSpacetimeMaternStationary kernel in the ExaGeoStat software package.
 * @details This file contains Catch2 unit tests that validate the functionality of the TestBivariateSpacetimeMaternStationary kernel
 * in the ExaGeoStat software package. The tests cover the generation of data using this kernel with various configurations.
 * @version 1.1.0
 * @author Mahmoud ElKarargy
 * @author Sameh Abdulah
 * @date 2023-05-10
**/

#include <catch2/catch_all.hpp>
#include <api/ExaGeoStat.hpp>
#include <hardware/ExaGeoStatHardware.hpp>

using namespace std;

using namespace exageostat::api;
using namespace exageostat::common;

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
        synthetic_data_configurations.SetTimeSlot(3);
        int dts = 4;
        synthetic_data_configurations.SetDenseTileSize(dts);
        synthetic_data_configurations.SetComputation(EXACT_DENSE);

        // initialize ExaGeoStat Hardware.
        auto hardware = ExaGeoStatHardware(EXACT_DENSE, 3, 0);

        int seed = 0;
        srand(seed);
        std::unique_ptr<exageostat::dataunits::ExaGeoStatData<double>> data;
        exageostat::api::ExaGeoStat<double>::ExaGeoStatLoadData(hardware, synthetic_data_configurations,
                                                                data);
        auto *CHAM_descriptorZ = data->GetDescriptorData()->GetDescriptor(exageostat::common::CHAMELEON_DESCRIPTOR,
                                                                          exageostat::common::DESCRIPTOR_Z).chameleon_desc;
        auto *A = (double *) CHAM_descriptorZ->mat;
        // Define the expected output
        double expected_output_data[] = {-1.272336, -2.516013, -1.171584, -2.513616, -1.168821, -2.531118, -1.200293,
                                         -1.960279, -1.919141, -2.005663, -2.279083, -2.006927, -0.807580, -1.719351,
                                         -1.532754, -1.733719, -0.935139, -1.752957, 0.125841, -1.722780, -0.162095,
                                         -1.738346, -0.866950, -1.753649};

        for (size_t i = 0; i < N * 3 * 2; i++) {
            double diff = A[i] - expected_output_data[i];
            REQUIRE(diff == Catch::Approx(0.0).margin(1e-6));
        }
    }
}

TEST_CASE("BivariateSpacetimeMaternStationary kernel test") {
    TEST_KERNEL_GENERATION_BivariateSpacetimeMaternStationary();

}