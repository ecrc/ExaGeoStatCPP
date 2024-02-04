
// Copyright (c) 2017-2023 King Abdullah University of Science and Technology,
// All rights reserved.
// ExaGeoStat is a software package, provided by King Abdullah University of Science and Technology (KAUST).

/**
 * @file TestBivariateMaternParsimonious.cpp
 * @brief Unit tests for the TestBivariateMaternParsimonious kernel in the ExaGeoStat software package.
 * @details This file contains Catch2 unit tests that validate the functionality of the TestBivariateMaternParsimonious kernel
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

void TEST_KERNEL_GENERATION_BivariateMaternParsimonious() {

    SECTION("BivariateMaternParsimonious")
    {
        // Create a new synthetic_data_configurations object with the provided command line arguments
        Configurations synthetic_data_configurations;

        int N = 16;

        synthetic_data_configurations.SetProblemSize(N);
        synthetic_data_configurations.SetKernelName("BivariateMaternParsimonious");

        vector<double> initial_theta{1, 1, 0.1, 0.5, 0.5, 0.1};
        synthetic_data_configurations.SetInitialTheta(initial_theta);


        int dts = 8;
        synthetic_data_configurations.SetDenseTileSize(dts);
        synthetic_data_configurations.SetComputation(EXACT_DENSE);

        // initialize ExaGeoStat Hardware.
        auto hardware = ExaGeoStatHardware(EXACT_DENSE, 3, 0);

        int seed = 0;
        srand(seed);
        std::unique_ptr<ExaGeoStatData<double>> data;
        exageostat::api::ExaGeoStat<double>::ExaGeoStatLoadData(hardware, synthetic_data_configurations,
                                                                data);
        auto *CHAM_descriptorZ = data->GetDescriptorData()->GetDescriptor(exageostat::common::CHAMELEON_DESCRIPTOR,
                                                                         exageostat::common::DESCRIPTOR_Z).chameleon_desc;
        auto *A = (double *) CHAM_descriptorZ->mat;
        // Define the expected output
        double expected_output_data[] = {-1.272336, -2.466950, 0.294719, -0.605327, 0.386028, -1.598090, 0.278897,
                                         0.489645, -1.508498, -0.965737, -1.884671, -0.058567, 1.024710, 0.598136,
                                         -1.257452, 0.124507};

        for (size_t i = 0; i < N; i++) {
            double diff = A[i] - expected_output_data[i];
            REQUIRE(diff == Catch::Approx(0.0).margin(1e-6));
        }
    }
}

TEST_CASE("Bivariate Matern Parsimonious kernel test") {
    TEST_KERNEL_GENERATION_BivariateMaternParsimonious();

}