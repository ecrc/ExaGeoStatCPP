
// Copyright (c) 2017-2023 King Abdullah University of Science and Technology,
// All rights reserved.
// ExaGeoStat is a software package, provided by King Abdullah University of Science and Technology (KAUST).

/**
 * @file TestUnivariateMaternNonStat.cpp
 * @brief Unit tests for the TestUnivariateMaternNonStat kernel in the ExaGeoStat software package.
 * @details This file contains Catch2 unit tests that validate the functionality of the TestUnivariateMaternNonStat kernel
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

void TEST_KERNEL_GENERATION_UnivariateMaternNonStat() {

    SECTION("UnivariateMaternNonStat")
    {

        // Create a new synthetic_data_configurations object with the provided command line arguments
        Configurations synthetic_data_configurations;

        int N = 16;
        synthetic_data_configurations.SetProblemSize(N);
        synthetic_data_configurations.SetKernelName("UnivariateMaternNonStat");
        synthetic_data_configurations.SetDimension(Dimension2D);

        vector<double> initial_theta{0.04, 1.57, 0.33, -1, 0.8, 0.1, -0.5, 0.5};
        synthetic_data_configurations.SetInitialTheta(initial_theta);

        int dts = 3;
        synthetic_data_configurations.SetDenseTileSize(dts);
        synthetic_data_configurations.SetComputation(EXACT_DENSE);
        // initialize ExaGeoStat Hardware.
        auto hardware = ExaGeoStatHardware(EXACT_DENSE, 3, 0);

        int seed = 0;
        srand(seed);
        exageostat::dataunits::ExaGeoStatData<double> data;
        exageostat::api::ExaGeoStat<double>::ExaGeoStatLoadData(hardware, synthetic_data_configurations,
                                                                data);
        auto *CHAM_descriptorZ = data.GetDescriptorData()->GetDescriptor(exageostat::common::CHAMELEON_DESCRIPTOR,
                                                                         exageostat::common::DESCRIPTOR_Z).chameleon_desc;
        auto *A = (double *) CHAM_descriptorZ->mat;
        // Define the expected output
        double expected_output_data[] = {
                -1.329875, -2.691617, 0.043584, -0.606902, -0.213657, -1.322967,
                -0.281561, 0.084918, -1.152730, -1.209422, -1.776260, -0.410195,
                0.221300, 0.454534, -0.742289, 0.135949
        };

        for (size_t i = 0; i < N; i++) {
            double diff = A[i] - expected_output_data[i];
            REQUIRE(diff == Catch::Approx(0.0).margin(1e-6));
        }
    }
}

TEST_CASE("UnivariateMaternNonStat kernel test") {
    TEST_KERNEL_GENERATION_UnivariateMaternNonStat();

}