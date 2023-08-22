
// Copyright (c) 2017-2023 King Abdullah University of Science and Technology,
// All rights reserved.
// ExaGeoStat is a software package, provided by King Abdullah University of Science and Technology (KAUST).

/**
 * @file TestUnivariateMaternNonGaussian.cpp
 * @brief Unit tests for the TestUnivariateMaternNonGaussian kernel in the ExaGeoStat software package.
 * @details This file contains Catch2 unit tests that validate the functionality of the TestUnivariateMaternNonGaussian kernel
 * in the ExaGeoStat software package. The tests cover the generation of data using this kernel with various configurations.
 * @version 1.0.0
 * @author Sameh Abdulah
 * @author Mahmoud ElKarargy
 * @date 2023-05-10
**/

#include <catch2/catch_all.hpp>
#include <api/ExaGeoStat.hpp>

using namespace std;

using namespace exageostat::configurations;
using namespace exageostat::api;
using namespace exageostat::common;
using namespace exageostat::hardware;

void TEST_KERNEL_GENERATION_UnivariateMaternNonGaussian() {

    SECTION("UnivariateMaternNonGaussian")
    {

        // Create a new synthetic_data_configurations object with the provided command line arguments
        Configurations synthetic_data_configurations;
        int N = 16;
        synthetic_data_configurations.SetProblemSize(N);
        synthetic_data_configurations.SetKernelName("UnivariateMaternNonGaussian");

        vector<double> lb{0.01, 0.01, -5, 0.1, -2, 0};
        synthetic_data_configurations.SetLowerBounds(lb);

        vector<double> ub{15, 5, 5, 5, 2, 2};
        synthetic_data_configurations.SetUpperBounds(ub);

        vector<double> initial_theta{7.0711, 1, 0, 2, 0, 0};
        synthetic_data_configurations.SetInitialTheta(initial_theta);

#ifdef EXAGEOSTAT_USE_CHAMELEON
        int dts = 3;
        synthetic_data_configurations.SetDenseTileSize(dts);
        synthetic_data_configurations.SetComputation(EXACT_DENSE);
        // Initialise ExaGeoStat Hardware.
        auto hardware = ExaGeoStatHardware(EXACT_DENSE, 3, 0);

        int seed = 0;
        srand(seed);
        exageostat::dataunits::ExaGeoStatData<double> data(synthetic_data_configurations.GetProblemSize(),
                                                           synthetic_data_configurations.GetDimension(), hardware);
        exageostat::api::ExaGeoStat<double>::ExaGeoStatGenerateData(hardware, synthetic_data_configurations, data);
        auto *CHAM_descriptorZ = data.GetDescriptorData()->GetDescriptor(exageostat::common::CHAMELEON_DESCRIPTOR,
                                                                         exageostat::common::DESCRIPTOR_Z).chameleon_desc;
        auto *A = (double *) CHAM_descriptorZ->mat;
        // Define the expected output
        double expected_output_data[] = {
                -2.544672, -3.384362, -2.318286, -3.072392, -3.574980, -4.060970,
                -3.272706, -3.435357, -2.721393, -3.245572, -3.737585, -3.587366,
                -3.032256, -2.840226, -3.605374, -2.957129
        };

        for (size_t i = 0; i < N; i++) {
            double diff = A[i] - expected_output_data[i];
            REQUIRE(diff ==Catch::Approx(0.0).margin(1e-6));
        }

#endif
    }
}

TEST_CASE("UnivariateMaternNonGaussian kernel test") {
TEST_KERNEL_GENERATION_UnivariateMaternNonGaussian();

}