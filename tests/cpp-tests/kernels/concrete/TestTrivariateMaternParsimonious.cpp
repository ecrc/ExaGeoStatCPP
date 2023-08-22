
// Copyright (c) 2017-2023 King Abdullah University of Science and Technology,
// All rights reserved.
// ExaGeoStat is a software package, provided by King Abdullah University of Science and Technology (KAUST).

/**
 * @file TestTrivariateMaternParsimonious.cpp
 * @brief Unit tests for the TestTrivariateMaternParsimonious kernel in the ExaGeoStat software package.
 * @details This file contains Catch2 unit tests that validate the functionality of the TestTrivariateMaternParsimonious kernel
 * in the ExaGeoStat software package. The tests cover the generation of data using this kernel with various configurations.
 * @version 1.0.0
 * @author Sameh Abdulah
 * @author Mahmoud ElKarargy
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

void TEST_KERNEL_GENERATION_TrivariateMaternParsimonious() {

    SECTION("TrivariateMaternParsimonious")
    {

        // Create a new synthetic_data_configurations object with the provided command line arguments
        Configurations synthetic_data_configurations;
        int N = 16;
        synthetic_data_configurations.SetProblemSize(N);
        synthetic_data_configurations.SetKernelName("TrivariateMaternParsimonious");

        vector<double> lb{0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1};
        synthetic_data_configurations.SetLowerBounds(lb);

        vector<double> ub{5, 5, 5, 5, 1, 1, 5, 5, 5, 5};
        synthetic_data_configurations.SetUpperBounds(ub);

        vector<double> initial_theta{1, 1, 1, 0.1, 0.5, 1, 1.5, 0.1, 0.1, 0};
        synthetic_data_configurations.SetInitialTheta(initial_theta);

#ifdef EXAGEOSTAT_USE_CHAMELEON
        int dts = 3;
        synthetic_data_configurations.SetDenseTileSize(dts);
        synthetic_data_configurations.SetComputation(EXACT_DENSE);
        // Initialise ExaGeoStat Hardware.
        auto hardware = ExaGeoStatHardware(EXACT_DENSE, 3, 0);

        int seed = 0;
        srand(seed);
                exageostat::dataunits::ExaGeoStatData<double> data(synthetic_data_configurations.GetProblemSize(), synthetic_data_configurations.GetDimension(), hardware);
        exageostat::api::ExaGeoStat<double>::ExaGeoStatGenerateData(hardware, synthetic_data_configurations, data);
        auto *CHAM_descriptorZ = data.GetDescriptorData()->GetDescriptor(exageostat::common::CHAMELEON_DESCRIPTOR,
                                                                         exageostat::common::DESCRIPTOR_Z).chameleon_desc;
        auto *A = (double *) CHAM_descriptorZ->mat;
        // Define the expected output
        double expected_output_data[] = {
                -1.272336, -2.460986, 0.530374, -0.388094, -0.720870, -0.938307, 0.168800, 0.588775, -1.514807,
                -0.834397, -1.520086, -0.643174, 1.076160, 0.491034, -1.035977, 0.489992, 0.749032, -1.720433,
                2.170364, -0.819514, 0.479739, -0.636219, -0.277047, -0.259929, -0.819523, -1.586966, -0.788818,
                -1.333566, -0.251557, -2.193209, -0.584865, -1.211615, 0.285566, -2.408562, -0.449549, -1.283447,
                0.638515, 0.422009, -0.434735, -0.542449, 1.873117, -2.245733, -0.493981, -0.039696, -1.486225,
                0.988614, 1.305675, -2.858326
        };

        for (size_t i = 0; i < N; i++) {
            double diff = A[i] - expected_output_data[i];
            REQUIRE(diff ==Catch::Approx(0.0).margin(1e-6));
        }


#endif
    }
}

TEST_CASE("TrivariateMaternParsimonious kernel test") {
TEST_KERNEL_GENERATION_TrivariateMaternParsimonious();

}