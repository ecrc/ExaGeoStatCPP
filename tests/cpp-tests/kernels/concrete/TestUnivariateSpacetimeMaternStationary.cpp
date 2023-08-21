
// Copyright (c) 2017-2023 King Abdullah University of Science and Technology,
// All rights reserved.
// ExaGeoStat is a software package, provided by King Abdullah University of Science and Technology (KAUST).

/**
 * @file TestUnivariateSpacetimeMaternStationary.cpp
 * @brief Unit tests for the TestUnivariateSpacetimeMaternStationary kernel in the ExaGeoStat software package.
 * @details This file contains Catch2 unit tests that validate the functionality of the TestUnivariateSpacetimeMaternStationary kernel
 * in the ExaGeoStat software package. The tests cover the generation of data using this kernel with various configurations.
 * @version 1.0.0
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

void TEST_KERNEL_GENERATION_UnivariateSpacetimeMaternStationary() {

    SECTION("UnivariateSpacetimeMaternStationary")
    {

        // Create a new synthetic_data_configurations object with the provided command line arguments
        Configurations synthetic_data_configurations;
        int N = 4;
        synthetic_data_configurations.SetProblemSize(N);
        synthetic_data_configurations.SetKernelName("UnivariateSpacetimeMaternStationary");
        synthetic_data_configurations.SetDimension(exageostat::common::DimensionST);
        synthetic_data_configurations.SetTimeSlot(5);

        vector<double> lb{0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0};
        synthetic_data_configurations.SetLowerBounds(lb);

        vector<double> ub{5, 5, 5, 5, 5, 5, 0};
        synthetic_data_configurations.SetUpperBounds(ub);

        vector<double> initial_theta{1, 1, 0.1, 0.5, 0.5, 0.1, 0};
        synthetic_data_configurations.SetInitialTheta(initial_theta);

        vector<double> target_theta{-1, -1, -1, -1, -1, -1, 0};
        synthetic_data_configurations.SetTargetTheta(target_theta);

#ifdef EXAGEOSTAT_USE_CHAMELEON
        int dts = 2;
        synthetic_data_configurations.SetDenseTileSize(dts);
        synthetic_data_configurations.SetComputation(EXACT_DENSE);

        // Initialise ExaGeoStat Hardware.
        auto hardware = ExaGeoStatHardware(EXACT_DENSE, 3, 0);

        int seed = 0;
        srand(seed);
        exageostat::dataunits::ExaGeoStatData<double> data(
                synthetic_data_configurations.GetProblemSize() * synthetic_data_configurations.GetTimeSlot(),
                synthetic_data_configurations.GetDimension(), hardware);
        exageostat::api::ExaGeoStat<double>::ExaGeoStatGenerateData(hardware, synthetic_data_configurations, data);
        auto *CHAM_descriptorZ = data.GetDescriptorData()->GetDescriptor(exageostat::common::CHAMELEON_DESCRIPTOR,
                                                                         exageostat::common::DESCRIPTOR_Z).chameleon_desc;
        auto *A = (double *) CHAM_descriptorZ->mat;
        // Define the expected output
        double expected_output_data[] = {
                -1.272336, -2.499643, 0.241533, -0.680865,
                -0.966917, -2.919659, 0.289843, -0.387418,
                -1.633527, -2.982286, -0.534392, -0.453970,
                -0.907212, -2.455697, -0.617580, -0.289998,
                -0.755687, -3.013465, 0.547427, -0.665510
        };

        for (size_t i = 0; i < N; i++) {
            double diff = A[i] - expected_output_data[i];
            REQUIRE(diff ==Catch::Approx(0.0).margin(1e-6));
        }
#endif
    }
}


TEST_CASE("UnivariateSpacetimeMaternStationary kernel test") {
TEST_KERNEL_GENERATION_UnivariateSpacetimeMaternStationary();

}