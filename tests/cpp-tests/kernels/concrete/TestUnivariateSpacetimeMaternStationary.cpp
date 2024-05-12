
// Copyright (c) 2017-2024 King Abdullah University of Science and Technology,
// All rights reserved.
// ExaGeoStat is a software package, provided by King Abdullah University of Science and Technology (KAUST).

/**
 * @file TestUnivariateSpacetimeMaternStationary.cpp
 * @brief Unit tests for the TestUnivariateSpacetimeMaternStationary kernel in the ExaGeoStat software package.
 * @details This file contains Catch2 unit tests that validate the functionality of the TestUnivariateSpacetimeMaternStationary kernel
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
using namespace exageostat::configurations;

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

        vector<double> initial_theta{1, 1, 0.1, 0.5, 0.5, 0.1, 0};
        synthetic_data_configurations.SetInitialTheta(initial_theta);

        int dts = 2;
        synthetic_data_configurations.SetDenseTileSize(dts);
        synthetic_data_configurations.SetComputation(EXACT_DENSE);

        // initialize ExaGeoStat Hardware.
        auto hardware = ExaGeoStatHardware(EXACT_DENSE, 3, 0);

        int seed = 0;
        srand(seed);
        std::unique_ptr<ExaGeoStatData<double>> data;
        exageostat::api::ExaGeoStat<double>::ExaGeoStatLoadData(synthetic_data_configurations, data);
        auto *CHAM_descriptorZ = data->GetDescriptorData()->GetDescriptor(exageostat::common::CHAMELEON_DESCRIPTOR,
                                                                          exageostat::common::DESCRIPTOR_Z).chameleon_desc;
        auto *A = (double *) CHAM_descriptorZ->mat;
        // Define the expected output
        double expected_output_data[] = {-1.272336, -2.600097, -0.482699, -0.533521, 0.008692, -1.750492, -0.453709,
                                         0.336176, -1.573801, -1.256633, -1.817694, -0.305688, 0.627641, 0.376389,
                                         -0.939680, 0.167822, 0.514814, -1.315864, 1.884674, -0.234791};

        for (size_t i = 0; i < N * 5; i++) {
            double diff = A[i] - expected_output_data[i];
            REQUIRE(diff == Catch::Approx(0.0).margin(1e-6));
        }
    }
}


TEST_CASE("UnivariateSpacetimeMaternStationary kernel test") {
    TEST_KERNEL_GENERATION_UnivariateSpacetimeMaternStationary();

}