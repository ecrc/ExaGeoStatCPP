
// Copyright (c) 2017-2024 King Abdullah University of Science and Technology,
// All rights reserved.
// ExaGeoStat is a software package, provided by King Abdullah University of Science and Technology (KAUST).

/**
 * @file TestUnivariateMaternDbeta.cpp
 * @brief Unit tests for the TestUnivariateMaternDbeta kernel in the ExaGeoStat software package.
 * @details This file contains Catch2 unit tests that validate the functionality of the TestUnivariateMaternDbeta kernel
 * in the ExaGeoStat software package. The tests cover the generation of data using this kernel with various configurations.
 * @version 1.1.0
 * @author Mahmoud ElKarargy
 * @author Sameh Abdulah
 * @date 2023-05-10
**/

#include <catch2/catch_all.hpp>
#include <api/ExaGeoStat.hpp>

using namespace std;

using namespace exageostat::common;
using namespace exageostat::configurations;

void TEST_KERNEL_GENERATION_UnivariateMaternDbeta() {

    SECTION("UnivariateMaternDbeta")
    {
        // Create a new synthetic_data_configurations object with the provided command line arguments
        Configurations synthetic_data_configurations;
        int N = 9;
        synthetic_data_configurations.SetSeed(0);
        synthetic_data_configurations.SetProblemSize(N);
        synthetic_data_configurations.SetKernelName("UnivariateMaternDbeta");
        int dts = 5;
        synthetic_data_configurations.SetDenseTileSize(dts);
        synthetic_data_configurations.SetComputation(EXACT_DENSE);

        vector<double> initial_theta{0.1, 0.1, 0.1};
        synthetic_data_configurations.SetInitialTheta(initial_theta);

        //// TODO: values are missing in the C version.

    }
}

TEST_CASE("UnivariateMaternDbeta kernel test") {
    TEST_KERNEL_GENERATION_UnivariateMaternDbeta();

}