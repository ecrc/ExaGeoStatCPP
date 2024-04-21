
// Copyright (c) 2017-2024 King Abdullah University of Science and Technology,
// All rights reserved.
// ExaGeoStat is a software package, provided by King Abdullah University of Science and Technology (KAUST).

/**
 * @file TestUnivariateMaternDdnuNu.cpp
 * @brief Unit tests for the TestUnivariateMaternDdnuNu kernel in the ExaGeoStat software package.
 * @details This file contains Catch2 unit tests that validate the functionality of the TestUnivariateMaternDdnuNu kernel
 * in the ExaGeoStat software package. The tests cover the generation of data using this kernel with various configurations.
 * @version 1.1.0
 * @author Mahmoud ElKarargy
 * @author Sameh Abdulah
 * @date 2023-05-10
**/

#include <catch2/catch_all.hpp>
#include <api/ExaGeoStat.hpp>

using namespace std;

using namespace exageostat::api;
using namespace exageostat::common;
using namespace exageostat::configurations;

void TEST_KERNEL_GENERATION_UnivariateMaternDdnuNu() {

    SECTION("UnivariateMaternDdnuNu")
    {

        // Create a new synthetic_data_configurations object with the provided command line arguments
        Configurations synthetic_data_configurations;
        int N = 9;
        synthetic_data_configurations.SetProblemSize(N);
        synthetic_data_configurations.SetKernelName("UnivariateMaternDdnuNu");

        vector<double> initial_theta{0.1, 0.1, 0.1};
        synthetic_data_configurations.SetInitialTheta(initial_theta);

        int dts = 5;
        synthetic_data_configurations.SetDenseTileSize(dts);
        synthetic_data_configurations.SetComputation(EXACT_DENSE);
        //// TODO: Missing values in C.
    }
}

TEST_CASE("UnivariateMaternDdnuNu kernel test") {
    TEST_KERNEL_GENERATION_UnivariateMaternDdnuNu();

}