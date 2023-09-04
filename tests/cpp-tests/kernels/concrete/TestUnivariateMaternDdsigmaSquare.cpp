
// Copyright (c) 2017-2023 King Abdullah University of Science and Technology,
// All rights reserved.
// ExaGeoStat is a software package, provided by King Abdullah University of Science and Technology (KAUST).

/**
 * @file TestUnivariateMaternDdsigmaSquare.cpp
 * @brief Unit tests for the TestUnivariateMaternDdsigmaSquare kernel in the ExaGeoStat software package.
 * @details This file contains Catch2 unit tests that validate the functionality of the TestUnivariateMaternDdsigmaSquare kernel
 * in the ExaGeoStat software package. The tests cover the generation of data using this kernel with various configurations.
 * @version 1.0.0
 * @author Mahmoud ElKarargy
 * @author Sameh Abdulah
 * @date 2023-05-10
**/

#include <catch2/catch_all.hpp>
#include <api/ExaGeoStat.hpp>

using namespace exageostat::configurations;
using namespace exageostat::common;
using namespace exageostat::api;

using namespace std;

void TEST_KERNEL_GENERATION_UnivariateMaternDdsigmaSquare() {

    SECTION("UnivariateMaternDdsigmaSquare")
    {
        // Create a new synthetic_data_configurations object with the provided command line arguments
        Configurations synthetic_data_configurations;
        int N = 9;
        synthetic_data_configurations.SetProblemSize(N);
        synthetic_data_configurations.SetKernelName("UnivariateMaternDdsigmaSquare");
        synthetic_data_configurations.SetDimension(Dimension2D);

        vector<double> initial_theta{0.1, 0.1, 0.1};
        synthetic_data_configurations.SetInitialTheta(initial_theta);

#ifdef EXAGEOSTAT_USE_CHAMELEON
        int dts = 5;
        synthetic_data_configurations.SetDenseTileSize(dts);
        synthetic_data_configurations.SetComputation(EXACT_DENSE);

        //// TODO: values are missing in the C version.
#endif
    }
}

TEST_CASE("Univariate Matern DdsigmaSquare kernel test") {
    TEST_KERNEL_GENERATION_UnivariateMaternDdsigmaSquare();

}