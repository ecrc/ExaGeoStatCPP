
// Copyright (c) 2017-2023 King Abdullah University of Science and Technology,
// All rights reserved.
// ExaGeoStat is a software package, provided by King Abdullah University of Science and Technology (KAUST).

/**
 * @file TestUnivariateMaternDdbetaBeta.cpp
 * @brief Unit tests for the TestUnivariateMaternDdbetaBeta kernel in the ExaGeoStat software package.
 * @details This file contains Catch2 unit tests that validate the functionality of the TestUnivariateMaternDdbetaBeta kernel
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

void TEST_KERNEL_GENERATION_UnivariateMaternDdbetaBeta() {

    SECTION("UnivariateMaternDdbetaBeta")
    {

        // Create a new synthetic_data_configurations object with the provided command line arguments
        Configurations synthetic_data_configurations;
        int N = 16;
        synthetic_data_configurations.SetProblemSize(N);
        synthetic_data_configurations.SetKernelName("UnivariateMaternDdbetaBeta");

        vector<double> lb{0.1, 0.1, 0.1};
        synthetic_data_configurations.SetLowerBounds(lb);

        vector<double> ub{5, 5, 5};
        synthetic_data_configurations.SetUpperBounds(ub);

        vector<double> initial_theta{0.1, 0.1, 0.5};
        synthetic_data_configurations.SetInitialTheta(initial_theta);

#ifdef EXAGEOSTAT_USE_CHAMELEON
        int dts = 8;
        synthetic_data_configurations.SetDenseTileSize(dts);
        synthetic_data_configurations.SetComputation(EXACT_DENSE);

        //// TODO: The values are missing in C
#endif
    }
}

TEST_CASE("UnivariateMaternDdbetaBeta kernel test") {
TEST_KERNEL_GENERATION_UnivariateMaternDdbetaBeta();

}