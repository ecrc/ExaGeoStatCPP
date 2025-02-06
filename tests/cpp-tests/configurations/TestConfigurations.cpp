
// Copyright (c) 2017-2024 King Abdullah University of Science and Technology,
// All rights reserved.
// ExaGeoStat is a software package, provided by King Abdullah University of Science and Technology (KAUST).

/**
 * @file synthetic_data.cpp
 * @brief Unit tests for the Configurations class in the ExaGeoStat software package.
 * @details This file contains Catch2 unit tests that validate the functionality of the Configurations class
 * in the ExaGeoStat software package. The tests cover various setters, getters, and value checks
 * for configuration parameters such as dimensions, P-GRID, kernel name, problem size, precision, and more.
 * Additionally, the tests include a copy-constructor test for the Configurations class.
 * @version 1.1.0
 * @author Mahmoud ElKarargy
 * @date 2023-01-31
**/

#include <catch2/catch_all.hpp>
#include <configurations/Configurations.hpp>
#include <configurations/Validator.hpp>

using namespace std;

using namespace exageostat::common;
using namespace exageostat::configurations;
using namespace exageostat::configurations::validator;

void TEST_ARGUMENT_INITIALIZATION() {

    const int argc = 17;
    char *argv[] = {
            const_cast<char *>("program_name"),
            const_cast<char *>("--N=16"),
            const_cast<char *>("--dts=8"),
            const_cast<char *>("--kernel=univariate_matern_stationary"),
            const_cast<char *>("--computation=exact"),
            const_cast<char *>("--precision=double"),
            const_cast<char *>("--initial_theta=1:0.1:0.5"),
            const_cast<char *>("--ub=5:5:5"),
            const_cast<char *>("--lb=0.1:0.1:0.1"),
            const_cast<char *>("--max_mle_iterations=5"),
            const_cast<char *>("--tolerance=8"),
            const_cast<char *>("--ZMiss=6"),
            const_cast<char *>("--mspe"),
            const_cast<char *>("--idw"),
            const_cast<char *>("--mloe-mmom"),
            const_cast<char *>("--fisher"),
            const_cast<char *>("--data_path=./dummy-path")
    };

    Configurations configurations;
    configurations.InitializeArguments(argc, argv);

    REQUIRE(configurations.GetProblemSize() == 16);
    REQUIRE(configurations.GetKernelName() == "UnivariateMaternStationary");
    REQUIRE(configurations.GetDenseTileSize() == 8);
    REQUIRE(configurations.GetPrecision() == DOUBLE);

    REQUIRE(configurations.GetDataPath() == string("./dummy-path"));

    REQUIRE(configurations.GetMaxMleIterations() == 5);
    REQUIRE(configurations.GetTolerance() == pow(10, -8));

    REQUIRE(configurations.GetIsMSPE() == true);
    REQUIRE(configurations.GetIsIDW() == true);
    REQUIRE(configurations.GetIsFisher() == true);
    REQUIRE(configurations.GetIsMLOEMMOM() == true);
    REQUIRE(configurations.GetUnknownObservationsNb() == 6);

}


void TEST_ARGUMENT_INITIALIZATION_PARSEC() {

    const int argc = 13;
    char *argv[] = {
            const_cast<char *>("program_name"),
            const_cast<char *>("--N=16"),
            const_cast<char *>("--dts=8"),
            const_cast<char *>("--precision=double"),
            const_cast<char *>("--band_dense=100"),
            const_cast<char *>("--objects-number=72"),
            const_cast<char *>("--adaptive_decision=1"),
            const_cast<char *>("--add_diagonal=10"),
            const_cast<char *>("--file_time_slot=1"),
            const_cast<char *>("--file-number=1"),
            const_cast<char *>("--enable-inverse"),
            const_cast<char *>("--mpiio"),
            const_cast<char *>("--data_path=./dummy-path"),
    };

    Configurations configurations;
    // Initialize configuration dictionary with only common arguments
    configurations.InitializeArguments(argc, argv);

    REQUIRE(configurations.GetProblemSize() == 16);
    REQUIRE(configurations.GetDenseTileSize() == 8);
    REQUIRE(configurations.GetPrecision() == DOUBLE);

    // Check Hicma-Parsec parameters
    REQUIRE(configurations.GetDenseBandDP() == 100);
    REQUIRE(configurations.GetObjectsNumber() == 72);
    REQUIRE(configurations.GetAdaptiveDecision() == 1);
    REQUIRE(configurations.GetDiagonalAddition() == 10);
    REQUIRE(configurations.GetTimeSlotPerFile() == 1);
    REQUIRE(configurations.GetFileNumber() == 1);
    REQUIRE(configurations.GetEnableInverse() == true);
    REQUIRE(configurations.GetMPIIO() == true);

    REQUIRE(configurations.GetDataPath() == string("./dummy-path"));
}


void TEST_SYNTHETIC_CONFIGURATIONS() {

    Configurations synthetic_data_configurations;
    int random_number = 512;

    SECTION("Dimensions setter/getter test")
    {
        synthetic_data_configurations.SetDimension(Dimension2D);
        REQUIRE(synthetic_data_configurations.GetDimension() == Dimension2D);
    }SECTION("Dimensions value checker test")
    {
        REQUIRE_THROWS_WITH(
                Validator::CheckDimensionValue("4D"),
                "Invalid value for Dimension. Please use 2D, 3D or ST.");
        Validator::CheckDimensionValue("2D");
    }

    SECTION("P-GRID setter/getter test")
    {
        REQUIRE(synthetic_data_configurations.GetPGrid() == 1);
        synthetic_data_configurations.SetPGrid(random_number);
        REQUIRE(synthetic_data_configurations.GetPGrid() == random_number);
    }SECTION("P-GRID value checker test")
    {
        REQUIRE_THROWS_WITH(
                Validator::CheckNumericalValue("K"),
                "Invalid value. Please use Numerical values only.");
        REQUIRE_THROWS_WITH(
                Validator::CheckNumericalValue("-100"),
                "Invalid value. Please use positive values");
        int test_nb = Validator::CheckNumericalValue("512");
        synthetic_data_configurations.SetPGrid(test_nb);
        REQUIRE(synthetic_data_configurations.GetPGrid() == 512);
    }SECTION("Kernel setter/getter test")
    {
        synthetic_data_configurations.SetKernelName("univariate_matern_stationary");
        REQUIRE(synthetic_data_configurations.GetKernelName() == "univariate_matern_stationary");
    }SECTION("Kernel checker value test")
    {
        REQUIRE_THROWS_WITH(
                Validator::CheckKernelValue("100"),
                "Invalid value for Kernel. Please check manual.");
        REQUIRE_THROWS_WITH(
                Validator::CheckKernelValue("univariate_matern_dnu%"),
                "Invalid value for Kernel. Please check manual.");
        Validator::CheckKernelValue("univariate_matern_dnu");
    }SECTION("Problem size setter/getter test")
    {
        synthetic_data_configurations.SetProblemSize(random_number);
        REQUIRE(synthetic_data_configurations.GetProblemSize() == random_number);
    }SECTION("Problem size checker value test")
    {
        REQUIRE_THROWS_WITH(
                Validator::CheckNumericalValue("K"),
                "Invalid value. Please use Numerical values only.");
        REQUIRE_THROWS_WITH(
                Validator::CheckNumericalValue("-100"),
                "Invalid value. Please use positive values");
        int test_nb = Validator::CheckNumericalValue("512");
        synthetic_data_configurations.SetProblemSize(test_nb);
        REQUIRE(synthetic_data_configurations.GetProblemSize() == 512);
    }
}

void TEST_COPY_CONSTRUCTOR() {
    SECTION("copy-constructor")
    {
        Configurations synthetic_data_configurations;
        synthetic_data_configurations.SetProblemSize(10);
        synthetic_data_configurations.SetKernelName("BivariateSpacetimeMaternStationary");
        synthetic_data_configurations.SetPrecision(MIXED);
        synthetic_data_configurations.SetLoggerPath("any/path");
        vector<double> lb{0.1, 0.1, 0.1};
        synthetic_data_configurations.SetLowerBounds(lb);
        Configurations copied_synthetic_conf(synthetic_data_configurations);
        REQUIRE(synthetic_data_configurations.GetProblemSize() == copied_synthetic_conf.GetProblemSize());
        REQUIRE(copied_synthetic_conf.GetProblemSize() == 10);
        REQUIRE(synthetic_data_configurations.GetKernelName() == copied_synthetic_conf.GetKernelName());
        REQUIRE(copied_synthetic_conf.GetKernelName() == "BivariateSpacetimeMaternStationary");
        REQUIRE(synthetic_data_configurations.GetLowerBounds() == copied_synthetic_conf.GetLowerBounds());
        REQUIRE(synthetic_data_configurations.GetPrecision() == copied_synthetic_conf.GetPrecision());
        REQUIRE(copied_synthetic_conf.GetLoggerPath() == "any/path");
        REQUIRE(synthetic_data_configurations.GetLoggerPath() == copied_synthetic_conf.GetLoggerPath());
    }
}

TEST_CASE("Configurations Tests") {
    TEST_SYNTHETIC_CONFIGURATIONS();
    TEST_COPY_CONSTRUCTOR();
    TEST_ARGUMENT_INITIALIZATION();
#if !DEFAULT_RUNTIME
    TEST_ARGUMENT_INITIALIZATION_PARSEC();
#endif
}
