
/*
 * Copyright (c) 2017-2023 King Abdullah University of Science and Technology,
 * Copyright (C) 2023 by Brightskies inc,
 * All rights reserved.
 * ExaGeoStat is a software package, provided by King Abdullah University of Science and Technology (KAUST).
 */

/**
* @file synthetic_data.cpp
* @version 1.0.0
* @author Sameh Abdulah
* @date 2023-01-31
**/
#include <libraries/catch/catch.hpp>
#include <configurations/data-generation/concrete/SyntheticDataConfigurations.hpp>
#include <iostream>

using namespace std;
using namespace exageostat::configurations::data_configurations;
using namespace exageostat::common;

void TEST_SYNTHETIC_CONFIGURATIONS() {

    SyntheticDataConfigurations synthetic_data_configurations;
    int random_number = 512;

    SECTION("Dimensions setter/getter test")
    {
        synthetic_data_configurations.SetDimension(Dimension2D);
        REQUIRE(synthetic_data_configurations.GetDimension() == Dimension2D);
    }SECTION("Dimensions value checker test")
    {
        REQUIRE_THROWS_WITH(
                synthetic_data_configurations.CheckDimensionValue("4D"),
                "Invalid value for Dimension. Please use 2D, 3D or ST.");
        SyntheticDataConfigurations::CheckDimensionValue("2D");
    }

    SECTION("P-GRID setter/getter test")
    {
        REQUIRE(synthetic_data_configurations.GetPGrid() == 1);
        synthetic_data_configurations.SetPGrid(random_number);
        REQUIRE(synthetic_data_configurations.GetPGrid() == random_number);
    }SECTION("P-GRID value checker test")
    {
        REQUIRE_THROWS_WITH(
                synthetic_data_configurations.CheckNumericalValue("K"),
                "Invalid value. Please use Numerical values only.");
        REQUIRE_THROWS_WITH(
                synthetic_data_configurations.CheckNumericalValue("-100"),
                "Invalid value. Please use positive values");
        int test_nb = SyntheticDataConfigurations::CheckNumericalValue("512");
        synthetic_data_configurations.SetPGrid(test_nb);
        REQUIRE(synthetic_data_configurations.GetPGrid() == 512);
    }
}

void TEST_DATA_CONFIGURATIONS() {
    SyntheticDataConfigurations synthetic_data_configurations;
    SECTION("Kernel setter/getter test")
    {
        REQUIRE(synthetic_data_configurations.GetKernel().empty());
        synthetic_data_configurations.SetKernel("univariate_matern_stationary");
        REQUIRE(synthetic_data_configurations.GetKernel() == "univariate_matern_stationary");
    }SECTION("Kernel checker value test")
    {
        REQUIRE_THROWS_WITH(
                synthetic_data_configurations.CheckKernelValue("100"),
                "Invalid value for Kernel. Please check manual.");
        REQUIRE_THROWS_WITH(
                synthetic_data_configurations.CheckKernelValue("univariate_matern_dnu%"),
                "Invalid value for Kernel. Please check manual.");
        synthetic_data_configurations.CheckKernelValue("univariate_matern_dnu");
    }
}

void TEST_CONFIGURATIONS() {
    SyntheticDataConfigurations synthetic_data_configurations;
    int random_number = 512;
    SECTION("Problem size setter/getter test")
    {
        REQUIRE(synthetic_data_configurations.GetProblemSize() == 0);
        synthetic_data_configurations.SetProblemSize(random_number);
        REQUIRE(synthetic_data_configurations.GetProblemSize() == random_number);
    }SECTION("Problem size checker value test")
    {
        REQUIRE_THROWS_WITH(
                synthetic_data_configurations.CheckNumericalValue("K"),
                "Invalid value. Please use Numerical values only.");
        REQUIRE_THROWS_WITH(
                synthetic_data_configurations.CheckNumericalValue("-100"),
                "Invalid value. Please use positive values");
        int test_nb = SyntheticDataConfigurations::CheckNumericalValue("512");
        synthetic_data_configurations.SetProblemSize(test_nb);
        REQUIRE(synthetic_data_configurations.GetProblemSize() == 512);
    }
}

TEST_CASE("Synthetic Data Configurations") {
    TEST_SYNTHETIC_CONFIGURATIONS();
    TEST_DATA_CONFIGURATIONS();
    TEST_CONFIGURATIONS();
}
