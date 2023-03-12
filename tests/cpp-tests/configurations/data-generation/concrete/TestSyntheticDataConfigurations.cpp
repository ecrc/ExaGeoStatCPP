
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
using namespace exageostat::dataunits;

void TEST_SYNTHETIC_CONFIGURATIONS() {

    // Init using unique ptr as when unique_ptr is destroyed/get out of scope, the resource is automatically claimed.
    unique_ptr<SyntheticDataConfigurations> syntheticDataConfigurations = make_unique<SyntheticDataConfigurations>();
    int random_number = 512;

    SECTION("Dimensions setter/getter test")
    {
        syntheticDataConfigurations->SetDimension(Dimension2D);
        REQUIRE(syntheticDataConfigurations->GetDimension() == Dimension2D);
    }SECTION("Dimensions value checker test")
    {
        REQUIRE_THROWS_WITH(
                syntheticDataConfigurations->CheckDimensionValue("4D"),
                "Invalid value for Dimension. Please use 2D, 3D or ST.");
        syntheticDataConfigurations->CheckDimensionValue("2D");
    }

    SECTION("P-GRID setter/getter test")
    {
        REQUIRE(syntheticDataConfigurations->GetPGrid() == 0);
        syntheticDataConfigurations->SetPGrid(random_number);
        REQUIRE(syntheticDataConfigurations->GetPGrid() == random_number);
    }SECTION("P-GRID value checker test")
    {
        REQUIRE_THROWS_WITH(
                syntheticDataConfigurations->CheckNumericalValue("K"),
                "Invalid value. Please use Numerical values only.");
        REQUIRE_THROWS_WITH(
                syntheticDataConfigurations->CheckNumericalValue("-100"),
                "Invalid value. Please use positive values");
        int test_nb = syntheticDataConfigurations->CheckNumericalValue("512");
        syntheticDataConfigurations->SetPGrid(test_nb);
        REQUIRE(syntheticDataConfigurations->GetPGrid() == 512);
    }
}

void TEST_DATA_CONFIGURATIONS() {
    unique_ptr<SyntheticDataConfigurations> syntheticDataConfigurations = make_unique<SyntheticDataConfigurations>();
    SECTION("Kernel setter/getter test")
    {
        REQUIRE(syntheticDataConfigurations->GetKernel() == "");
        syntheticDataConfigurations->SetKernel("univariate_matern_stationary");
        REQUIRE(syntheticDataConfigurations->GetKernel() == "univariate_matern_stationary");
    }SECTION("Kernel checker value test")
    {
        REQUIRE_THROWS_WITH(
                syntheticDataConfigurations->CheckKernelValue("100"),
                "Invalid value for Kernel. Please check manual.");
        REQUIRE_THROWS_WITH(
                syntheticDataConfigurations->CheckKernelValue("univariate_matern_dnu%"),
                "Invalid value for Kernel. Please check manual.");
        syntheticDataConfigurations->CheckKernelValue("univariate_matern_dnu");
    }
}

void TEST_CONFIGURATIONS() {
    unique_ptr<SyntheticDataConfigurations> syntheticDataConfigurations = make_unique<SyntheticDataConfigurations>();
    int random_number = 512;
    SECTION("Problem size setter/getter test")
    {
        REQUIRE(syntheticDataConfigurations->GetProblemSize() == 0);
        syntheticDataConfigurations->SetProblemSize(random_number);
        REQUIRE(syntheticDataConfigurations->GetProblemSize() == random_number);
    }SECTION("Problem size checker value test")
    {
        REQUIRE_THROWS_WITH(
                syntheticDataConfigurations->CheckNumericalValue("K"),
                "Invalid value. Please use Numerical values only.");
        REQUIRE_THROWS_WITH(
                syntheticDataConfigurations->CheckNumericalValue("-100"),
                "Invalid value. Please use positive values");
        int test_nb = syntheticDataConfigurations->CheckNumericalValue("512");
        syntheticDataConfigurations->SetProblemSize(test_nb);
        REQUIRE(syntheticDataConfigurations->GetProblemSize() == 512);
    }
}

TEST_CASE("Synthetic Data Configurations") {
    TEST_SYNTHETIC_CONFIGURATIONS();
    TEST_DATA_CONFIGURATIONS();
    TEST_CONFIGURATIONS();
}
