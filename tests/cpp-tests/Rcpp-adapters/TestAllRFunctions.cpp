// Copyright (c) 2017-2023 King Abdullah University of Science and Technology,
// All rights reserved.
// ExaGeoStat is a software package, provided by King Abdullah University of Science and Technology (KAUST).

/**
 * @file TestAllRFunctions.cpp
 * @brief Test suite for R/Rcpp adapters in C++.
 * @details This file contains tests for R methods adapted for use in C++ within the ExaGeoStat software package.
 * It includes tests for initializing arguments, loading data, modeling data, and predicting data using Rcpp adapters.
 * The test suite specifically verifies the integration and functionality of R methods through the Rcpp interface
 * within a C++ environment, ensuring that statistical models and algorithms are correctly initialized,
 * data is accurately loaded and processed, and predictions are properly executed.
 * @version 1.1.0
 * @author Mahmoud ElKarargy
 * @date 2024-02-09
**/

#include <catch2/catch_all.hpp>

#include <Rcpp-adapters/FunctionsAdapter.hpp>

using namespace std;

using namespace exageostat::adapters;

void TEST_ALL_R_METHODS() {
    SECTION("R METHODS"){

        auto configurations = R_InitializeArguments(16, "univariate_matern_stationary", {8, 0}, {1, 1}, 1, "exact", "double", {1, 0}, 0, 1, {1,0.1,0.5}, { {0.1,0.1,0.1}, {5,5,5}}, {-1,-1,-1}, "standard", "2D", 10, 4, {5,1,1,1,1});
        auto hardware = ExaGeoStatHardware("exact", 1, 0);
        auto data_source = new ExaGeoStatData<double>(16, "2D");
        auto exageostat_data = R_ExaGeoStatLoadData(&hardware, configurations, data_source);
        R_ExaGeoStatModelData(&hardware, configurations, exageostat_data);
        R_ExaGeoStatPredictData(&hardware, configurations, exageostat_data);

        delete exageostat_data;
        delete configurations;
    }
}

TEST_CASE("Test R/Rcpp adapters in C++") {
    TEST_ALL_R_METHODS();
}