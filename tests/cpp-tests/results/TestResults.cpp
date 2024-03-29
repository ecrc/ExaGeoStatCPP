// Copyright (c) 2017-2023 King Abdullah University of Science and Technology,
// All rights reserved.
// ExaGeoStat is a software package, provided by King Abdullah University of Science and Technology (KAUST).

/**
 * @file TestResults.cpp
 * @brief Unit tests for the Results class in the ExaGeoStat software package.
 * @details This file contains Catch2 unit tests that validate the functionality of the Results class
 * @version 1.1.0
 * @author Mahmoud ElKarargy
 * @date 2024-01-24
**/

#include <catch2/catch_all.hpp>

#include <results/Results.hpp>

using namespace exageostat::results;

void TEST_SINGLETON_RESULTS() {
    //Test singleton instance of results
    auto results1 = Results::GetInstance();
    auto results2 = Results::GetInstance();
    REQUIRE(results1 == results2);
}

void TEST_SETTERS_AND_GETTERS() {
    
    auto results_instacne = Results::GetInstance();
    
    SECTION("Total Modeling Execution Time Setter/Getter") {
        results_instacne->SetTotalModelingExecutionTime(1.0);
        REQUIRE(results_instacne->GetTotalModelingExecutionTime() == 1.0);

    }SECTION("Total Modeling Flops Setter/Getter") {
        results_instacne->SetTotalModelingFlops(1.0);
        REQUIRE(results_instacne->GetTotalModelingFlops() == 1.0);

    }SECTION("Avg Modeling Execution Time Setter/Getter") {
        results_instacne->SetMLEIterations(0);
        REQUIRE_THROWS(results_instacne->GetAverageModelingExecutionTime());

        results_instacne->SetMLEIterations(2);
        results_instacne->SetTotalModelingExecutionTime(4.0);
        REQUIRE(results_instacne->GetAverageModelingExecutionTime() == 2.0);

    }SECTION("Avg Modeling Flops Setter/Getter") {
        results_instacne->SetMLEIterations(0);
        REQUIRE_THROWS(results_instacne->GetAverageModelingFlops());

        results_instacne->SetMLEIterations(2);
        results_instacne->SetTotalModelingFlops(4.0);
        REQUIRE(results_instacne->GetAverageModelingFlops() == 2.0);
    }
}

TEST_CASE("Test Results") {
    TEST_SINGLETON_RESULTS();
    TEST_SETTERS_AND_GETTERS();
}
