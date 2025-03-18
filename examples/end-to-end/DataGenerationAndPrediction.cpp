
// Copyright (c) 2017-2024 King Abdullah University of Science and Technology,
// All rights reserved.
// ExaGeoStat is a software package, provided by King Abdullah University of Science and Technology (KAUST).

/**
 * @file DataGenerationAndPrediction.cpp
 * @brief This program This program either generates synthetic data using the ExaGeoStat library, or reads an CSV file containing real data, then predicts missing measurements using the ExaGeoStat library.
 * @details The program takes command line arguments to configure the data generation.
 * @version 1.1.0
 * @author Mahmoud ElKarargy
 * @date 2024-03-03
**/

#include <api/ExaGeoStat.hpp>

using namespace exageostat::api;
using namespace exageostat::configurations;

/**
 * @brief Main entry point for the Data Generation & Data Modeling program.
 * @details This function either generates synthetic data using the ExaGeoStat library, or reads an CSV file containing real data, models it, and predicts missing values.
 * @param[in] argc The number of command line arguments.
 * @param[in] argv An array of command line argument strings.
 * @return An integer indicating the success or failure of the program. A return value of 0 indicates success, while any non-zero value indicates failure.
 *
 */
int main(int argc, char **argv) {

    // Create a new configurations object.
    Configurations configurations;
    // Initialize the arguments with the provided command line arguments
    configurations.InitializeArguments(argc, argv);
    // Initialize the ExaGeoStat Hardware
    auto hardware = ExaGeoStatHardware(configurations);
    // Load data by either read from file or create synthetic data.
    std::unique_ptr<ExaGeoStatData<double>> data;
    ExaGeoStat<double>::ExaGeoStatLoadData(configurations, data);
    // Prediction module
    ExaGeoStat<double>::ExaGeoStatPrediction(configurations, data);

    return 0;
}
