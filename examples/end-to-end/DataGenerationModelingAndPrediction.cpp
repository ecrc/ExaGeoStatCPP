
// Copyright (c) 2017-2023 King Abdullah University of Science and Technology,
// All rights reserved.
// ExaGeoStat is a software package, provided by King Abdullah University of Science and Technology (KAUST).

/**
 * @file DataGenerationModelingAndPrediction.cpp
 * @brief This program This program either generates synthetic data using the ExaGeoStat library, or reads an CSV file containing real data, performs data modeling on loaded data, then predicts missing measurements using the ExaGeoStat library.
 * @details The program takes command line arguments to configure the data generation.
 * @version 1.0.0
 * @author Mahmoud ElKarargy
 * @date 2023-06-21
**/

#include <common/Utils.hpp>
#include <api/ExaGeoStat.hpp>

using namespace exageostat::configurations;
using namespace exageostat::api;
using namespace exageostat::hardware;
using namespace exageostat::dataunits;

/**
 * @brief Main entry point for the Data Generation & Data Modeling program.
 * @details This function either generates synthetic data using the ExaGeoStat library, or reads an CSV file containing real data, models it, and predicts missing values.
 * @param[in] argc The number of command line arguments.
 * @param[in] argv An array of command line argument strings.
 * @return An integer indicating the success or failure of the program.
 */
int main(int argc, char **argv) {

    // Create a new configurations object.
    Configurations configurations;
    //  Initialize the arguments with the provided command line arguments
    configurations.InitializeArguments(argc, argv);
    LOGGER("** Initialise ExaGeoStat hardware **")
    auto hardware = ExaGeoStatHardware(configurations.GetComputation(), configurations.GetCoresNumber(),
                                       configurations.GetGPUsNumbers());
    LOGGER("** Create ExaGeoStat data **")
    ExaGeoStatData<double> data;
    LOGGER("** ExaGeoStat data generation **")
    ExaGeoStat<double>::ExaGeoStatLoadData(hardware, configurations, data);
    LOGGER("** ExaGeoStat data Modeling **")
    ExaGeoStat<double>::ExaGeoStatDataModeling(hardware, configurations, data);
    LOGGER("** ExaGeoStat data Prediction **")
    ExaGeoStat<double>::ExaGeoStatPrediction(hardware, configurations, data);
    LOGGER("** All example stages have been completed successfully ** ")
    return 0;
}
