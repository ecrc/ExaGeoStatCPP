
// Copyright (c) 2017-2023 King Abdullah University of Science and Technology,
// All rights reserved.
// ExaGeoStat is a software package, provided by King Abdullah University of Science and Technology (KAUST).

/**
 * @file DataGeneration.cpp
 * @brief This program generates synthetic data using the ExaGeoStat library.
 * @details The program takes command line arguments to configure the data generation.
 * @version 1.0.0
 * @author Sameh Abdulah
 * @author Mahmoud ElKarargy
 * @date 2023-05-30
**/

#include <api/ExaGeoStat.hpp>
#include <configurations/Configurations.hpp>
#include <hardware/ExaGeoStatHardware.hpp>

using namespace std;

using namespace exageostat::configurations;
using namespace exageostat::api;
using namespace exageostat::common;
using namespace exageostat::hardware;

/**
 * @brief Main entry point for the DataGeneration program.
 * @details This function generates synthetic data using the ExaGeoStat library.
 * @param[in] argc The number of command line arguments.
 * @param[in] argv An array of command line argument strings.
 * @return An integer indicating the success or failure of the program.
 */
int main(int argc, char **argv) {

    // Create a new configurations object.
    Configurations configurations;
    //  Initialize the arguments with the provided command line arguments
    configurations.InitializeArguments(argc, argv);
    cout << "** Initialise ExaGeoStat hardware ** " << endl;
    auto hardware = ExaGeoStatHardware(EXACT_DENSE, configurations.GetCoresNumber(),
                                       configurations.GetGPUsNumbers()); // Or you could use configurations.GetComputation().
    cout << "** Create ExaGeoStat data ** " << endl;
    exageostat::dataunits::ExaGeoStatData<double> data(configurations.GetProblemSize(), configurations.GetDimension(),
                                                       hardware);
    cout << "** Generate ExaGeoStat data ** " << endl;
    ExaGeoStat<double>::ExaGeoStatGenerateData(hardware, configurations, data);
    cout << "** Finalize data generation ** " << endl;

    return 0;
}