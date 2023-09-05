
// Copyright (c) 2017-2023 King Abdullah University of Science and Technology,
// All rights reserved.
// ExaGeoStat is a software package, provided by King Abdullah University of Science and Technology (KAUST).

/**
 * @file DataGenerationAndModeling.cpp
 * @brief This program generates synthetic data then performs data modeling on generated data using the ExaGeoStat library.
 * @details The program takes command line arguments to configure the data generation.
 * @version 1.0.0
 * @author Mahmoud ElKarargy
 * @date 2023-06-21
**/

#include <configurations/Configurations.hpp>
#include <api/ExaGeoStat.hpp>

using namespace std;

using namespace exageostat::configurations;
using namespace exageostat::api;
using namespace exageostat::common;
using namespace exageostat::hardware;
using namespace exageostat::dataunits;

/**
 * @brief Main entry point for the Data Generation & Data Modeling program.
 * @details This function generates synthetic data using the ExaGeoStat library and models it.
 * @param[in] argc The number of command line arguments.
 * @param[in] argv An array of command line argument strings.
 * @return An integer indicating the success or failure of the program.
 */
int main(int argc, char **argv) {

    // Create a new configurations object. it needs to be a heap variable
    Configurations configurations;
    //  Initialize the arguments with the provided command line arguments
    configurations.InitializeArguments(argc, argv);
    cout << "** initialize ExaGeoStat hardware ** " << endl;
    auto hardware = ExaGeoStatHardware(configurations.GetComputation(), configurations.GetCoresNumber(),
                                       configurations.GetGPUsNumbers()); // Or you could use configurations.GetComputation().
    cout << "** Create ExaGeoStat data ** " << endl;
    ExaGeoStatData<double> data(configurations.GetProblemSize(), configurations.GetDimension(), hardware);
    cout << "** ExaGeoStat data generation** " << endl;
    ExaGeoStat<double>::ExaGeoStatGenerateData(hardware, configurations, data);
    cout << "** ExaGeoStat data Modeling** " << endl;
    ExaGeoStat<double>::ExaGeoStatDataModeling(hardware, configurations, data);
    cout << "** Finalize data generation and modeling ** " << endl;

    return 0;
}

