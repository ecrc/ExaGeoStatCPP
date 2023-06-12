// Copyright (c) 2017-2023 King Abdullah University of Science and Technology,
// Copyright (C) 2023 by Brightskies inc,
// All rights reserved.
// ExaGeoStat is a software package, provided by King Abdullah University of Science and Technology (KAUST).

/**
 * @file DataGeneration.cpp
 * @brief 
 * @version 1.0.0
 * @author Sameh Abdulah
 * @date 2023-05-30
**/
#include <configurations/data-generation/concrete/SyntheticDataConfigurations.hpp>
#include <api/ExaGeoStat.hpp>

using namespace exageostat::configurations::data_configurations;
using namespace exageostat::api;

int main(int argc, char **argv) {

    // Create a new synthetic_data_configurations object with the provided command line arguments
    SyntheticDataConfigurations synthetic_data_configurations;
    synthetic_data_configurations.InitializeArguments(argc, argv);

    std::cout << "** Initialise ExaGeoStat hardware ** " << std::endl;
    // Initialise ExaGeoStat hardware with the selected number of cores and  gpus.
    ExaGeoStat<double>::ExaGeoStatInitializeHardware(&synthetic_data_configurations);

    std::cout << "** Generate ExaGeoStat data ** " << std::endl;
    ExaGeoStat<double>::ExaGeoStatGenerateData(&synthetic_data_configurations);

    //// TODO: Avoid sending configurations to finalize
    std::cout << "** Finalize ExaGeoStat hardware ** " << std::endl;
    // Finalise ExaGeoStat context.
    ExaGeoStat<double>::ExaGeoStatFinalizeHardware(&synthetic_data_configurations);

}