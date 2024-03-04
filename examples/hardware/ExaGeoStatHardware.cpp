
// Copyright (c) 2017-2023 King Abdullah University of Science and Technology,
// All rights reserved.
// ExaGeoStat is a software package, provided by King Abdullah University of Science and Technology (KAUST).

/**
 * @file TestExaGeoStatHardware.cpp
 * @briefExample for the ExaGeoStatHardware class in the ExaGeoStat software package.
 * @version 1.1.0
 * @author Mahmoud ElKarargy
 * @date 2024-02-14
**/

#include <configurations/Configurations.hpp>
#include <hardware/ExaGeoStatHardware.hpp>
#include <utilities/Logger.hpp>

using namespace exageostat::common;

int main(int argc, char **argv) {

    LOGGER("** Example of Hardware **")

    // Initialize Configuration
    Configurations configuration;
    configuration.InitializeArguments(argc,argv);

    // Initialize Hardware
    auto hardware = ExaGeoStatHardware(configuration.GetComputation(), configuration.GetCoresNumber(), configuration.GetGPUsNumbers());
}