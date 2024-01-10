
// Copyright (c) 2017-2023 King Abdullah University of Science and Technology,
// All rights reserved.
// ExaGeoStat is a software package, provided by King Abdullah University of Science and Technology (KAUST).

/**
 * @file SyntheticLocationsGeneration.cpp
 * @brief This file contains the main function for generating synthetic Locations for ExaGeoStat
 * @version 1.0.1
 * @author Mahmoud ElKarargy
 * @date 2023-03-04
**/

#include <iostream>

#include <data-generators/DataGenerator.hpp>

using namespace std;

using namespace exageostat::configurations;
using namespace exageostat::generators;
using namespace exageostat::common;
using namespace exageostat::hardware;
using namespace exageostat::kernels;

/**
 * @brief The main function of the program.
 * @details This function generates synthetic data for ExaGeoStat using the provided command line arguments.
 * @param[in] argc The number of command line arguments.
 * @param[in] argv The command line arguments.
 * @return The status code of the program.
 */

int main(int argc, char **argv) {

    // Create a new synthetic_data_configurations object with the provided command line arguments
    Configurations synthetic_data_configurations;
    synthetic_data_configurations.InitializeArguments(argc, argv);
    synthetic_data_configurations.InitializeDataGenerationArguments();

    // initialize ExaGeoStat Hardware.
    auto hardware = ExaGeoStatHardware(synthetic_data_configurations.GetComputation(),
                                       synthetic_data_configurations.GetCoresNumber(),
                                       synthetic_data_configurations.GetGPUsNumbers());

    Kernel<double> *pKernel = exageostat::plugins::PluginRegistry<Kernel<double>>::Create(
            synthetic_data_configurations.GetKernelName(), synthetic_data_configurations.GetTimeSlot());

    // Create a unique pointer to a DataGenerator object
    unique_ptr<DataGenerator<double>> synthetic_generator = DataGenerator<double>::CreateGenerator(
            synthetic_data_configurations);

    // Initialize the locations of the generated data
    auto data = synthetic_generator->CreateData(synthetic_data_configurations, hardware, *pKernel);
    // Define a struct to hold pointers to the x, y, and z coordinates of the generated data
    struct DataPointers {
        double *x;
        double *y;
        double *z;
    } data_pointers{};

    // Set the pointers in the DataPointers struct to the location coordinates of the generated data
    data_pointers.x = data->GetLocations()->GetLocationX();
    data_pointers.y = data->GetLocations()->GetLocationY();
    data_pointers.z = data->GetLocations()->GetLocationZ();

    // Print the generated location coordinates
    LOGGER("Generated Data ...")
    int timeSlot;
    if (synthetic_data_configurations.GetDimension() != DimensionST) {
        timeSlot = 1;
    } else {
        timeSlot = synthetic_data_configurations.GetTimeSlot();
    }
    for (auto i = 0; i < synthetic_data_configurations.GetProblemSize() * timeSlot; i++) {
        LOGGER_PRECISION("X: " << data_pointers.x[i] << " Y: " << data_pointers.y[i], 18)
        if (synthetic_data_configurations.GetDimension() != Dimension2D) {
            LOGGER_PRECISION(" Z: " << data_pointers.z[i], 18)
        }
        LOGGER_PRECISION(" Measurements: " << ((double *) data->GetDescriptorData()->GetDescriptor(
                exageostat::common::CHAMELEON_DESCRIPTOR, exageostat::common::DESCRIPTOR_Z).chameleon_desc->mat)[i]
                                           << "\n", 18)
    }

    delete pKernel;

    return 0;
}