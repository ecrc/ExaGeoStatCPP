
// Copyright (c) 2017-2023 King Abdullah University of Science and Technology,
// All rights reserved.
// ExaGeoStat is a software package, provided by King Abdullah University of Science and Technology (KAUST).

/**
 * @file CSVLoader.cpp
 * @brief Example of the CSVLoader class
 * @version 1.0.1
 * @author Mahmoud ElKarargy
 * @date 2024-02-18
**/

#include <data-loader/concrete/CSVLoader.hpp>
#include <data-generators/DataGenerator.hpp>

using namespace std;

using namespace exageostat::kernels;
using namespace exageostat::common;
using namespace exageostat::generators;
using namespace exageostat::dataLoader::csv;

int main(int argc, char **argv) {
    LOGGER("** Example of CSV Loader **")

    // Create and Initialize a new configurations object.
    Configurations configurations;
    configurations.InitializeArguments(argc, argv);

    // Generate Data and Log it into file
    configurations.SetLogger(true);
    configurations.InitializeDataGenerationArguments();

    // Initialize ExaGeoStat Hardware and Kernel.
    auto hardware = ExaGeoStatHardware(configurations.GetComputation(), configurations.GetCoresNumber(),
                                       configurations.GetGPUsNumbers());

    Kernel<double> *kernel = exageostat::plugins::PluginRegistry<Kernel<double>>::Create(
            configurations.GetKernelName(), configurations.GetTimeSlot());
    int kernel_variables = kernel->GetVariablesNumber();

    // Create a unique pointer to a DataGenerator object and generate data
    configurations.SetIsSynthetic(true);
    unique_ptr<DataGenerator<double>> synthetic_generator = DataGenerator<double>::CreateGenerator(configurations);
    auto data = synthetic_generator->CreateData(configurations, hardware, *kernel);

    // Read csv file using data loader
    vector<double> measurements_vector;
    vector<double> x_locations;
    vector<double> y_locations;
    vector<double> z_locations;

    auto loader = CSVLoader<double>::GetInstance();
    loader->ReadData(configurations, measurements_vector, x_locations, y_locations, z_locations, kernel_variables);

    // Print loaded data
    LOGGER("Data Loaded:")
    for (int i=0;i<x_locations.size();i++) {
        LOGGER_PRECISION("X: " << x_locations[i] << " Y: " << y_locations[i], 18)
        if (configurations.GetDimension() != Dimension2D) {
            LOGGER_PRECISION(" Z: " << z_locations[i], 18)
        }
        LOGGER_PRECISION(" z: " << measurements_vector[i] << "\n", 18)
    }

    delete kernel;
}
