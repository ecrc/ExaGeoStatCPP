
/*
 * Copyright (c) 2017-2023 King Abdullah University of Science and Technology,
 * All rights reserved.
 * ExaGeoStat is a software package, provided by King Abdullah University of Science and Technology (KAUST).
 */

/**
 * @file DataModeling.cpp
 * @brief This program models data using the ExaGeoStat library.
 * @details The program takes command line arguments to configure the data generation.
 * @version 1.0.0
 * @author Sameh Abdulah
 * @date 2023-06-21
**/
#include <cmath>

#include <configurations/data-modeling/DataModelingConfigurations.hpp>
#include <api/ExaGeoStat.hpp>

using namespace exageostat::configurations::data_modeling;
using namespace exageostat::api;
using namespace exageostat::dataunits;
using namespace exageostat::common;

/**
 * @brief Main entry point for the Data Modeling program.
 * @details This function models data using the ExaGeoStat library.
 * @param argc The number of command line arguments.
 * @param argv An array of command line argument strings.
 * @return An integer indicating the success or failure of the program.
 */
int main(int argc, char **argv) {
    //Data Setup
    std::cout << "** ExaGeoStat Initiation ** " << std::endl;

    int problem_size = 9;
    Dimension dimension = Dimension2D;
    auto* x = new double[9]{0.257389, 0.456062, 0.797269, 0.242161, 0.440742, 0.276432, 0.493965, 0.953933, 0.86952};
    auto* y = new double[9]{0.138506, 0.238193, 0.170245, 0.579583, 0.514397, 0.752682, 0.867704, 0.610986, 0.891279};

    auto *data = new Locations(problem_size, dimension);
    data->SetLocationX(x);
    data->SetLocationY(y);

    for(int i = 0; i < problem_size; i++){
        std::cout << "[" << i << "] X: " << data->GetLocationX()[i] << " Y: " << data->GetLocationY()[i] << std::endl;
    }

    // Create a new data_modeling_configurations object with the provided command line arguments
    DataModelingConfigurations data_modeling_configurations(argc, argv);

    // Initialise ExaGeoStat hardware with the selected number of cores and  gpus.
    ExaGeoStat<double>::ExaGeoStatInitializeHardware(&data_modeling_configurations);

    std::cout << "** ExaGeoStat Data Modeling ** " << std::endl;
    ExaGeoStat<double>::ExaGeoStatDataModeling(&data_modeling_configurations, data);

    std::cout << "** Finalize ExaGeoStat hardware ** " << std::endl;
    // Finalise ExaGeoStat context.
    ExaGeoStat<double>::ExaGeoStatFinalizeHardware(&data_modeling_configurations);

    return 0;
}
