
// Copyright (c) 2017-2023 King Abdullah University of Science and Technology,
// All rights reserved.
// ExaGeoStat is a software package, provided by King Abdullah University of Science and Technology (KAUST).

/**
 * @file DataModeling.cpp
 * @brief This program models data using the ExaGeoStat library.
 * @details The program takes command line arguments and example variables to configure the data modeling.
 * @version 1.0.1
 * @author Mahmoud ElKarargy
 * @date 2023-06-21
**/

#include <iostream>

#include <common/Utils.hpp>
#include <api/ExaGeoStat.hpp>

using namespace exageostat::api;
using namespace exageostat::dataunits;
using namespace exageostat::configurations;
using namespace exageostat::hardware;

/**
 * @brief Main entry point for the Data Modeling program.
 * @details This example illustrates the process of data modeling using the ExaGeoStat library's CHAMELEON descriptor framework.
 * It involves configuring parameters such as problem size and computation mode, initializing hardware resources, setting up matrices for descriptors,
 * and creating location information. The ExaGeoStatDataModeling function is then called to perform geo statistical analysis. The example showcases
 * the library's efficiency in handling large spatial datasets while efficiently utilizing hardware resources..
 * @param[in] argc The number of command line arguments.
 * @param[in] argv An array of command line argument strings.
 * @return An integer indicating the success or failure of the program.
 */
int main(int argc, char **argv) {
    // Create a new data_modeling_configurations object with the provided command line arguments and example variables
    Configurations configurations;
    configurations.InitializeArguments(argc, argv);
    /**
     * Since this example is currently relying on user inputs instead of reading files, The following points are important to know:
     * The N and dts have to match with the Location X, Y and Z_values you're going to provide.
     * You have to provide Locations and Z values in order to use Modeling without generation.
     */
    int N = 16;
    int dts = 8;
    configurations.SetProblemSize(N);
    configurations.SetDenseTileSize(dts);

    // initialize ExaGeoStat hardware with the selected number of cores and  gpus.
    auto hardware = ExaGeoStatHardware(configurations.GetComputation(), configurations.GetCoresNumber(),
                                       configurations.GetGPUsNumbers());

    //Data Setup
    std::unique_ptr<ExaGeoStatData<double>> data = std::make_unique<ExaGeoStatData<double>>(configurations.GetProblemSize(), configurations.GetDimension());

    // Initiating the matrix of the CHAMELEON Descriptor Z.
    auto *z_matrix = new double[N]{-1.272336140360187606, -2.590699695867695773, 0.512142584178685967,
                                   -0.163880452049749520, 0.313503633252489700, -1.474410682226017677,
                                   0.161705025505231914, 0.623389205185149065, -1.341858445399783495,
                                   -1.054282062428600009, -1.669383221392507943, 0.219170645803740793,
                                   0.971213790000161170, 0.538973474182433021, -0.752828466476077041,
                                   0.290822066007430102};
    //creating locations x and y.
    auto *location_x = new double[N]{0.193041886015106440, 0.330556191348134576, 0.181612878614480805,
                                     0.370473792629892440, 0.652140077821011688, 0.806332494087129037,
                                     0.553322652018005678, 0.800961318379491916, 0.207324330510414295,
                                     0.347951476310368490, 0.092042420080872822, 0.465445944914930965,
                                     0.528267338063630132, 0.974792095826657490, 0.552452887769893985,
                                     0.877592126344701295};

    auto *location_y = new double[N]{0.103883421072709245, 0.135790035858701447, 0.434683756771190977,
                                     0.400778210116731537, 0.168459601739528508, 0.105195696955825133,
                                     0.396398870832379624, 0.296757457846952011, 0.564507515068284116,
                                     0.627679865720607300, 0.928648813611047563, 0.958236057068741931,
                                     0.573571374074921758, 0.568657969024185528, 0.935835812924391552,
                                     0.942824444953078489};

    data->GetLocations()->SetLocationX(*location_x, N);
    data->GetLocations()->SetLocationY(*location_y, N);

    // Modeling module.
    ExaGeoStat<double>::ExaGeoStatDataModeling(hardware, configurations, data, z_matrix);

    // Freeing the allocated memory.
    delete[] z_matrix;
    delete[] location_x;
    delete[] location_y;

    return 0;
}
