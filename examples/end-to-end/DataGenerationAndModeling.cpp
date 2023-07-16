
/*
 * Copyright (c) 2017-2023 King Abdullah University of Science and Technology,
 * All rights reserved.
 * ExaGeoStat is a software package, provided by King Abdullah University of Science and Technology (KAUST).
 */

/**
 * @file DataGenerationAndModeling.cpp
 * @brief This program generates synthetic data then performs data modeling on generated data using the ExaGeoStat library.
 * @details The program takes command line arguments to configure the data generation.
 * @version 1.0.0
 * @author Sameh Abdulah
 * @date 2023-06-21
**/

#include <configurations/data-generation/concrete/SyntheticDataConfigurations.hpp>
#include <configurations/data-modeling/DataModelingConfigurations.hpp>
#include <api/ExaGeoStat.hpp>

using namespace exageostat::configurations::data_configurations;
using namespace exageostat::configurations::data_modeling;
using namespace exageostat::configurations;
using namespace exageostat::api;

/**
 * @brief Main entry point for the Data Generation & Data Modeling program.
 * @details This function generates synthetic data using the ExaGeoStat library and models it.
 * @param argc The number of command line arguments.
 * @param argv An array of command line argument strings.
 * @return An integer indicating the success or failure of the program.
 */
int main(int argc, char **argv) {

    // Creating Synthetic data configurations.
    Configurations *configurations = new SyntheticDataConfigurations(argc, argv);

    std::cout << "** Initialise ExaGeoStat hardware ** " << std::endl;
    ExaGeoStat<double>::ExaGeoStatInitializeHardware(configurations);

    std::cout << "** ExaGeoStat Data Generation ** " << std::endl;
    exageostat::dataunits::Locations *data = ExaGeoStat<double>::ExaGeoStatGenerateData((SyntheticDataConfigurations *) configurations);
    std::cout << "X: "<<"\t\t" << "Y" <<std::endl;
//    std::cout << "msize = "<< data->GetSize() <<std::endl;
//    for(int i = 0; i < data->GetSize()-1; i++){
//    std::cout
//    <<data->GetLocationX()[i]
//    <<  "\t\t" << data->GetLocationY()[i]
//    <<"  ";}

    // Changing ExaGeoStat Configurations to Data Modeling Configurations
//    configurations = new DataModelingConfigurations(*(DataModelingConfigurations *) configurations);
    std::cout << "** ExaGeoStat Data Modeling ** " << std::endl;
////    ExaGeoStat<double>::ExaGeoStatDataModeling(&&configurations.GetConfigurations(), ??data??);
    auto *dataModelingConfigurations = dynamic_cast<DataModelingConfigurations *>(configurations);
    ExaGeoStat<double>::ExaGeoStatDataModeling(dataModelingConfigurations, data);

    std::cout << "** Finalize ExaGeoStat hardware ** " << std::endl;
    ExaGeoStat<double>::ExaGeoStatFinalizeHardware(configurations);

    delete configurations;
    return 0;
}
