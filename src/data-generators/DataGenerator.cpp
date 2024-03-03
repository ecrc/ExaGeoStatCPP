
// Copyright (c) 2017-2023 King Abdullah University of Science and Technology,
// All rights reserved.
// ExaGeoStat is a software package, provided by King Abdullah University of Science and Technology (KAUST).

/**
 * @file DataGenerator.cpp
 * @brief Implementation of DataGenerator class
 * @version 1.0.0
 * @author Mahmoud ElKarargy
 * @date 2023-02-14
**/

#include <data-generators/DataGenerator.hpp>
#include <data-generators/concrete/SyntheticGenerator.hpp>
#include <data-loader/concrete/CSVLoader.hpp>
#include <results/Results.hpp>

using namespace exageostat::generators;
using namespace exageostat::dataLoader::csv;
using namespace exageostat::generators::synthetic;
using namespace exageostat::dataunits;
using namespace exageostat::configurations;
using namespace exageostat::common;

template<typename T>
std::unique_ptr<DataGenerator<T>> DataGenerator<T>::CreateGenerator(Configurations &apConfigurations) {

    // Check the used Data generation method, whether it's synthetic or real.
    if(apConfigurations.GetIsSynthetic() && apConfigurations.GetIsCSV()){
        throw std::domain_error("Please activate either the synthetic or the CSV file for data generation, but not both.");
    }
    if(!apConfigurations.GetIsSynthetic() && !apConfigurations.GetIsCSV()){
        throw std::domain_error("Please activate either the synthetic or the CSV file for data generation");
    }
    if(apConfigurations.GetIsSynthetic()){
        aDataSourceType = SYNTHETIC;
    }
    else if(apConfigurations.GetIsCSV()){
        aDataSourceType = CSV_FILE;
    }
    results::Results::GetInstance()->SetIsSynthetic(apConfigurations.GetIsSynthetic());

    // Return DataGenerator unique pointer of Synthetic type
    //// TODO: In case of other file support, Then we can create another layer for the factory creation depending on the file size.
    if (aDataSourceType == SYNTHETIC) {
        return std::unique_ptr<DataGenerator<T>>(SyntheticGenerator<T>::GetInstance());
    } else if (aDataSourceType == CSV_FILE) {
        return std::unique_ptr<DataGenerator<T>>(CSVLoader<T>::GetInstance());
    } else {
        throw std::runtime_error("Data Loading for this file type is unsupported for now");
    }
}

template<typename T>
DataGenerator<T>::~DataGenerator() {
    // Return DataGenerator unique pointer of Synthetic type
    if (aDataSourceType == SYNTHETIC) {
        SyntheticGenerator<T>::GetInstance()->ReleaseInstance();
    } else if (aDataSourceType == CSV_FILE) {
        CSVLoader<T>::GetInstance()->ReleaseInstance();
    } else {
        std::cerr << "Data Loading for this file type is unsupported for now" << std::endl;
        std::exit(1);
    }
}

template<typename T> DataSourceType DataGenerator<T>::aDataSourceType = common::SYNTHETIC;
