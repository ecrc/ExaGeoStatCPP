
// Copyright (c) 2017-2024 King Abdullah University of Science and Technology,
// All rights reserved.
// ExaGeoStat is a software package, provided by King Abdullah University of Science and Technology (KAUST).

/**
 * @file DataGenerator.cpp
 * @brief Implementation of DataGenerator class
 * @version 1.1.0
 * @author Mahmoud ElKarargy
 * @date 2024-02-04
**/

#include <data-generators/DataGenerator.hpp>
#include <data-generators/concrete/SyntheticGenerator.hpp>
#include <data-loader/concrete/CSVLoader.hpp>

using namespace exageostat::generators;
using namespace exageostat::dataLoader::csv;
using namespace exageostat::generators::synthetic;
using namespace exageostat::common;
using namespace exageostat::results;

template<typename T>
std::unique_ptr<DataGenerator<T>> DataGenerator<T>::CreateGenerator(configurations::Configurations &apConfigurations) {

    // Check the used Data generation method, whether it's synthetic or real.
#if DEFAULT_RUNTIME
    aDataSourceType = apConfigurations.GetIsSynthetic() ? SYNTHETIC : CSV_FILE;
#else
    aDataSourceType = CSV_FILE;
#endif
    // Return DataGenerator unique pointer of Synthetic type
    if (aDataSourceType == SYNTHETIC) {
        Results::GetInstance()->SetIsSynthetic(true);
        return std::unique_ptr<DataGenerator<T>>(SyntheticGenerator<T>::GetInstance());
    } else if (aDataSourceType == CSV_FILE) {
        Results::GetInstance()->SetIsSynthetic(false);
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
        exit(1);
    }
}

template<typename T> DataSourceType DataGenerator<T>::aDataSourceType = SYNTHETIC;
