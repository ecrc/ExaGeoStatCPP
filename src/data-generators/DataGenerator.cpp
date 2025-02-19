
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
#include <data-loader/DataLoader.hpp>

using namespace exageostat::common;
using namespace exageostat::generators;
using namespace exageostat::dataLoader;
using namespace exageostat::generators::synthetic;

template<typename T>
std::unique_ptr<DataGenerator<T>> DataGenerator<T>::CreateGenerator(configurations::Configurations &apConfigurations) {
    if (apConfigurations.GetIsSynthetic()){
        aIsSynthetic = true;
        return std::unique_ptr<DataGenerator<T>>(SyntheticGenerator<T>::GetInstance());
    } else {
        aIsSynthetic = false;
        return DataLoader<T>::CreateDataLoader(apConfigurations);
    }
}

template<typename T>
DataGenerator<T>::~DataGenerator() {
    if (aIsSynthetic) {
        SyntheticGenerator<T>::GetInstance()->ReleaseInstance();
    } else {
        DataLoader<T>::ReleaseDataLoader();
    }
}

template<typename T> bool DataGenerator<T>::aIsSynthetic = true;
