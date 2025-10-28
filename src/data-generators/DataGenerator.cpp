
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
#if DEFAULT_RUNTIME
#include <data-generators/concrete/SyntheticGenerator.hpp>
#endif
#include <data-loader/DataLoader.hpp>

using namespace exageostat::common;
using namespace exageostat::generators;
using namespace exageostat::dataLoader;
#if DEFAULT_RUNTIME
using namespace exageostat::generators::synthetic;
#endif

template<typename T>
std::unique_ptr<DataGenerator<T>> DataGenerator<T>::CreateGenerator(configurations::Configurations &apConfigurations) {
#if DEFAULT_RUNTIME
    if (apConfigurations.GetIsSynthetic()){
        aIsSynthetic = true;
        return std::unique_ptr<DataGenerator<T>>(SyntheticGenerator<T>::GetInstance());
    } else {
        aIsSynthetic = false;
        return DataLoader<T>::CreateDataLoader(apConfigurations);
    }
#else
    aIsSynthetic = false;
    return DataLoader<T>::CreateDataLoader(apConfigurations);
#endif
}

template<typename T>
DataGenerator<T>::~DataGenerator() {
#if DEFAULT_RUNTIME
    if (aIsSynthetic) {
        SyntheticGenerator<T>::GetInstance()->ReleaseInstance();
    } else {
        DataLoader<T>::ReleaseDataLoader();
    }
#else
    // For PaRSEC runtime, clean up DataLoader (ParsecLoader)
    if (!aIsSynthetic) {
        DataLoader<T>::ReleaseDataLoader();
    }
#endif
}

template<typename T> bool DataGenerator<T>::aIsSynthetic = true;
