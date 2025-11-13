
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
#else
#include <data-generators/concrete/runtime/ParsecGenerator.hpp>
#endif
#include <data-loader/DataLoader.hpp>

using namespace exageostat::common;
using namespace exageostat::generators;
using namespace exageostat::dataLoader;
#if DEFAULT_RUNTIME
using namespace exageostat::generators::synthetic;
#else
using namespace exageostat::generators::synthetic::parsec;
#endif

template<typename T>
std::unique_ptr<DataGenerator<T>> DataGenerator<T>::CreateGenerator(configurations::Configurations &apConfigurations) {
#if DEFAULT_RUNTIME
    if (apConfigurations.GetIsSynthetic()){
        aIsSynthetic = true;
        return SyntheticGenerator<T>::CreateSyntheticGenerator();
    } else {
        aIsSynthetic = false;
        return DataLoader<T>::CreateDataLoader(apConfigurations);
    }
#else
    // PaRSEC runtime - check if climate emulator or general operations
    if (apConfigurations.GetIsClimateEmulator()) {
        aIsSynthetic = false;
        return DataLoader<T>::CreateDataLoader(apConfigurations);
    } else {
        aIsSynthetic = true;
        return std::unique_ptr<DataGenerator<T>>(ParsecGenerator<T>::GetInstance());
    }
#endif
}

template<typename T>
DataGenerator<T>::~DataGenerator() {
#if DEFAULT_RUNTIME
    if (aIsSynthetic) {
        SyntheticGenerator<T>::ReleaseSyntheticGenerator();
    } else {
        DataLoader<T>::ReleaseDataLoader();
    }
#else
    // PaRSEC runtime cleanup
    if (aIsSynthetic) {
        ParsecGenerator<T>::GetInstance()->ReleaseInstance();
    } else {
        DataLoader<T>::ReleaseDataLoader();
    }
#endif
}

template<typename T> bool DataGenerator<T>::aIsSynthetic = true;
