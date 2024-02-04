
// Copyright (c) 2017-2023 King Abdullah University of Science and Technology,
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
#include <data-generators/concrete/CSVDataGenerator.hpp>
#include <results/Results.hpp>

using namespace exageostat::generators;
using namespace exageostat::generators::synthetic;
using namespace exageostat::generators::csv;
using namespace exageostat::dataunits;

template<typename T>
std::unique_ptr<DataGenerator<T>> DataGenerator<T>::CreateGenerator(Configurations &apConfigurations) {

    // Check the used Data generation method, whether it's synthetic or real.
    mIsSynthetic = apConfigurations.GetIsSynthetic();
    results::Results::GetInstance()->SetIsSynthetic(mIsSynthetic);

    // Return DataGenerator unique pointer of Synthetic type
    if (mIsSynthetic) {
        return std::unique_ptr<DataGenerator<T>>(SyntheticGenerator<T>::GetInstance());
    } else {
        return std::unique_ptr<DataGenerator<T>>(CSVDataGenerator<T>::GetInstance());
    }
}

template<typename T>
DataGenerator<T>::~DataGenerator() {
    // Return DataGenerator unique pointer of Synthetic type
    if (mIsSynthetic) {
        SyntheticGenerator<T>::GetInstance()->ReleaseInstance();
    } else {
        CSVDataGenerator<T>::GetInstance()->ReleaseInstance();
    }
}

template<typename T> bool DataGenerator<T>::mIsSynthetic = true;
