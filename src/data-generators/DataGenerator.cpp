
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
#include <data-generators/concrete/CSVDataGenerator.hpp>
#include <results/Results.hpp>

using namespace exageostat::generators;
using namespace exageostat::generators::synthetic;
using namespace exageostat::generators::csv;
using namespace exageostat::dataunits;
using namespace exageostat::configurations;

template<typename T>
std::unique_ptr<DataGenerator<T>> DataGenerator<T>::CreateGenerator(Configurations &apConfigurations) {

    // Check the used Data generation method, whether it's synthetic or real.
    mIsSynthetic = apConfigurations.GetIsSynthetic();
    mIsCSV = apConfigurations.GetIsCSV();
    results::Results::GetInstance()->SetIsSynthetic(mIsSynthetic);

    // Return DataGenerator unique pointer of Synthetic type
    if (mIsSynthetic) {
        return std::unique_ptr<DataGenerator<T>>(SyntheticGenerator<T>::GetInstance());
    } else if (mIsCSV) {
        return std::unique_ptr<DataGenerator<T>>(CSVDataGenerator<T>::GetInstance());
    } else {
        throw std::runtime_error("Data Loading for this file type is unsupported for now");
    }
}

template<typename T>
DataGenerator<T>::~DataGenerator() {
    // Return DataGenerator unique pointer of Synthetic type
    if (mIsSynthetic) {
        SyntheticGenerator<T>::GetInstance()->ReleaseInstance();
    } else if (mIsCSV) {
        CSVDataGenerator<T>::GetInstance()->ReleaseInstance();
    } else {
        std::cerr << "Data Loading for this file type is unsupported for now" << std::endl;
        std::exit(1);
    }
}

template<typename T> bool DataGenerator<T>::mIsSynthetic = true;
template<typename T> bool DataGenerator<T>::mIsCSV = false;
