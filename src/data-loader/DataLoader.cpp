
// Copyright (c) 2017-2024 King Abdullah University of Science and Technology,
// All rights reserved.
// ExaGeoStat is a software package, provided by King Abdullah University of Science and Technology (KAUST).

/**
 * @file SyntheticGenerator.cpp
 * @brief Implementation of the SyntheticGenerator class
 * @version 1.1.0
 * @author Mahmoud ElKarargy
 * @author Sameh Abdulah
 * @date 2023-02-14
**/

#if !DEFAULT_RUNTIME
#include <data-loader/concrete/ParsecLoader.hpp>
#else
#include <data-loader/concrete/CSVLoader.hpp>
#endif

using namespace std;

using namespace exageostat::dataLoader;

template<typename T>
std::unique_ptr<ExaGeoStatData<T>>
DataLoader<T>::CreateData(configurations::Configurations &aConfigurations, exageostat::kernels::Kernel<T> &aKernel) {

    auto data = this->LoadData(aConfigurations, aKernel);
    return data;
}

template<typename T>
std::unique_ptr<DataLoader<T>>
DataLoader<T>::CreateDataLoader(exageostat::configurations::Configurations &apConfigurations){

#if DEFAULT_RUNTIME
    return std::unique_ptr<DataLoader<T>>(csv::CSVLoader<T>::GetInstance());
#else
    return std::unique_ptr<DataLoader<T>>(parsec::ParsecLoader<T>::GetInstance());
#endif
}

template<typename T>
void DataLoader<T>::ReleaseDataLoader() {

#if DEFAULT_RUNTIME
    csv::CSVLoader<T>::GetInstance()->ReleaseInstance();
#else
    parsec::ParsecLoader<T>::GetInstance()->ReleaseInstance();
#endif
}