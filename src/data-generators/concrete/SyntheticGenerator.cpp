
// Copyright (c) 2017-2024 King Abdullah University of Science and Technology,
// All rights reserved.
// ExaGeoStat is a software package, provided by King Abdullah University of Science and Technology (KAUST).

/**
 * @file SyntheticGenerator.cpp
 * @brief Implementation of the SyntheticGenerator class
 * @version 1.1.0
 * @author Mahmoud ElKarargy
 * @author Sameh Abdulah
 * @date 2025-02-18
**/

#include <data-generators/concrete/SyntheticGenerator.hpp>
#if !DEFAULT_RUNTIME
#include <data-generators/concrete/runtime/ParsecGenerator.hpp>
#else
#include <data-generators/concrete/runtime/StarpuGenerator.hpp>
#endif

using namespace exageostat::configurations;
using namespace exageostat::generators::synthetic;

template<typename T>
std::unique_ptr<ExaGeoStatData<T>>
SyntheticGenerator<T>::CreateData(Configurations &aConfigurations, exageostat::kernels::Kernel<T> &aKernel) {
    
    auto data = this->CreateSyntheticData(aConfigurations, aKernel);
    return data;
}

template<typename T>
std::unique_ptr<SyntheticGenerator<T>> SyntheticGenerator<T>::CreateSyntheticGenerator() {

#if DEFAULT_RUNTIME
    return std::unique_ptr<SyntheticGenerator<T>>(starpu::StarpuGenerator<T>::GetInstance());
#else
    return std::unique_ptr<SyntheticGenerator<T>>(parsec::ParsecGenerator<T>::GetInstance());
#endif
}

template<typename T>
void SyntheticGenerator<T>::ReleaseSyntheticGenerator() {

#if DEFAULT_RUNTIME
    starpu::StarpuGenerator<T>::GetInstance()->ReleaseInstance();
#else
    parsec::ParsecGenerator<T>::GetInstance()->ReleaseInstance();
#endif
}
