// Copyright (c) 2017-2023 King Abdullah University of Science and Technology,
// Copyright (C) 2023 by Brightskies inc,
// All rights reserved.
// ExaGeoStat is a software package, provided by King Abdullah University of Science and Technology (KAUST).

/**
 * @file ExaGeoStat.cpp
 * @brief 
 * @version 1.0.0
 * @author Sameh Abdulah
 * @date 2023-05-30
**/
#include <api/ExaGeoStat.hpp>
#include <linear-algebra-solvers/LinearAlgebraFactory.hpp>
#include <linear-algebra-solvers/LinearAlgebraMethods.hpp>

using namespace exageostat::api;
using namespace exageostat::configurations;
using namespace exageostat::linearAlgebra;

template<typename T> void ExaGeoStat<T>::ExaGeoStatInitializeHardware(Configurations *apConfigurations) {

    auto linearAlgebraSolver = LinearAlgebraFactory<T>::CreateLinearAlgebraSolver(apConfigurations->GetComputation());
    //// TODO: Is there a better way to avoid this?
    linearAlgebraSolver->SetConfigurations(apConfigurations);

    linearAlgebraSolver->ExaGeoStatInitContext(apConfigurations->GetCoresNumber(), apConfigurations->GetGPUsNumber());
    delete linearAlgebraSolver;
}

template<typename T>
void ExaGeoStat<T>::ExaGeoStatFinalizeHardware(Configurations *apConfigurations) {

    auto linearAlgebraSolver = LinearAlgebraFactory<T>::CreateLinearAlgebraSolver(apConfigurations->GetComputation());
    linearAlgebraSolver->SetConfigurations(apConfigurations);

    linearAlgebraSolver->ExaGeoStatFinalizeContext();
    delete linearAlgebraSolver;
}