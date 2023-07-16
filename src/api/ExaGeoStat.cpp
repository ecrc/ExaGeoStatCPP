// Copyright (c) 2017-2023 King Abdullah University of Science and Technology,
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
#include <data-generators/DataGenerator.hpp>

#include <nlopt.hpp>

using namespace std;

using namespace exageostat::api;
using namespace exageostat::configurations;
using namespace exageostat::linearAlgebra;
using namespace exageostat::generators;
using namespace exageostat::dataunits;

template<typename T>
void ExaGeoStat<T>::ExaGeoStatInitializeHardware() {

    auto linearAlgebraSolver = LinearAlgebraFactory<T>::CreateLinearAlgebraSolver(
            Configurations::GetConfigurations()->GetComputation());
    linearAlgebraSolver->ExaGeoStatInitContext();
    delete linearAlgebraSolver;
}

template<typename T>
void ExaGeoStat<T>::ExaGeoStatFinalizeHardware() {

    auto linearAlgebraSolver = LinearAlgebraFactory<T>::CreateLinearAlgebraSolver(
            Configurations::GetConfigurations()->GetComputation());
    linearAlgebraSolver->DestroyDescriptors();
    linearAlgebraSolver->ExaGeoStatFinalizeContext();
    delete linearAlgebraSolver;
}

template<typename T>
Locations *ExaGeoStat<T>::ExaGeoStatGenerateData() {
    // Create a unique pointer to a DataGenerator object
    unique_ptr<DataGenerator<double>> data_generator;
    // Create the DataGenerator object
    data_generator = data_generator->CreateGenerator();

    data_generator->GenerateLocations();
    data_generator->GenerateDescriptors();
    data_generator->GenerateObservations();
    return data_generator->GetLocations();
}

static double ExaGeoStatMleTile(unsigned n, const double* theta, double* grad, void* apDataLocations) {

}

template<typename T>
void ExaGeoStat<T>::ExaGeoStatDataModeling(dataunits::Locations *apDataLocations) {

    auto linear_algebra_solver = LinearAlgebraFactory<T>::CreateLinearAlgebraSolver(
            Configurations::GetConfigurations()->GetComputation());
    linear_algebra_solver->InitiateDescriptors();

    nlopt_opt opt;
    double opt_f;

    nlopt_set_max_objective(opt, ExaGeoStatMleTile, apDataLocations);
    nlopt_optimize(opt, (double *) Configurations::GetConfigurations()->GetStartingTheta().data(), &opt_f);
//    linear_algebra_solver->ExaGeoStatMleTile(apDataLocations);
}
