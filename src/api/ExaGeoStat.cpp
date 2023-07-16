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

using namespace std;

using namespace exageostat::api;
using namespace exageostat::configurations;
using namespace exageostat::linearAlgebra;
using namespace exageostat::generators;
using namespace exageostat::dataunits;
using namespace exageostat::configurations::data_configurations;
using namespace exageostat::configurations::data_modeling;

template<typename T>
void ExaGeoStat<T>::ExaGeoStatInitializeHardware(Configurations *apConfigurations) {

    auto linearAlgebraSolver = LinearAlgebraFactory<T>::CreateLinearAlgebraSolver(apConfigurations->GetComputation());
    linearAlgebraSolver->SetConfigurations(apConfigurations);

    linearAlgebraSolver->ExaGeoStatInitContext(apConfigurations->GetCoresNumber(), apConfigurations->GetGPUsNumber());
    delete linearAlgebraSolver;
}

template<typename T>
void ExaGeoStat<T>::ExaGeoStatFinalizeHardware(Configurations *apConfigurations) {

    auto linearAlgebraSolver = LinearAlgebraFactory<T>::CreateLinearAlgebraSolver(apConfigurations->GetComputation());
    linearAlgebraSolver->SetConfigurations(apConfigurations);

    linearAlgebraSolver->DestoryDescriptors();
    linearAlgebraSolver->ExaGeoStatFinalizeContext();
    delete linearAlgebraSolver;
}

template<typename T>
Locations* ExaGeoStat<T>::ExaGeoStatGenerateData(SyntheticDataConfigurations *apConfigurations) {
    // Create a unique pointer to a DataGenerator object
    unique_ptr<DataGenerator<double>> data_generator;
    // Create the DataGenerator object
    data_generator = data_generator->CreateGenerator(
            (data_configurations::SyntheticDataConfigurations *) apConfigurations);

    data_generator->GenerateLocations();
    data_generator->GenerateDescriptors();
    data_generator->GenerateObservations();
    return data_generator->GetLocations();
}

template<typename T>
 void ExaGeoStat<T>::ExaGeoStatDataModeling(DataModelingConfigurations *apConfigurations, dataunits::Locations *apDataLocations){

     auto linear_algebra_solver = LinearAlgebraFactory<T>::CreateLinearAlgebraSolver(apConfigurations->GetComputation());
     linear_algebra_solver->SetConfigurations(apConfigurations);
     linear_algebra_solver->InitiateDescriptors();
     linear_algebra_solver->ExageostatMleTile(apConfigurations->GetInitialTheta(), apDataLocations, (configurations::data_modeling::DataModelingConfigurations*) apConfigurations);
 }

