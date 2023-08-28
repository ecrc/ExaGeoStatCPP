
// Copyright (c) 2017-2023 King Abdullah University of Science and Technology,
// All rights reserved.
// ExaGeoStat is a software package, provided by King Abdullah University of Science and Technology (KAUST).

/**
 * @file ExaGeoStat.cpp
 * @brief High-Level Wrapper class containing the static API for ExaGeoStat operations.
 * @version 1.0.0
 * @author Sameh Abdulah
 * @author Mahmoud ElKarargy
 * @date 2023-05-30
**/

#include <api/ExaGeoStat.hpp>
#include <configurations/Configurations.hpp>
#include <data-generators/DataGenerator.hpp>
#include <linear-algebra-solvers/LinearAlgebraFactory.hpp>
#include <data-units/ModelingDataHolders.hpp>

using namespace std;
using namespace nlopt;

using namespace exageostat::api;
using namespace exageostat::configurations;
using namespace exageostat::linearAlgebra;
using namespace exageostat::generators;
using namespace exageostat::dataunits;
using namespace exageostat::hardware;


template<typename T>
void ExaGeoStat<T>::ExaGeoStatGenerateData(const ExaGeoStatHardware &aHardware, Configurations &aConfigurations,
                                           ExaGeoStatData<T> &aData) {

    // Add the data generation arguments.
    aConfigurations.InitializeDataGenerationArguments();
    // Create a unique pointer to a DataGenerator object
    unique_ptr<DataGenerator<T>> data_generator = DataGenerator<T>::CreateGenerator(aConfigurations);

    aData.SetLocations(*data_generator->CreateLocationsData(aConfigurations));

    auto linear_algebra_solver = LinearAlgebraFactory<T>::CreateLinearAlgebraSolver(aConfigurations.GetComputation());
#ifdef EXAGEOSTAT_USE_CHAMELEON
    linear_algebra_solver->GenerateSyntheticData(aConfigurations, aHardware, aData, common::CHAMELEON_DESCRIPTOR);
#endif
#ifdef EXAGEOSTAT_USE_HICMA
    linear_algebra_solver->GenerateSyntheticData(aConfigurations, aHardware, aData, common::HICMA_DESCRIPTOR);
#endif
}

template<typename T>
void ExaGeoStat<T>::ExaGeoStatDataModeling(const ExaGeoStatHardware &aHardware, Configurations &aConfigurations,
                                           ExaGeoStatData<T> &aData, T *apMeasurementsMatrix) {

    // Add the data modeling arguments.
    aConfigurations.InitializeDataModelingArguments();
    int max_number_of_iterations = aConfigurations.GetMaxMleIterations();

    // Setting struct of data to pass to the modeling.
    auto modeling_data = new mModelingData(&aData, &aConfigurations, &aHardware, apMeasurementsMatrix);

    // Create a kernel object depending on which kernel the user is going to use.
    kernels::Kernel<T> *kernel = exageostat::plugins::PluginRegistry<kernels::Kernel<T>>::Create(
            aConfigurations.GetKernelName());

    // Create nlopt
    double opt_f;
    opt optimizing_function(nlopt::LN_BOBYQA, kernel->GetParametersNumbers());
    delete kernel;

    // Initialize problem's bound.
    optimizing_function.set_lower_bounds(aConfigurations.GetLowerBounds());
    optimizing_function.set_upper_bounds(aConfigurations.GetUpperBounds());
    optimizing_function.set_ftol_abs(pow(10, -1 * aConfigurations.GetTolerance()));

    // Set max iterations value.
    optimizing_function.set_maxeval(max_number_of_iterations);

    optimizing_function.set_max_objective(ExaGeoStatMleTileAPI, (void *) modeling_data);
    // Optimize mle using nlopt.
    optimizing_function.optimize(aConfigurations.GetStartingTheta(), opt_f);

    delete modeling_data;
}

template<typename T>
double
ExaGeoStat<T>::ExaGeoStatMleTileAPI(const std::vector<double> &aTheta, std::vector<double> &aGrad, void *apInfo) {
    auto config = ((mModelingData<T> *) apInfo)->mpConfiguration;
    auto data = ((mModelingData<T> *) apInfo)->mpData;
    auto hardware = ((mModelingData<T> *) apInfo)->mpHardware;
    auto measurements = ((mModelingData<T> *) apInfo)->mpMeasurementsMatrix;

    auto linear_algebra_solver = linearAlgebra::LinearAlgebraFactory<T>::CreateLinearAlgebraSolver(
            config->GetComputation());

    return linear_algebra_solver->ExaGeoStatMleTile(*hardware, *data, *config, aTheta.data(), measurements);
}
