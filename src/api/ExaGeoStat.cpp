
// Copyright (c) 2017-2023 King Abdullah University of Science and Technology,
// All rights reserved.
// ExaGeoStat is a software package, provided by King Abdullah University of Science and Technology (KAUST).

/**
 * @file ExaGeoStat.cpp
 * @brief High-Level Wrapper class containing the static API for ExaGeoStat operations.
 * @version 1.0.0
 * @author Mahmoud ElKarargy
 * @date 2023-05-30
**/

#include <api/ExaGeoStat.hpp>
#include <data-generators/DataGenerator.hpp>
#include <data-units/ModelingDataHolders.hpp>
#include <prediction/Prediction.hpp>

using namespace std;
using namespace nlopt;

using namespace exageostat::api;
using namespace exageostat::configurations;
using namespace exageostat::generators;
using namespace exageostat::dataunits;
using namespace exageostat::hardware;
using namespace exageostat::prediction;


template<typename T>
void ExaGeoStat<T>::ExaGeoStatGenerateData(const ExaGeoStatHardware &aHardware, Configurations &aConfigurations,
                                           ExaGeoStatData<T> &aData) {
    // Register and create a kernel object
    kernels::Kernel<T> *pKernel = plugins::PluginRegistry<kernels::Kernel<T>>::Create(aConfigurations.GetKernelName());
    // Add the data generation arguments.
    aConfigurations.InitializeDataGenerationArguments();
    // Create a unique pointer to a DataGenerator object
    unique_ptr<DataGenerator<T>> data_generator = DataGenerator<T>::CreateGenerator(aConfigurations);
    aData.SetLocations(*data_generator->CreateLocationsData(aConfigurations));
    // We Generate date with only Exact computation. This is a pre-request.
    auto linear_algebra_solver = linearAlgebra::LinearAlgebraFactory<T>::CreateLinearAlgebraSolver(common::EXACT_DENSE);
    linear_algebra_solver->GenerateSyntheticData(aConfigurations, aHardware, aData, *pKernel);
    delete pKernel;
}

template<typename T>
T ExaGeoStat<T>::ExaGeoStatDataModeling(const ExaGeoStatHardware &aHardware, Configurations &aConfigurations,
                                        ExaGeoStatData<T> &aData, T *apMeasurementsMatrix) {

    // Register and create a kernel object
    kernels::Kernel<T> *pKernel = plugins::PluginRegistry<kernels::Kernel<T>>::Create(aConfigurations.GetKernelName());
    // Add the data modeling arguments.
    aConfigurations.InitializeDataModelingArguments();

    int parameters_number = pKernel->GetParametersNumbers();
    int max_number_of_iterations = aConfigurations.GetMaxMleIterations();
    // Setting struct of data to pass to the modeling.
    auto modeling_data = new mModelingData(aData, aConfigurations, aHardware, *apMeasurementsMatrix, *pKernel);
    // Create nlopt
    double opt_f;
    opt optimizing_function(nlopt::LN_BOBYQA, parameters_number);
    // Initialize problem's bound.
    optimizing_function.set_lower_bounds(aConfigurations.GetLowerBounds());
    optimizing_function.set_upper_bounds(aConfigurations.GetUpperBounds());
    optimizing_function.set_ftol_abs(pow(10, -1 * aConfigurations.GetTolerance()));
    // Set max iterations value.
    optimizing_function.set_maxeval(max_number_of_iterations);
    optimizing_function.set_max_objective(ExaGeoStatMLETileAPI, (void *) modeling_data);
    // Optimize mle using nlopt.
    optimizing_function.optimize(aConfigurations.GetStartingTheta(), opt_f);

    delete pKernel;
    delete modeling_data;
    return optimizing_function.last_optimum_value();
}

template<typename T>
double
ExaGeoStat<T>::ExaGeoStatMLETileAPI(const std::vector<double> &aTheta, std::vector<double> &aGrad, void *apInfo) {
    auto config = ((mModelingData<T> *) apInfo)->mpConfiguration;
    auto data = ((mModelingData<T> *) apInfo)->mpData;
    auto hardware = ((mModelingData<T> *) apInfo)->mpHardware;
    auto measurements = ((mModelingData<T> *) apInfo)->mpMeasurementsMatrix;
    auto kernel = ((mModelingData<T> *) apInfo)->mpKernel;

    // We do Date Modeling with any computation.
    auto linear_algebra_solver = linearAlgebra::LinearAlgebraFactory<T>::CreateLinearAlgebraSolver(
            config->GetComputation());
    return linear_algebra_solver->ExaGeoStatMLETile(*hardware, *data, *config, aTheta.data(), measurements, *kernel);
}


template<typename T>
void ExaGeoStat<T>::ExaGeoStatPrediction(const ExaGeoStatHardware &aHardware, Configurations &aConfigurations,
                                         ExaGeoStatData<T> &aData, T *apMeasurementsMatrix) {

    Prediction<T> predictor;
    // Register and create a kernel object
    kernels::Kernel<T> *pKernel = plugins::PluginRegistry<kernels::Kernel<T>>::Create(aConfigurations.GetKernelName());
    // Add the data prediction arguments.
    aConfigurations.InitializeDataPredictionArguments();
    predictor.PredictMissingData(aHardware, aData, aConfigurations, apMeasurementsMatrix, *pKernel);
    delete pKernel;
}