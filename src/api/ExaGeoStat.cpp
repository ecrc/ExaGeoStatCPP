
// Copyright (c) 2017-2024 King Abdullah University of Science and Technology,
// All rights reserved.
// ExaGeoStat is a software package, provided by King Abdullah University of Science and Technology (KAUST).

/**
 * @file ExaGeoStat.cpp
 * @brief High-Level Wrapper class containing the static API for ExaGeoStat operations.
 * @version 1.1.0
 * @author Mahmoud ElKarargy
 * @date 2024-02-04
**/

#include <api/ExaGeoStat.hpp>
#include <data-generators/DataGenerator.hpp>
#include <prediction/Prediction.hpp>
#include <runtime-solver/RuntimeSolverFactory.hpp>

using namespace std;

using namespace exageostat::api;
using namespace exageostat::generators;
using namespace exageostat::dataunits;
using namespace exageostat::configurations;

template<typename T>
void ExaGeoStat<T>::ExaGeoStatLoadData(Configurations &aConfigurations, std::unique_ptr<ExaGeoStatData<T>> &aData) {

    aConfigurations.PrintSummary();
    LOGGER("** ExaGeoStat data generation/loading **")
    // Register and create a kernel object
    kernels::Kernel<T> *pKernel = plugins::PluginRegistry<kernels::Kernel<T>>::Create(aConfigurations.GetKernelName(),
                                                                                      aConfigurations.GetTimeSlot());
    // Add the data generation arguments.
    aConfigurations.InitializeDataGenerationArguments();
    // Create a unique pointer to a DataGenerator object
    unique_ptr<DataGenerator<T>> data_generator = DataGenerator<T>::CreateGenerator(aConfigurations);
    aData = data_generator->CreateData(aConfigurations, *pKernel);
    delete pKernel;
    LOGGER("\t*Data generation/loading finished*")

}

template<typename T>
T ExaGeoStat<T>::ExaGeoStatDataModeling(Configurations &aConfigurations, std::unique_ptr<ExaGeoStatData<T>> &aData,
                                        T *apMeasurementsMatrix) {

    aConfigurations.PrintSummary();
    LOGGER("** ExaGeoStat data Modeling **")
    // Register and create a kernel object
    kernels::Kernel<T> *pKernel = plugins::PluginRegistry<kernels::Kernel<T>>::Create(aConfigurations.GetKernelName(),
                                                                                      aConfigurations.GetTimeSlot());
    // Add the data modeling arguments.
    aConfigurations.InitializeDataModelingArguments();

    // We do Date Modeling with any computation.
    auto runtime_solver = runtimesolver::RuntimeSolverFactory<T>::CreateRuntimeSolver();
    T result = runtime_solver->ModelingOperations(aData, aConfigurations, apMeasurementsMatrix, *pKernel);
    delete pKernel;
    return result;
}


template<typename T>
void ExaGeoStat<T>::ExaGeoStatPrediction(Configurations &aConfigurations, std::unique_ptr<ExaGeoStatData<T>> &aData,
                                         T *apMeasurementsMatrix, Locations<T> *apTrainLocations,
                                         Locations<T> *apTestLocations) {

    aConfigurations.PrintSummary();
    LOGGER("** ExaGeoStat data Prediction **")
    // Register and create a kernel object
    kernels::Kernel<T> *pKernel = plugins::PluginRegistry<kernels::Kernel<T>>::Create(aConfigurations.GetKernelName(),
                                                                                      aConfigurations.GetTimeSlot());
    // Add the data prediction arguments.
    aConfigurations.InitializeDataPredictionArguments();
    prediction::Prediction<T>::PredictMissingData(aData, aConfigurations, apMeasurementsMatrix, *pKernel,
                                                  apTrainLocations, apTestLocations);
    delete pKernel;
}

