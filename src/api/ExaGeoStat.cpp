
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
#include <data-units/ModelingDataHolders.hpp>
#include <prediction/Prediction.hpp>
#include <data-transformer/DataTransformer.hpp>

using namespace std;
using namespace nlopt;

using namespace exageostat::api;
using namespace exageostat::generators;
using namespace exageostat::dataunits;
using namespace exageostat::configurations;
using namespace exageostat::transformers;

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

    int parameters_number = pKernel->GetParametersNumbers();
    int max_number_of_iterations = aConfigurations.GetMaxMleIterations();
    // Setting struct of data to pass to the modeling.
    auto modeling_data = new mModelingData(aData, aConfigurations, *apMeasurementsMatrix, *pKernel);
    // Create nlopt
    double opt_f;
    opt optimizing_function(nlopt::LN_BOBYQA, parameters_number);
    // Initialize problem's bound.
    optimizing_function.set_lower_bounds(aConfigurations.GetLowerBounds());
    optimizing_function.set_upper_bounds(aConfigurations.GetUpperBounds());
    optimizing_function.set_ftol_abs(aConfigurations.GetTolerance());
    // Set max iterations value.
    optimizing_function.set_maxeval(max_number_of_iterations);
    // TODO: ON API level it should be consistent regardless of the runtime being implemented
#if DEFAULT_RUNTIME
    optimizing_function.set_max_objective(ModelingAPI, (void *) modeling_data);
#endif
    // Optimize mle using nlopt.
    optimizing_function.optimize(aConfigurations.GetStartingTheta(), opt_f);
    aConfigurations.SetEstimatedTheta(aConfigurations.GetStartingTheta());

    auto theta = aConfigurations.GetStartingTheta();

    LOGGER("--> Final Theta Values (", true)
    for (int i = 0; i < parameters_number; i++) {
        LOGGER_PRECISION(theta[i])
        if (i != parameters_number - 1) {
            LOGGER_PRECISION(", ")
        }
    }
    LOGGER_PRECISION(")")
    LOGGER("")

    delete pKernel;
    delete modeling_data;
    return optimizing_function.last_optimum_value();
}

template<typename T>
double ExaGeoStat<T>::ModelingAPI(const std::vector<double> &aTheta, std::vector<double> &aGrad, void *apInfo) {

    auto config = ((mModelingData<T> *) apInfo)->mpConfiguration;
    auto data = ((mModelingData<T> *) apInfo)->mpData;
    auto measurements = ((mModelingData<T> *) apInfo)->mpMeasurementsMatrix;
    auto kernel = ((mModelingData<T> *) apInfo)->mpKernel;

    // We do Date Modeling with any computation.
    auto linear_algebra_solver = linearAlgebra::LinearAlgebraFactory<T>::CreateLinearAlgebraSolver(
            config->GetComputation());

    return linear_algebra_solver->ModelingOperations(*data, *config, aTheta.data(), measurements, *kernel);
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

#if !DEFAULT_RUNTIME
template<typename T>
void ExaGeoStat<T>::ExaGeoStatTransformData(Configurations &aConfigurations, std::unique_ptr<ExaGeoStatData<T>> &aData){

    DataTransformer<T>::ForwardSphericalHarmonicsTransform(aConfigurations.GetDenseTileSize(), aData);
    DataTransformer<T>::ForwardReshape(aConfigurations, aData);
    DataTransformer<T>::InverseSphericalHarmonicsTransform(aConfigurations.GetDenseTileSize(), aData);
}

template<typename T>
void ExaGeoStat<T>::ExaGeoStatDataAnalyzer(Configurations &aConfigurations, std::unique_ptr<ExaGeoStatData<T>> &aData){

}
#endif