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
#include <lapacke.h>

#include <nlopt.hpp>

#include <api/ExaGeoStat.hpp>
#include <linear-algebra-solvers/LinearAlgebraFactory.hpp>
#include <data-generators/DataGenerator.hpp>

using namespace std;

using namespace exageostat::api;
using namespace exageostat::configurations;
using namespace exageostat::linearAlgebra;
using namespace exageostat::generators;
using namespace exageostat::dataunits;

template<typename T>
void ExaGeoStat<T>::ExaGeoStatInitializeHardware(common::Computation aComputation, int aCoreNumber, int aGpuNumber) {

    auto linearAlgebraSolver = LinearAlgebraFactory<T>::CreateLinearAlgebraSolver(aComputation);
    linearAlgebraSolver->ExaGeoStatInitContext(aCoreNumber, aGpuNumber);
    delete linearAlgebraSolver;
}

template<typename T>
void ExaGeoStat<T>::ExaGeoStatFinalizeHardware(common::Computation aComputation, DescriptorData<T> *apDescriptorData) {

    auto linearAlgebraSolver = LinearAlgebraFactory<T>::CreateLinearAlgebraSolver(aComputation);
//    linearAlgebraSolver->DestroyDescriptors(apDescriptorData);
    linearAlgebraSolver->ExaGeoStatFinalizeContext();
    delete linearAlgebraSolver;
}

template<typename T>
void ExaGeoStat<T>::ExaGeoStatInitOptimizer(nlopt_opt *apOpt, double *apLowerBound, double *apUpperBound, double aTolerance)
{
    //initalizing opt library
    nlopt_set_lower_bounds(*apOpt, apLowerBound);
    nlopt_set_upper_bounds(*apOpt, apUpperBound);
    nlopt_set_ftol_abs(*apOpt, aTolerance);
}

template<typename T>
ExaGeoStatData<T> *ExaGeoStat<T>::ExaGeoStatGenerateData(Configurations *apConfigurations) {

    // Add the data generation arguments.
    apConfigurations->InitializeDataGenerationArguments();
    // Create a unique pointer to a DataGenerator object
    unique_ptr<DataGenerator<T>> data_generator;
    // Create the DataGenerator object
    data_generator = data_generator->CreateGenerator(apConfigurations);
    // Create a new object of ExaGeoStat Data to store the locations and Descriptors
    auto *data = new ExaGeoStatData<T>(apConfigurations->GetProblemSize(), apConfigurations->GetDimension());
    data_generator->GenerateLocations();
    // Copy the locations from generator to the data.
    auto pLocation = new Locations<T>(*data_generator->GetLocations());
    data->SetLocations(pLocation);
    data_generator->GetLocations()->SetLocationX(nullptr);


    auto linear_algebra_solver = LinearAlgebraFactory<T>::CreateLinearAlgebraSolver(apConfigurations->GetComputation());

    linear_algebra_solver->InitiateDescriptors(apConfigurations, data->GetDescriptorData());

#ifdef EXAGEOSTAT_USE_CHAMELEON
    BaseDescriptor descC = data->GetDescriptorData()->GetDescriptor(common::CHAMELEON_DESCRIPTOR, common::DESCRIPTOR_C);
#endif
#ifdef EXAGEOSTAT_USE_HiCMA
    auto *descC = data->GetDescriptorData()->GetDescriptor(common::HICMA_DESCRIPTOR, common::DESCRIPTOR_C);
#endif

    int iseed[4] = {apConfigurations->GetSeed(), apConfigurations->GetSeed(), apConfigurations->GetSeed(), 1};
    //nomral random generation of e -- ei~N(0, 1) to generate Z
    T* Nrand = (T* ) malloc(apConfigurations->GetProblemSize() * 1 * sizeof(T));
    LAPACKE_dlarnv(3, iseed, apConfigurations->GetProblemSize() * 1, (double *)Nrand);


    linear_algebra_solver->GenerateObservationsVector(apConfigurations, data->GetDescriptorData(), descC,
                                                      data->GetLocations(), data->GetLocations(), nullptr, 0,
                                                      data_generator->GetKernel(), &Nrand[0]);

    return data;
}

template<typename T>
double ExaGeoStat<T>::ExaGeoStatMleTileAPI(unsigned aN, const double *apTheta, double *apGrad, void *apInfo) {
    auto config = ((mModelingData *) apInfo)->Configuration;
    auto data = ((mModelingData *) apInfo)->Data;
    auto linear_algebra_solver = linearAlgebra::LinearAlgebraFactory<T>::CreateLinearAlgebraSolver(config->GetComputation());
    linear_algebra_solver->ExaGeoStatMleTile(data, config, apTheta);
}

template<typename T>
void ExaGeoStat<T>::ExaGeoStatDataModeling(Configurations *apConfigurations, ExaGeoStatData<T> *apData) {

    apConfigurations->InitializeDataModelingArguments();

    int max_number_of_iterations = apConfigurations->GetMaxMleIterations();

    // Add the data modeling arguments.
    nlopt_opt opt;

    //// @warning: this shouldn't be how the num of params passed.
    kernels::Kernel<T> *kernel = exageostat::plugins::PluginRegistry<kernels::Kernel<T>>::Create(
            apConfigurations->GetKernelName());

    opt = nlopt_create(NLOPT_LN_BOBYQA,  kernel->GetParametersNumbers());

    delete kernel;

    ExaGeoStatInitOptimizer(&opt, apConfigurations->GetLowerBounds().data(), apConfigurations->GetUpperBounds().data(),
                   apConfigurations->GetTolerance());
    nlopt_set_maxeval(opt, max_number_of_iterations);


    auto linear_algebra_solver = LinearAlgebraFactory<double>::CreateLinearAlgebraSolver(apConfigurations->GetComputation());
    linear_algebra_solver->ExaGeoStatSequenceWait(apData->GetDescriptorData()->GetSequence());
    delete linear_algebra_solver;

    PrintSummary(apConfigurations->GetProblemSize(), apConfigurations->GetCoresNumber(),
                  apConfigurations->GetGPUsNumbers(), apConfigurations->GetDenseTileSize(), apConfigurations->GetComputation(),
                  apConfigurations->GetPGrid(), apConfigurations->GetQGrid(), apConfigurations->GetPrecision());

    //// TODO: that's the most challenging part, read more about NLOPT and redesign this part
    auto modeling_data = new mModelingData();
    modeling_data->Configuration= apConfigurations;
    modeling_data->Data = apData;
    nlopt_set_max_objective(opt, ExaGeoStatMleTileAPI, (void *) modeling_data);

    double opt_f = 0;
    double * theta = apConfigurations->GetStartingTheta().data();
    nlopt_optimize(opt, theta, &opt_f);

    delete modeling_data;
    nlopt_destroy(opt);
}
