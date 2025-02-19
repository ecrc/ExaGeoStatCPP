
// Copyright (c) 2017-2024 King Abdullah University of Science and Technology,
// All rights reserved.
// ExaGeoStat is a software package, provided by King Abdullah University of Science and Technology (KAUST).

/**
 * @file StarpuRuntimeSolver.cpp
 * @brief This file contains the implementation of StarpuRuntimeSolver class.
 * @details StarpuRuntimeSolver is a concrete implementation of the RuntimeSolversMethods class.
 * @version 2.0.0
 * @author Mahmoud ElKarargy
 * @author Sameh Abdulah
 * @author Qinglei Cao
 * @date 2024-11-04
**/

#include <runtime-solver/concrete/StarpuRuntimeSolver.hpp>
#include <data-units/ModelingDataHolders.hpp>
#include <utilities/Logger.hpp>
#include <linear-algebra-solvers/LinearAlgebraFactory.hpp>

using namespace exageostat::common;
using namespace exageostat::configurations;
using namespace exageostat::runtimesolver;
using namespace exageostat::dataunits;
using namespace nlopt;

template<typename T>
T StarpuRuntimeSolver<T>::ModelingOperations(std::unique_ptr<ExaGeoStatData<T>> &aData, Configurations &aConfigurations,
                                              T *apMeasurementsMatrix, const kernels::Kernel<T> &aKernel) {

    int parameters_number = aKernel.GetParametersNumbers();
    int max_number_of_iterations = aConfigurations.GetMaxMleIterations();
    // Setting struct of data to pass to the modeling.
    auto modeling_data = new mModelingData(aData, aConfigurations, *apMeasurementsMatrix, aKernel);
    // Create nlopt
    double opt_f;
    opt optimizing_function(nlopt::LN_BOBYQA, parameters_number);
    // Initialize problem's bound.
    optimizing_function.set_lower_bounds(aConfigurations.GetLowerBounds());
    optimizing_function.set_upper_bounds(aConfigurations.GetUpperBounds());
    optimizing_function.set_ftol_abs(aConfigurations.GetTolerance());
    // Set max iterations value.
    optimizing_function.set_maxeval(max_number_of_iterations);
    optimizing_function.set_max_objective(DataModelingAPI, (void *) modeling_data);

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

    delete modeling_data;
    return optimizing_function.last_optimum_value();

}

template<typename T>
double StarpuRuntimeSolver<T>::DataModelingAPI(const std::vector<double> &aTheta, std::vector<double> &aGrad, void *apInfo) {

    auto config = ((mModelingData<T> *) apInfo)->mpConfiguration;
    auto data = ((mModelingData<T> *) apInfo)->mpData;
    auto measurements = ((mModelingData<T> *) apInfo)->mpMeasurementsMatrix;
    auto kernel = ((mModelingData<T> *) apInfo)->mpKernel;

    // We do Date Modeling with any computation.
    auto linear_algebra_solver = linearAlgebra::LinearAlgebraFactory<T>::CreateLinearAlgebraSolver(config->GetComputation());
    return linear_algebra_solver->ExaGeoStatMLETile(*data, *config, aTheta.data(), measurements, *kernel);
}
