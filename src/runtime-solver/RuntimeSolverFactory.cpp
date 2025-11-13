
// Copyright (c) 2017-2024 King Abdullah University of Science and Technology,
// All rights reserved.
// ExaGeoStat is a software package, provided by King Abdullah University of Science and Technology (KAUST).

/**
 * @file RuntimeSolverFactory.cpp
 * @brief Implementation of the RuntimeSolverFactory class for creating runtime solvers for different runtime systems using StarPU or PaRSEC libraries.
 * The factory creates a unique pointer to a concrete implementation of the RuntimeSolverMethods class based on the runtime specified.
 * @version 2.0.0
 * @author Mahmoud ElKarargy
 * @date 2024-11-04
**/

#include <runtime-solver/RuntimeSolverFactory.hpp>

#if DEFAULT_RUNTIME
#include <runtime-solver/concrete/StarpuRuntimeSolver.hpp>
#else
#include <runtime-solver/concrete/ParsecRuntimeSolver.hpp>
#endif

using namespace exageostat::runtimesolver;
using namespace exageostat::common;

template<typename T>
std::unique_ptr<RuntimeSolverMethods<T>> RuntimeSolverFactory<T>::CreateRuntimeSolver() {

    // Check which Runtime is used
#if DEFAULT_RUNTIME
    return std::make_unique<StarpuRuntimeSolver<T>>();
#else
    return std::make_unique<ParsecRuntimeSolver<T>>();
#endif

}
