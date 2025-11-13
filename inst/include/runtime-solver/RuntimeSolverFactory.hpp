
// Copyright (c) 2017-2024 King Abdullah University of Science and Technology,
// All rights reserved.
// ExaGeoStat is a software package, provided by King Abdullah University of Science and Technology (KAUST).

/**
 * @file RuntimeSolverFactory.hpp
 * @brief Header file for the RuntimeSolverFactory class, which creates runtime solvers based on the configured runtime.
 * @version 2.0.0
 * @author Mahmoud ElKarargy
 * @date 2024-11-04
**/

#ifndef EXAGEOSTATCPP_RUNTIMESOLVERFACTORY_HPP
#define EXAGEOSTATCPP_RUNTIMESOLVERFACTORY_HPP

#include <memory>

#include <common/Definitions.hpp>
#include <data-units/DescriptorData.hpp>
#include <runtime-solver/RuntimeSolverMethods.hpp>

namespace exageostat::runtimesolver {

    /**
     * @class RuntimeSolverFactory
     * @brief A class that creates linear algebra solvers based on the input computation type.
     * @tparam T Data Type: float or double.
     *
     */
    template<typename T>
    class RuntimeSolverFactory {
    public:

        /**
         * @brief Creates a linear algebra solver based on the input computation type.
         * @return Pointer to the created linear algebra solver.
         *
         */
        static std::unique_ptr<RuntimeSolverMethods<T>> CreateRuntimeSolver();
    };

    /**
    * @brief Instantiates the Runtime Solver Factory class for float and double types.
    * @tparam T Data Type: float or double
    *
    */
    EXAGEOSTAT_INSTANTIATE_CLASS(RuntimeSolverFactory)

}//namespace exageostat

#endif //EXAGEOSTATCPP_RUNTIMESOLVERFACTORY_HPP