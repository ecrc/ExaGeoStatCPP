
// Copyright (c) 2017-2024 King Abdullah University of Science and Technology,
// All rights reserved.
// ExaGeoStat is a software package, provided by King Abdullah University of Science and Technology (KAUST).

/**
 * @file StarpuRuntimeSolver.hpp
 * @brief This file contains the declaration of StarpuRuntimeSolver class.
 * @details StarpuRuntimeSolver is a concrete implementation of the RuntimeSolverMethods class.
 * @version 2.0.0
 * @author Mahmoud ElKarargy
 * @author Sameh Abdulah
 * @date 2024-11-04
**/

#ifndef EXAGEOSTATCPP_STARPURUNTIMESOLVER_HPP
#define EXAGEOSTATCPP_STARPURUNTIMESOLVER_HPP

#include <runtime-solver/RuntimeSolverMethods.hpp>
#include <nlopt.hpp>

namespace exageostat::runtimesolver {

    /**
     * @brief StarpuRuntimeSolver is a concrete implementation of RuntimeSolverMethods class for dense or diagonal-super tile matrices.
     * @tparam T Data Type: float or double
     *
     */
    template<typename T>
    class StarpuRuntimeSolver : public RuntimeSolverMethods<T> {
    public:

        /**
         * @brief Calculates the log likelihood value of a given value theta.
         * @copydoc RuntimeSolverMethods::ModelingOperations()
         *
         */
        T ModelingOperations(std::unique_ptr <ExaGeoStatData<T>> &aData,
                             configurations::Configurations &aConfigurations, T *apMeasurementsMatrix, const kernels::Kernel <T> &aKernel) override;


        /**
         * @brief Objective function used in optimization, and following the NLOPT objective function format.
         * @param[in] aTheta An array of length n containing the current point in the parameter space.
         * @param[in] aGrad  An array of length n where you can optionally return the gradient of the objective function.
         * @param[in] apInfo pointer containing needed configurations and data.
         * @return double MLE results.
         *
         */
        static double DataModelingAPI(const std::vector<double> &aTheta, std::vector<double> &aGrad, void *apInfo);

    };

    /**
    * @brief Instantiates the Starpu Runtime Solver class for float and double types.
    * @tparam T Data Type: float or double
    *
    */
    EXAGEOSTAT_INSTANTIATE_CLASS(StarpuRuntimeSolver)
}//namespace exageostat

#endif //EXAGEOSTATCPP_STARPURUNTIMESOLVER_HPP