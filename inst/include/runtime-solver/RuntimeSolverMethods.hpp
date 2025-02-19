
// Copyright (c) 2017-2024 King Abdullah University of Science and Technology,
// All rights reserved.
// ExaGeoStat is a software package, provided by King Abdullah University of Science and Technology (KAUST).

/**
 * @file RuntimeSolverMethods.hpp
 * @brief Header file for the RuntimeSolverMethods class, which defines the interface for runtime solvers.
 * @version 2.0.0
 * @author Mahmoud ElKarargy
 * @author Sameh Abdulah
 * @date 2024-11-04
**/

#ifndef EXAGEOSTATCPP_RUNTIMESOLVERMETHODS_HPP
#define EXAGEOSTATCPP_RUNTIMESOLVERMETHODS_HPP

#include <data-units/ExaGeoStatData.hpp>
#include <kernels/Kernel.hpp>

namespace exageostat::runtimesolver {

    /**
     * @class RuntimeSolverMethods
     * @brief A class that defines the interface for linear algebra solvers.
     * @tparam T Data Type: float or double.
     *
     */
    template<typename T>
    class RuntimeSolverMethods {
    public:

        /**
         * @brief Virtual destructor to allow calls to the correct concrete destructor.
         *
         */
        virtual ~RuntimeSolverMethods() = default;

        /**
         * @brief The Gateway for the Modeling Operation
         * @param[in,out] aData DescriptorData object to be populated with descriptors and data.
         * @param[in] aConfigurations Configurations object containing relevant settings.
         * @param[in] apMeasurementsMatrix measurements matrix to be stored in DescZ.
         * @param[in] aKernel Reference to the kernel object to use.
         * @return log likelihood value
         *
         */
        virtual T ModelingOperations(std::unique_ptr<ExaGeoStatData<T>> &aData, configurations::Configurations &aConfigurations,
                                     T *apMeasurementsMatrix, const kernels::Kernel<T> &aKernel) = 0;

    };
}//namespace exageostat

#endif //EXAGEOSTATCPP_RUNTIMESOLVERMETHODS_HPP