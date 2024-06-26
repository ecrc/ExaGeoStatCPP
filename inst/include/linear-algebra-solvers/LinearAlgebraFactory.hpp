
// Copyright (c) 2017-2024 King Abdullah University of Science and Technology,
// All rights reserved.
// ExaGeoStat is a software package, provided by King Abdullah University of Science and Technology (KAUST).

/**
 * @file LinearAlgebraFactory.hpp
 * @brief Header file for the LinearAlgebraFactory class, which creates linear algebra solvers based on the input computation type.
 * @version 1.1.0
 * @author Mahmoud ElKarargy
 * @date 2023-03-20
**/

#ifndef EXAGEOSTATCPP_LINEARALGEBRAFACTORY_HPP
#define EXAGEOSTATCPP_LINEARALGEBRAFACTORY_HPP

#include <memory>

#include <common/Definitions.hpp>
#include <data-units/DescriptorData.hpp>
#include <linear-algebra-solvers/LinearAlgebraMethods.hpp>

namespace exageostat::linearAlgebra {

    /**
     * @class LinearAlgebraFactory
     * @brief A class that creates linear algebra solvers based on the input computation type.
     * @tparam T Data Type: float or double.
     *
     */
    template<typename T>
    class LinearAlgebraFactory {
    public:

        /**
         * @brief Creates a linear algebra solver based on the input computation type.
         * @param[in] aComputation The computation type to create the solver for.
         * @return Pointer to the created linear algebra solver.
         *
         */
        static std::unique_ptr<LinearAlgebraMethods<T>> CreateLinearAlgebraSolver(common::Computation aComputation);
    };

    /**
    * @brief Instantiates the Linear Algebra factory class for float and double types.
    * @tparam T Data Type: float or double
    *
    */
    EXAGEOSTAT_INSTANTIATE_CLASS(LinearAlgebraFactory)

}//namespace exageostat

#endif //EXAGEOSTATCPP_LINEARALGEBRAFACTORY_HPP