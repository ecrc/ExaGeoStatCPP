
// Copyright (c) 2017-2023 King Abdullah University of Science and Technology,
// Copyright (C) 2023 by Brightskies inc,
// All rights reserved.
// ExaGeoStat is a software package, provided by King Abdullah University of Science and Technology (KAUST).

/**
 * @file LinearAlgebraFactory.hpp
 *
 * @version 1.0.0
 * @author Sameh Abdulah
 * @date 2023-03-20
**/

#ifndef EXAGEOSTATCPP_LINEARALGEBRAFACTORY_HPP
#define EXAGEOSTATCPP_LINEARALGEBRAFACTORY_HPP

#include <common/Definitions.hpp>
#include <linear-algebra-solvers/LinearAlgebraMethods.hpp>
#include <memory>

#ifdef EXAGEOSTAT_USE_CHAMELEON
#include <linear-algebra-solvers/concrete/dense/ChameleonImplementationDense.hpp>
#include <linear-algebra-solvers/concrete/diagonal-super-tile/ChameleonImplementationDST.hpp>
#endif
#ifdef EXAGEOSTAT_USE_HiCMA
#include <linear-algebra-solvers/concrete/tile-low-rank/HicmaImplementation.hpp>
#endif

namespace exageostat {
    namespace linearAlgebra {

        template<typename T>
        class LinearAlgebraFactory {
        public:

            /**
             * @brief
             * Factory creation, Whether it's Synthetic or Real data.
             *
             * @param[in] apConfigurations
             *  Pointer to Synthetic data Configurations.
             *
             * @return DataGenerator
             * Unique Pointer to the created type of Data Generators.
             */
            static std::unique_ptr<LinearAlgebraMethods<T>>
            CreateLinearAlgebraSolver(common::Computation aComputation);
        };

        EXAGEOSTAT_INSTANTIATE_CLASS(LinearAlgebraFactory)

    }//namespace linearAlgebra
}//namespace exageostat



#endif //EXAGEOSTATCPP_LINEARALGEBRAFACTORY_HPP
