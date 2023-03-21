
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

#include <memory>
#include <common/Definitions.hpp>

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
            static std::unique_ptr<LinearAlgebraFactory<T>>
            CreateLinearAlgebraSolver(common::Computation aComputation);
        };

        EXAGEOSTAT_INSTANTIATE_CLASS(LinearAlgebraFactory)

    }//namespace linearAlgebra
}//namespace exageostat



#endif //EXAGEOSTATCPP_LINEARALGEBRAFACTORY_HPP
