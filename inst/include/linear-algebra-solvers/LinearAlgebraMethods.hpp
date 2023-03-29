
// Copyright (c) 2017-2023 King Abdullah University of Science and Technology,
// Copyright (C) 2023 by Brightskies inc,
// All rights reserved.
// ExaGeoStat is a software package, provided by King Abdullah University of Science and Technology (KAUST).

/**
 * @file AllocateDescriptors.hpp
 *
 * @version 1.0.0
 * @author Sameh Abdulah
 * @date 2023-03-20
**/

#ifndef EXAGEOSTATCPP_LINEARALGEBRAMETHODS_HPP
#define EXAGEOSTATCPP_LINEARALGEBRAMETHODS_HPP

#include <common/Definitions.hpp>
#include <configurations/Configurations.hpp>
#include <linear-algebra-solvers/concrete/MatrixAllocation.hpp>
#include <gsl/gsl_errno.h>
#include <vector>

namespace exageostat {
    namespace linearAlgebra {

        template<typename T>
        class LinearAlgebraMethods {
        public:
            virtual void InitiateDescriptors() = 0;
            virtual void ExaGeoStatInitContext(const int &apCoresNumber, const int &apGPUs) = 0;

            /**
             * @brief
             * Configuration map setter.
             *
             * @param apConfigurations
             * Argument pointer to Synthetic Data generation configuration map
             *
             */
            void
            SetConfigurations(configurations::Configurations *apConfigurations){
                this->mpConfigurations = apConfigurations;
            }

            /**
             * @brief
             * Virtual destructor to allow calls to the correct concrete destructor.
             *
             */
            virtual ~LinearAlgebraMethods() = default;

        protected:
            //// Used configurations map.
            configurations::Configurations *mpConfigurations = nullptr;
        };

        EXAGEOSTAT_INSTANTIATE_CLASS(LinearAlgebraMethods)

    }//namespace linearAlgebra
}//namespace exageostat

#endif //EXAGEOSTATCPP_LINEARALGEBRAMETHODS_HPP
