
// Copyright (c) 2017-2023 King Abdullah University of Science and Technology,
// Copyright (C) 2023 by Brightskies inc,
// All rights reserved.
// ExaGeoStat is a software package, provided by King Abdullah University of Science and Technology (KAUST).

/**
 * @file AllocateDescriptors.hpp
 * @brief Header file for the LinearAlgebraMethods class, which defines the interface for linear algebra solvers.
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
#include <vector>
extern "C"{
#include <gsl/gsl_errno.h>
}

namespace exageostat {
    namespace linearAlgebra {

        /**
         * @class LinearAlgebraMethods
         * @brief A class that defines the interface for linear algebra solvers.
         * @tparam T The data type of the linear algebra solver.
         */
        template<typename T>
        class LinearAlgebraMethods {
        public:
            /**
             * @brief Initializes the descriptors necessary for the linear algebra solver.
             */
            virtual void InitiateDescriptors() = 0;

            /**
             * @brief Initializes the context for the linear algebra solver with the specified number of cores and GPUs.
             *
             * @param[in] apCoresNumber The number of cores to use for the solver.
             * @param[in] apGPUs The number of GPUs to use for the solver.
             */
            virtual void ExaGeoStatInitContext(const int &apCoresNumber, const int &apGPUs) = 0;

            /**
             * @brief Finalizes the context for the linear algebra solver.
             */
            virtual void ExaGeoStatFinalizeContext() = 0;
            virtual void testKernelfornow() = 0;

            /**
             * @brief Sets the configurations for the linear algebra solver.
             *
             * @param[in] apConfigurations A pointer to the configurations for the solver.
             */
            void SetConfigurations(configurations::Configurations *apConfigurations){
                this->mpConfigurations = apConfigurations;
            }

            /**
             * @brief Virtual destructor to allow calls to the correct concrete destructor.
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