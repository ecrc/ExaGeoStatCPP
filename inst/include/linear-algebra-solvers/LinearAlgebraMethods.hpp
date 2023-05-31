
// Copyright (c) 2017-2023 King Abdullah University of Science and Technology,
// Copyright (C) 2023 by Brightskies inc,
// All rights reserved.
// ExaGeoStat is a software package, provided by King Abdullah University of Science and Technology (KAUST).

/**
 * @file AllocateDescriptors.hpp
 * @brief Header file for the LinearAlgebraMethods class, which defines the interface for linear algebra solvers.
 * @version 1.0.0
 * @author Sameh Abdulah
 * @date 2023-03-20
 *
 *  This header file defines the abstract class LinearAlgebraMethods, which provides an interface for linear algebra solvers.
 *  The purpose of this interface is to allow different concrete linear algebra solvers to be interchangeable,
 *  so that they can be used interchangeably by other parts of the software system that rely on linear algebra.
 *
**/

#ifndef EXAGEOSTATCPP_LINEARALGEBRAMETHODS_HPP
#define EXAGEOSTATCPP_LINEARALGEBRAMETHODS_HPP

#include <common/Definitions.hpp>
#include <configurations/Configurations.hpp>
#include <linear-algebra-solvers/concrete/MatrixAllocation.hpp>
#include <vector>

extern "C" {
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
             * @brief Virtual destructor to allow calls to the correct concrete destructor.
             */
            virtual ~LinearAlgebraMethods() = default;

            /**
             * @brief Initializes the descriptors necessary for the linear algebra solver.
             */
            virtual void InitiateDescriptors() = 0;

            /**
             * @brief Computes the covariance matrix.
             *
             * @param[in] descA Pointer to the descriptor for the covariance matrix.
             * @param[in] uplo Specifies whether the upper or lower triangular part of the covariance matrix is stored.
             * @param[in] apLocation1 Pointer to the first set of locations.
             * @param[in] apLocation2 Pointer to the second set of locations.
             * @param[in] apLocation3 Pointer to the third set of locations.
             * @param[in] aLocalTheta Pointer to the local theta values.
             * @param[in] aDistanceMetric Specifies the distance metric to use.
             * @param[in] apKernel Pointer to the kernel function to use.
             */
            virtual void
            CovarianceMatrixCodelet(void *descA, int uplo, dataunits::Locations *apLocation1, dataunits::Locations *apLocation2,
                                    dataunits::Locations *apLocation3, double *aLocalTheta, int aDistanceMetric, exageostat::kernels::Kernel * apKernel) = 0;

            /**
             * @brief Generates the observations vector.
             *
             * @param[in] descA Pointer to the descriptor for the observations vector.
             * @param[in] apLocation1 Pointer to the first set of locations.
             * @param[in] apLocation2 Pointer to the second set of locations.
             * @param[in] apLocation3 Pointer to the third set of locations.
             * @param[in] aLocalTheta Pointer to the local theta values.
             * @param[in] aDistanceMetric Specifies the distance metric to use.
             * @param[in] apKernel Pointer to the kernel function to use.
             */
            virtual void GenerateObservationsVector(void *descA, dataunits::Locations *apLocation1, dataunits::Locations *apLocation2,
                                                    dataunits::Locations *apLocation3, std::vector<double> aLocalTheta, int aDistanceMetric, exageostat::kernels::Kernel * apKernel) = 0;
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

            /**
             * @brief Sets the configurations for the linear algebra solver.
             *
             * @param[in] apConfigurations A pointer to the configurations for the solver.
             */
            void SetConfigurations(configurations::Configurations *apConfigurations) {
                this->mpConfigurations = apConfigurations;
            }

            /**
             * @brief Gets the matrix.
             *
             * @return Pointer to the matrix.
             */
            double *GetMatrix() {
                return this->apMatrix;
            }


        protected:
            //// Used configurations map.
            configurations::Configurations *mpConfigurations = nullptr;
            //// Used Matrix
            double *apMatrix = nullptr;
        };

        /**
         * @brief Instantiates the LinearAlgebraMethods class for float and double types.
         */
        EXAGEOSTAT_INSTANTIATE_CLASS(LinearAlgebraMethods)

    }//namespace linearAlgebra
}//namespace exageostat

#endif //EXAGEOSTATCPP_LINEARALGEBRAMETHODS_HPP