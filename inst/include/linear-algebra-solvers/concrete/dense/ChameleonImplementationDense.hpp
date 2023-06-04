
// Copyright (c) 2017-2023 King Abdullah University of Science and Technology,
// Copyright (C) 2023 by Brightskies inc,
// All rights reserved.
// ExaGeoStat is a software package, provided by King Abdullah University of Science and Technology (KAUST).

/**
 * @file ChameleonAllocateDescriptors.hpp
 * @brief This file contains the declaration of ChameleonImplementationDense class.
 * ChameleonImplementationDense is a concrete implementation of LinearAlgebraMethods class for dense matrices.
 * @version 1.0.0
 * @date 2023-03-20
**/

#ifndef EXAGEOSTATCPP_CHAMELEONIMPLEMENTATIONDENSE_HPP
#define EXAGEOSTATCPP_CHAMELEONIMPLEMENTATIONDENSE_HPP

#include <linear-algebra-solvers/LinearAlgebraMethods.hpp>

namespace exageostat {
    namespace linearAlgebra {
        namespace dense {

            /**
             * @brief
             * ChameleonImplementationDense is a concrete implementation of LinearAlgebraMethods class for dense matrices.
             * @tparam T Type of matrix elements.
             */
            template<typename T>
            class ChameleonImplementationDense : public LinearAlgebraMethods<T> {
            public:

                /**
                 * @brief
                 * Default constructor.
                 */
                explicit ChameleonImplementationDense() = default;

                /**
                 * @brief
                 * Virtual destructor to allow calls to the correct concrete destructor.
                 */
                virtual ~ChameleonImplementationDense() = default;

                /**
                 * @brief
                 * Initializes the descriptors needed for the Chameleon solver.
                 */
                void InitiateDescriptors() override;
                void DestoryDescriptors() override;
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
                void
                CovarianceMatrixCodelet(void *descA, int uplo, dataunits::Locations *apLocation1, dataunits::Locations *apLocation2,
                                        dataunits::Locations *apLocation3, double* aLocalTheta, int aDistanceMetric,
                                        exageostat::kernels::Kernel *apKernel) override;

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
                void GenerateObservationsVector(void *descA, dataunits::Locations *apLocation1, dataunits::Locations *apLocation2,
                                                     dataunits::Locations *apLocation3, std::vector<double> aLocalTheta, int aDistanceMetric, exageostat::kernels::Kernel * apKernel) override;
                /**
                 * @brief
                 * Initializes the context needed for the Chameleon solver.
                 *
                 * @param apCoresNumber Number of cores to allocate.
                 * @param apGPUs Number of GPUs to allocate.
                 */
                void ExaGeoStatInitContext(const int &apCoresNumber, const int &apGPUs) override;

                /**
                 * @brief
                 * Finalizes the context needed for the Chameleon solver.
                 */
                void ExaGeoStatFinalizeContext() override;

                void EXAGEOSTAT_Zcpy(CHAM_desc_t *apDescA, double* apDoubleVector, RUNTIME_sequence_t *apSequence, RUNTIME_request_t *apRequest);

            private:
                //// Used context
                static void *apContext;
            };

            /**
             * @brief Instantiates the LinearAlgebraMethods class for float and double types.
             */
            EXAGEOSTAT_INSTANTIATE_CLASS(ChameleonImplementationDense)
        }//namespace dense
    }//namespace linearAlgebra
}//namespace exageostat

#endif //EXAGEOSTATCPP_CHAMELEONIMPLEMENTATIONDENSE_HPP