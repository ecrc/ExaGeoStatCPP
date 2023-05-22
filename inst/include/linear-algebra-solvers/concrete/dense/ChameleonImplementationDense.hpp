
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
                 * Initializes the descriptors needed for the Chameleon solver.
                 */
                void InitiateDescriptors() override;

                void
                CovarianceMatrixCodelet(void *descA, int uplo, dataunits::Locations *apLocation1, dataunits::Locations *apLocation2,
                                        dataunits::Locations *apLocation3, std::vector<double> aLocalTheta, int aDistanceMetric,
                                        exageostat::kernels::Kernel *apKernel) override;

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

                void EXAGEOSTAT_Zcpy(CHAM_desc_t *apDescA, double* apDoubleVector, RUNTIME_sequence_t *apSequence, RUNTIME_request_t *apRequest);
            };



            EXAGEOSTAT_INSTANTIATE_CLASS(ChameleonImplementationDense)

        }//namespace dense
    }//namespace linearAlgebra
}//namespace exageostat

#endif //EXAGEOSTATCPP_CHAMELEONIMPLEMENTATIONDENSE_HPP