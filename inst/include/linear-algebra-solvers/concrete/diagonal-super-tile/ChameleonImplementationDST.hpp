
// Copyright (c) 2017-2023 King Abdullah University of Science and Technology,
// Copyright (C) 2023 by Brightskies inc,
// All rights reserved.
// ExaGeoStat is a software package, provided by King Abdullah University of Science and Technology (KAUST).

/**
 * @file ChameleonImplementation.hpp
 * @brief This file contains the declaration of ChameleonImplementationDST class.
 * ChameleonImplementationDST is a concrete implementation of LinearAlgebraMethods class for diagonal super tile matrices.
 * @version 1.0.0
 * @date 2023-03-26
**/

#ifndef EXAGEOSTATCPP_CHAMELEONIMPLEMENTATIONDST_HPP
#define EXAGEOSTATCPP_CHAMELEONIMPLEMENTATIONDST_HPP

#include <linear-algebra-solvers/LinearAlgebraMethods.hpp>

namespace exageostat {
    namespace linearAlgebra {
        namespace diagonalSuperTile {

            /**
             * @brief
             * ChameleonImplementationDST is a concrete implementation of LinearAlgebraMethods class for diagonal super tile matrices.
             *
             * @tparam T Type of matrix elements.
             */
            template<typename T>
            class ChameleonImplementationDST : public LinearAlgebraMethods<T>{
            public:

                /**
                 * @brief
                 * Initializes the descriptors needed for the Chameleon solver.
                 */
                void InitiateDescriptors() override;

                void *EXAGEOSTAT_DATA_GET_ADDRESS(const void *A, int type, int m, int n) override;

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
                explicit ChameleonImplementationDST() = default;

                /**
                 * @brief
                 * Virtual destructor to allow calls to the correct concrete destructor.
                 */
                virtual ~ChameleonImplementationDST() = default;

            };

            EXAGEOSTAT_INSTANTIATE_CLASS(ChameleonImplementationDST)

        }//namespace diagonalSuperTile
    }//namespace linearAlgebra
}//namespace exageostat

#endif //EXAGEOSTATCPP_CHAMELEONIMPLEMENTATIONDST_HPP