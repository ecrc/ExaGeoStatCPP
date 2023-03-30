
// Copyright (c) 2017-2023 King Abdullah University of Science and Technology,
// Copyright (C) 2023 by Brightskies inc,
// All rights reserved.
// ExaGeoStat is a software package, provided by King Abdullah University of Science and Technology (KAUST).

/**
 * @file ChameleonImplementation.hpp
 *
 * @version 1.0.0
 * @author Sameh Abdulah
 * @date 2023-03-26
**/

#ifndef EXAGEOSTATCPP_CHAMELEONIMPLEMENTATIONDST_HPP
#define EXAGEOSTATCPP_CHAMELEONIMPLEMENTATIONDST_HPP

#include <linear-algebra-solvers/LinearAlgebraMethods.hpp>

// SUPPORT ONLY DOUBLE FOR NOW.
namespace exageostat {
    namespace linearAlgebra {
        namespace diagonalSuperTile {

            template<typename T>
            class ChameleonImplementationDST : public LinearAlgebraMethods<T>{
            public:

                void InitiateDescriptors() override;

                void ExaGeoStatInitContext(const int &apCoresNumber, const int &apGPUs) override;
                /**
                 * @brief
                 * Default constructor.
                 *
                 */
                explicit ChameleonImplementationDST() = default;

                /**
                 * @brief
                 * Virtual destructor to allow calls to the correct concrete destructor.
                 *
                 */
                virtual ~ChameleonImplementationDST() = default;

            };

            EXAGEOSTAT_INSTANTIATE_CLASS(ChameleonImplementationDST)

        }//namespace dense
    }//namespace linearAlgebra
}//namespace exageostat

#endif //EXAGEOSTATCPP_CHAMELEONIMPLEMENTATIONDST_HPP
