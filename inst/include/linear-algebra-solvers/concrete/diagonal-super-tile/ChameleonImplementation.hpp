
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

#ifndef EXAGEOSTATCPP_CHAMELEONIMPLEMENTATION_HPP
#define EXAGEOSTATCPP_CHAMELEONIMPLEMENTATION_HPP

#include <linear-algebra-solvers/LinearAlgebraMethods.hpp>
// SUPPORT ONLY DOUBLE FOR NOW.
namespace exageostat {
    namespace linearAlgebra {
        namespace diagonalSuperTile {

            template<typename T>
            class ChameleonImplementation : public LinearAlgebraMethods<T>{
            public:

                void InitiateDescriptors() override;

                /**
                 * @brief
                 * Default constructor.
                 *
                 */
                explicit ChameleonImplementation() = default;

                /**
                 * @brief
                 * Virtual destructor to allow calls to the correct concrete destructor.
                 *
                 */
                virtual ~ChameleonImplementation() = default;

            };

            EXAGEOSTAT_INSTANTIATE_CLASS(ChameleonImplementation)

        }//namespace dense
    }//namespace linearAlgebra
}//namespace exageostat

#endif //EXAGEOSTATCPP_CHAMELEONIMPLEMENTATION_HPP
