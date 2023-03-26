
// Copyright (c) 2017-2023 King Abdullah University of Science and Technology,
// Copyright (C) 2023 by Brightskies inc,
// All rights reserved.
// ExaGeoStat is a software package, provided by King Abdullah University of Science and Technology (KAUST).

/**
 * @file HicmaImplementation.hpp
 *
 * @version 1.0.0
 * @author Sameh Abdulah
 * @date 2023-03-26
**/

#ifndef EXAGEOSTATCPP_HICMAIMPLEMENTATION_HPP
#define EXAGEOSTATCPP_HICMAIMPLEMENTATION_HPP

#include <linear-algebra-solvers/LinearAlgebraMethods.hpp>

namespace exageostat {
    namespace linearAlgebra {
        namespace tileLowRank {

            template<typename T>
            class HicmaImplementation : public LinearAlgebraMethods<T>{
            public:

                void InitiateDescriptors() override;

                /**
                 * @brief
                 * Default constructor.
                 *
                 */
                explicit HicmaImplementation() = default;

                /**
                 * @brief
                 * Virtual destructor to allow calls to the correct concrete destructor.
                 *
                 */
                virtual ~HicmaImplementation() = default;

            };

            EXAGEOSTAT_INSTANTIATE_CLASS(HicmaImplementation);

        }//namespace dense
    }//namespace linearAlgebra
}//namespace exageostat

#endif //EXAGEOSTATCPP_HICMAIMPLEMENTATION_HPP
