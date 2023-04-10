
// Copyright (c) 2017-2023 King Abdullah University of Science and Technology,
// Copyright (C) 2023 by Brightskies inc,
// All rights reserved.
// ExaGeoStat is a software package, provided by King Abdullah University of Science and Technology (KAUST).

/**
 * @file ChameleonAllocateDescriptors.hpp
 *
 * @version 1.0.0
 * @author Sameh Abdulah
 * @date 2023-03-20
**/

#ifndef EXAGEOSTATCPP_CHAMELEONIMPLEMENTATIONDENSE_HPP
#define EXAGEOSTATCPP_CHAMELEONIMPLEMENTATIONDENSE_HPP

#include <linear-algebra-solvers/LinearAlgebraMethods.hpp>

namespace exageostat {
    namespace linearAlgebra {
        namespace dense {

            template<typename T>
            class ChameleonImplementationDense : public LinearAlgebraMethods<T>{
            public:

                void InitiateDescriptors() override;
                void ExaGeoStatInitContext(const int &apCoresNumber, const int &apGPUs) override;
                void ExaGeoStatFinalizeContext() override;

                /**
                 * @brief
                 * Default constructor.
                 *
                 */
                explicit ChameleonImplementationDense() = default;

                /**
                 * @brief
                 * Virtual destructor to allow calls to the correct concrete destructor.
                 *
                 */
                virtual ~ChameleonImplementationDense() = default;

            };

            EXAGEOSTAT_INSTANTIATE_CLASS(ChameleonImplementationDense)

        }//namespace dense
    }//namespace linearAlgebra
}//namespace exageostat

#endif //EXAGEOSTATCPP_CHAMELEONIMPLEMENTATIONDENSE_HPP
