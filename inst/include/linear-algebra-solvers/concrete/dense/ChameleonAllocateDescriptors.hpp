
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

#ifndef EXAGEOSTATCPP_CHAMELEONALLOCATEDESCRIPTORS_HPP
#define EXAGEOSTATCPP_CHAMELEONALLOCATEDESCRIPTORS_HPP

#include <linear-algebra-solvers/AllocateDescriptors.hpp>

namespace exageostat {
    namespace linearAlgebra {
        namespace dense {
            class ChameleonAllocateDescriptors : public AllocateDescriptors{
            public:
                void InitiateDescriptors() override;

                /**
                 * @brief
                 * Default constructor.
                 *
                 */
                ChameleonAllocateDescriptors();

                /**
                 * @brief
                 * Virtual destructor to allow calls to the correct concrete destructor.
                 *
                 */
                virtual ~ChameleonAllocateDescriptors() = default;

            private:
            };
        }//namespace dense
    }//namespace linearAlgebra
}//namespace exageostat

#endif //EXAGEOSTATCPP_CHAMELEONALLOCATEDESCRIPTORS_HPP
