
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

            template<typename T>
            class ChameleonAllocateDescriptors : public AllocateDescriptors<T>{
            public:

                void InitiateDescriptors(dataunits::Precision aPrecision) override;

                /**
                 * @brief
                 * Default constructor.
                 *
                 */
                explicit ChameleonAllocateDescriptors();

                /**
                 * @brief
                 * Virtual destructor to allow calls to the correct concrete destructor.
                 *
                 */
                virtual ~ChameleonAllocateDescriptors() = default;

                void CreateDescriptors(T aPrecision);

            private:
                /// Used Precision.
                dataunits::Precision mPrecision;
            };
            template class ChameleonAllocateDescriptors<float>;
            template class ChameleonAllocateDescriptors<double>;
        }//namespace dense
    }//namespace linearAlgebra
}//namespace exageostat

#endif //EXAGEOSTATCPP_CHAMELEONALLOCATEDESCRIPTORS_HPP
