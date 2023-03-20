
// Copyright (c) 2017-2023 King Abdullah University of Science and Technology,
// Copyright (C) 2023 by Brightskies inc,
// All rights reserved.
// ExaGeoStat is a software package, provided by King Abdullah University of Science and Technology (KAUST).

/**
 * @file AllocateDescriptors.hpp
 *
 * @version 1.0.0
 * @author Sameh Abdulah
 * @date 2023-03-20
**/

#ifndef EXAGEOSTATCPP_ALLOCATEDESCRIPTORS_HPP
#define EXAGEOSTATCPP_ALLOCATEDESCRIPTORS_HPP

namespace exageostat {
    namespace linearAlgebra {

        class AllocateDescriptors {
        public:
            virtual void InitiateDescriptors() = 0;

        private:
        };
    }//namespace linearAlgebra
}//namespace exageostat

#endif //EXAGEOSTATCPP_ALLOCATEDESCRIPTORS_HPP
