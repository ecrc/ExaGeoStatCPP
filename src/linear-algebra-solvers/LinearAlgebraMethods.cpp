
// Copyright (c) 2017-2023 King Abdullah University of Science and Technology,
// Copyright (C) 2023 by Brightskies inc,
// All rights reserved.
// ExaGeoStat is a software package, provided by King Abdullah University of Science and Technology (KAUST).

/**
 * @file AllocateDescriptors.cpp
 *
 * @version 1.0.0
 * @author Sameh Abdulah
 * @date 2023-03-20
**/

#include <linear-algebra-solvers/LinearAlgebraMethods.hpp>
#ifdef EXAGEOSTAT_USE_CHAMELEON
    #include <chameleon/struct.h>
    #include <chameleon.h>
    #include <control/context.h>
#endif

using namespace exageostat::linearAlgebra;
using namespace exageostat::configurations;

template<typename T>
void LinearAlgebraMethods<T>::SetConfigurations(Configurations *apConfigurations) {
    this->mpConfigurations = apConfigurations;
}

//template<typename T>
//void LinearAlgebraMethods<T>::ExaGeoStatInitContext(int *apCoresNumber, int *apGPUs) {
//
//#ifdef EXAGEOSTAT_USE_CHAMELEON
//    CHAM_context_t *chamctxt;
//    chamctxt = chameleon_context_self();
//    if (chamctxt != NULL) {
//        printf("Another instance of Chameleon is already running...!");
//    } else {
//        CHAMELEON_user_tag_size(31, 26);
//        CHAMELEON_Init(*apCoresNumber, *apGPUs);
//    }
//#endif
//}