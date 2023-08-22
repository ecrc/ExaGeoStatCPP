
// Copyright (c) 2017-2023 King Abdullah University of Science and Technology,
// All rights reserved.
// ExaGeoStat is a software package, provided by King Abdullah University of Science and Technology (KAUST).

/**
 * @file ExaGeoStatHardware.cpp
 * @brief Contains the implementation of the ExaGeoStatHardware class.
 * @version 1.0.0
 * @author Sameh Abdulah
 * @author Mahmoud ElKarargy
 * @date 2023-08-07
**/

#ifdef EXAGEOSTAT_USE_CHAMELEON
extern "C" {
#include <chameleon.h>
#include <control/context.h>
}
#endif

#ifdef EXAGEOSTAT_USE_HiCMA
extern "C"{
#include <hicma.h>
#include <control/hicma_context.h>
}
#endif

#include <hardware/ExaGeoStatHardware.hpp>

using namespace exageostat::hardware;

ExaGeoStatHardware::ExaGeoStatHardware(common::Computation aComputation, int aCoreNumber, int aGpuNumber) {

    this->mComputation = aComputation;

    // Init hardware using Hicma
    if (aComputation == common::TILE_LOW_RANK) {
#ifdef EXAGEOSTAT_USE_HiCMA
        if (!this->mpContext) {
            HICMA_user_tag_size(31, 26);
            HICMA_Init(aCoreNumber, aGpuNumber);
            this->mpContext = hicma_context_self();
        }
#else
        throw std::runtime_error("You need to enable Hicma to use TLR computation!");
#endif
    }
        // Init hardware using Chameleon
    else {
#ifdef EXAGEOSTAT_USE_CHAMELEON
        if (!this->mpContext) {
            CHAMELEON_user_tag_size(31, 26);
            CHAMELEON_Init(aCoreNumber, aGpuNumber)
            this->mpContext = chameleon_context_self();
        }
#else
        throw std::runtime_error("You need to enable Chameleon to use Dense or DST computation!");
#endif
    }
}

ExaGeoStatHardware::~ExaGeoStatHardware() {
    // Init hardware using Hicma
    if (this->mComputation == common::TILE_LOW_RANK) {
#ifdef EXAGEOSTAT_USE_HiCMA
        if (!this->mpContext) {
            std::cout << "No initialised context of HiCMA, Please use 'ExaGeoStatHardware::ExaGeoStatHardware(aComputation, CoreNumber, aGpuNumber);'" << std::endl;
        } else {
            HICMA_Finalize();
            this->mpContext = nullptr;
        }
#endif
    }
        // Init hardware using Chameleon
    else {
#ifdef EXAGEOSTAT_USE_CHAMELEON
        if (!this->mpContext) {
            std::cout << "No initialised context of Chameleon, Please Initialise a hardware first" << std::endl;
        } else {
            CHAMELEON_Finalize()
            this->mpContext = nullptr;
        }
#endif
    }
}

void *ExaGeoStatHardware::GetContext() {
    if (!this->mpContext) {
        throw std::runtime_error("Hardware is not initialised!");
    }
    return this->mpContext;
}
