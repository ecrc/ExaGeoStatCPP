
// Copyright (c) 2017-2023 King Abdullah University of Science and Technology,
// All rights reserved.
// ExaGeoStat is a software package, provided by King Abdullah University of Science and Technology (KAUST).

/**
 * @file ExaGeoStatHardware.cpp
 * @brief Contains the implementation of the ExaGeoStatHardware class.
 * @version 1.0.1
 * @author Mahmoud ElKarargy
 * @author Sameh Abdulah
 * @date 2024-01-24
**/

#include <linear-algebra-solvers/concrete/ChameleonHeaders.hpp>
#include <linear-algebra-solvers/concrete/HicmaHeaders.hpp>
#include <hardware/ExaGeoStatHardware.hpp>
#include <results/Results.hpp>
#include <helpers/CommunicatorMPI.hpp>
#include <common/Utils.hpp>

using namespace exageostat::hardware;

ExaGeoStatHardware::ExaGeoStatHardware(const common::Computation &aComputation, const int &aCoreNumber,
                                       const int &aGpuNumber) {
    LOGGER("** Initialise ExaGeoStat hardware **")
    this->mComputation = aComputation;
    int tag_width = 31, tag_sep = 26;

    // Init hardware using Chameleon
    if (!mpChameleonContext) {
        CHAMELEON_user_tag_size(tag_width, tag_sep);
        CHAMELEON_Init(aCoreNumber, aGpuNumber)
        mpChameleonContext = chameleon_context_self();
    }

    // Init hardware using Hicma
    if (aComputation == common::TILE_LOW_RANK) {
#ifdef USE_HICMA
        if (!mpHicmaContext) {
            HICMA_user_tag_size(tag_width, tag_sep);
            HICMA_Init(aCoreNumber, aGpuNumber);
            mpHicmaContext = hicma_context_self();
        }
#else
        throw std::runtime_error("You need to enable HiCMA to use TLR computation!");
#endif
    }
    helpers::CommunicatorMPI::GetInstance()->SetHardwareInitialization();
}

ExaGeoStatHardware::~ExaGeoStatHardware() {
    // finalize hardware using Chameleon
    if (mpChameleonContext) {
        CHAMELEON_Finalize()
        mpChameleonContext = nullptr;
    }
    // finalize hardware using HiCMA
    if (this->mComputation == common::TILE_LOW_RANK) {
#ifdef USE_HICMA
        if (mpHicmaContext) {
            HICMA_Finalize();
            mpHicmaContext = nullptr;
        }
#endif
    }
    helpers::CommunicatorMPI::GetInstance()->RemoveHardwareInitialization();
    results::Results::GetInstance()->PrintEndSummary();
}

void *ExaGeoStatHardware::GetHicmaContext() {
    if (!mpHicmaContext) {
        throw std::runtime_error("HiCMA Hardware is not initialized!");
    }
    return mpHicmaContext;
}

void *ExaGeoStatHardware::GetChameleonContext() {
    if (!mpChameleonContext) {
        throw std::runtime_error("Chameleon Hardware is not initialized!");
    }
    return mpChameleonContext;
}

void *ExaGeoStatHardware::GetContext(common::Computation aComputation) {
    if (aComputation == common::EXACT_DENSE || aComputation == common::DIAGONAL_APPROX) {
        return GetChameleonContext();
    }
    if (aComputation == common::TILE_LOW_RANK) {
        return GetHicmaContext();
    }
    return nullptr;
}

void * ExaGeoStatHardware::mpChameleonContext = nullptr;
void * ExaGeoStatHardware::mpHicmaContext = nullptr;