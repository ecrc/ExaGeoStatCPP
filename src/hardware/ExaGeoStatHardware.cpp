
// Copyright (c) 2017-2023 King Abdullah University of Science and Technology,
// All rights reserved.
// ExaGeoStat is a software package, provided by King Abdullah University of Science and Technology (KAUST).

/**
 * @file ExaGeoStatHardware.cpp
 * @brief Contains the implementation of the ExaGeoStatHardware class.
 * @version 1.1.0
 * @author Mahmoud ElKarargy
 * @author Sameh Abdulah
 * @date 2024-02-04
**/

#include <linear-algebra-solvers/concrete/ChameleonHeaders.hpp>
#include <linear-algebra-solvers/concrete/HicmaHeaders.hpp>
#include <hardware/ExaGeoStatHardware.hpp>
#include <results/Results.hpp>
#include <helpers/CommunicatorMPI.hpp>
#include <utilities/Logger.hpp>
#include <utilities/EnumStringParser.hpp>

using namespace exageostat::common;

ExaGeoStatHardware::ExaGeoStatHardware(const Computation &aComputation, const int &aCoreNumber,
                                       const int &aGpuNumber) {
    InitHardware(aComputation, aCoreNumber, aGpuNumber);
}

// Constructor for R
ExaGeoStatHardware::ExaGeoStatHardware(const std::string &aComputation, const int &aCoreNumber, const int &aGpuNumber) {

    InitHardware(GetInputComputation(aComputation), aCoreNumber, aGpuNumber);
}

void ExaGeoStatHardware::InitHardware(const Computation &aComputation, const int &aCoreNumber,
                                      const int &aGpuNumber) {
    LOGGER("** Initialize ExaGeoStat hardware **")
    int tag_width = 31, tag_sep = 26;

    // Init hardware using Chameleon
    if (!mpChameleonContext) {
        CHAMELEON_user_tag_size(tag_width, tag_sep);
        CHAMELEON_Init(aCoreNumber, aGpuNumber)
        mpChameleonContext = chameleon_context_self();
    }

    // Init hardware using HiCMA
    if (aComputation == TILE_LOW_RANK) {
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
    exageostat::helpers::CommunicatorMPI::GetInstance()->SetHardwareInitialization();
}

void ExaGeoStatHardware::FinalizeHardware(){
    // finalize hardware using Chameleon
    if (mpChameleonContext) {
        CHAMELEON_Finalize()
        mpChameleonContext = nullptr;
    }
    // finalize hardware using HiCMA
#ifdef USE_HICMA
    if (mpHicmaContext) {
        HICMA_Finalize();
        mpHicmaContext = nullptr;
    }
#endif
}

ExaGeoStatHardware::~ExaGeoStatHardware() {

    exageostat::results::Results::GetInstance()->PrintEndSummary();
    // finalize hardware using Chameleon
    if (mpChameleonContext) {
        CHAMELEON_Finalize()
        mpChameleonContext = nullptr;
    }
    // finalize hardware using HiCMA
#ifdef USE_HICMA
        if (mpHicmaContext) {
            HICMA_Finalize();
            mpHicmaContext = nullptr;
        }
#endif
    exageostat::helpers::CommunicatorMPI::GetInstance()->RemoveHardwareInitialization();
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

void *ExaGeoStatHardware::GetContext(Computation aComputation) {
    if (aComputation == EXACT_DENSE || aComputation == DIAGONAL_APPROX) {
        return GetChameleonContext();
    }
    if (aComputation == TILE_LOW_RANK) {
        return GetHicmaContext();
    }
    return nullptr;
}

void * ExaGeoStatHardware::mpChameleonContext = nullptr;
void * ExaGeoStatHardware::mpHicmaContext = nullptr;