
// Copyright (c) 2017-2024 King Abdullah University of Science and Technology,
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
using namespace exageostat::results;

ExaGeoStatHardware::ExaGeoStatHardware(const Computation &aComputation, const int &aCoreNumber, const int &aGpuNumber,
                                       const int &aP, const int &aQ) {
    InitHardware(aComputation, aCoreNumber, aGpuNumber, aP, aQ);
}

// Constructor for R
ExaGeoStatHardware::ExaGeoStatHardware(const std::string &aComputation, const int &aCoreNumber, const int &aGpuNumber,
                                       const int &aP, const int &aQ) {
    InitHardware(GetInputComputation(aComputation), aCoreNumber, aGpuNumber, aP, aQ);
}

void ExaGeoStatHardware::InitHardware(const Computation &aComputation, const int &aCoreNumber, const int &aGpuNumber,
                                      const int &aP, const int &aQ) {

    SetPGrid(aP);
    SetQGrid(aQ);
    int tag_width = 31, tag_sep = 40;

    // Init hardware using Chameleon
    if (!mpChameleonContext) {
#ifdef USE_MPI
        // Due to a bug in Chameleon if CHAMELEON_user_tag_size is called twice with MPI an error happens due to a static variable in chameleon that doesn't change.
        if(!mIsMPIInit){
            CHAMELEON_user_tag_size(tag_width, tag_sep);
            mIsMPIInit = true;
        }
#else
        CHAMELEON_user_tag_size(tag_width, tag_sep);
#endif

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
    LOGGER("** Initialize ExaGeoStat hardware **")
}

void ExaGeoStatHardware::FinalizeHardware() {
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

ExaGeoStatHardware::~ExaGeoStatHardware() {

    Results::GetInstance()->PrintEndSummary();
    FinalizeHardware();
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

int ExaGeoStatHardware::GetPGrid() {
    return mPGrid;
}

int ExaGeoStatHardware::GetQGrid() {
    return mQGrid;
}

void ExaGeoStatHardware::SetPGrid(int aP) {
    mPGrid = aP;
}

void ExaGeoStatHardware::SetQGrid(int aQ) {
    mQGrid = aQ;
}

void *ExaGeoStatHardware::mpChameleonContext = nullptr;
void *ExaGeoStatHardware::mpHicmaContext = nullptr;
int ExaGeoStatHardware::mPGrid = 1;
int ExaGeoStatHardware::mQGrid = 1;
bool ExaGeoStatHardware::mIsMPIInit = false;
