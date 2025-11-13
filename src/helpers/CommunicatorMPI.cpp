
// Copyright (c) 2017-2024 King Abdullah University of Science and Technology,
// All rights reserved.
// ExaGeoStat is a software package, provided by King Abdullah University of Science and Technology (KAUST).

/**
 * @file CommunicatorMPI.cpp
 * @brief Defines the CommunicatorMPI class for MPI rank communication.
 * @version 1.1.0
 * @author Sameh Abdulah
 * @date 2023-11-10
**/
#include <iostream>
#include <helpers/CommunicatorMPI.hpp>
#if DEFAULT_RUNTIME
#include <linear-algebra-solvers/concrete/ChameleonHeaders.hpp>
#endif
#include <hardware/ExaGeoStatHardware.hpp>

using namespace exageostat::helpers;

CommunicatorMPI *CommunicatorMPI::GetInstance() {
    if (mpInstance == nullptr) {
        mpInstance = new CommunicatorMPI();
    }
    return mpInstance;
}

int CommunicatorMPI::GetRank() const {
#ifdef USE_MPI
    if (!mIsHardwareInitialized) {
        return 0;
    }
    #if DEFAULT_RUNTIME
    else {
        return CHAMELEON_Comm_rank();
    }
    #else
    else {
        return ExaGeoStatHardware::GetParsecMPIRank();
    }
    #endif
#endif
    return 0;
}

void CommunicatorMPI::SetHardwareInitialization() {
    mIsHardwareInitialized = true;
}

void CommunicatorMPI::RemoveHardwareInitialization() {
    mIsHardwareInitialized = false;
}

CommunicatorMPI *CommunicatorMPI::mpInstance = nullptr;
