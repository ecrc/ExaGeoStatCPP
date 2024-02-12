
// Copyright (c) 2017-2023 King Abdullah University of Science and Technology,
// All rights reserved.
// ExaGeoStat is a software package, provided by King Abdullah University of Science and Technology (KAUST).

/**
 * @file CommunicatorMPI.cpp
 * @brief Defines the CommunicatorMPI class for MPI rank communication.
 * @version 1.0.1
 * @author Sameh Abdulah
 * @date 2023-11-10
**/

#include <helpers/CommunicatorMPI.hpp>
#include <linear-algebra-solvers/concrete/ChameleonHeaders.hpp>

using namespace exageostat::helpers;

CommunicatorMPI *CommunicatorMPI::GetInstance() {
    if (mpInstance == nullptr) {
        mpInstance = new CommunicatorMPI();
    }
    return mpInstance;
}

int CommunicatorMPI::GetRank() {
#ifdef USE_MPI
    if(!mIsHardwareInitialized){
        return 0;
    }
    else{
        return CHAMELEON_Comm_rank();
    }
#else
    return 0;
#endif
}

void CommunicatorMPI::SetHardwareInitialization() {
    mIsHardwareInitialized = true;
}
void CommunicatorMPI::RemoveHardwareInitialization() {
    mIsHardwareInitialized = false;
}

CommunicatorMPI *CommunicatorMPI::mpInstance = nullptr;
