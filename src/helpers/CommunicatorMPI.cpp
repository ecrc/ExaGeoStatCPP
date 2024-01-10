
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

bool CommunicatorMPI::GetRank() const {

#ifdef USE_MPI
    if(!this->mIsHardwareInitialized){
        return false;
    }
    else{
        if(CHAMELEON_Comm_rank() == 0){
            return true;
        }
        return false;
    }
#else
    return true;
#endif
}

void CommunicatorMPI::SetHardwareInitialization() {
    this->mIsHardwareInitialized = true;
}

CommunicatorMPI *CommunicatorMPI::mpInstance = nullptr;

void CommunicatorMPI::RemoveHardwareInitialization() {
    this->mIsHardwareInitialized = false;
}
