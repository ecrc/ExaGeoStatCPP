
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

ExaGeoStatHardware::ExaGeoStatHardware(const exageostat::common::Computation &aComputation, const int &aCoreNumber,
                                       const int &aGpuNumber) {
    this->InitHardware(aComputation, aCoreNumber, aGpuNumber);
}

// Constructor for R
ExaGeoStatHardware::ExaGeoStatHardware(const std::string &aComputation, const int &aCoreNumber, const int &aGpuNumber) {

    this->InitHardware(GetInputComputation(aComputation), aCoreNumber, aGpuNumber);
}

void ExaGeoStatHardware::InitHardware(const exageostat::common::Computation &aComputation, const int &aCoreNumber,
                                      const int &aGpuNumber) {
    LOGGER("** Initialise ExaGeoStat hardware **")
    int tag_width = 31, tag_sep = 26;
    // Init hardware using Chameleon
    if (!this->mpChameleonContext) {
        CHAMELEON_user_tag_size(tag_width, tag_sep);
        CHAMELEON_Init(aCoreNumber, aGpuNumber)
        this->mpChameleonContext = chameleon_context_self();
    }

    // Init hardware using HiCMA
    if (aComputation == exageostat::common::TILE_LOW_RANK) {
#ifdef USE_HICMA
        if (!this->mpHicmaContext) {
            HICMA_user_tag_size(tag_width, tag_sep);
            HICMA_Init(aCoreNumber, aGpuNumber);
            this->mpHicmaContext = hicma_context_self();
        }
#else
        throw std::runtime_error("You need to enable Hicma to use TLR computation!");
#endif
    }
    exageostat::helpers::CommunicatorMPI::GetInstance()->SetHardwareInitialization();
}

ExaGeoStatHardware::~ExaGeoStatHardware() {
    // finalize hardware using HiCMA
    // finalize hardware using Chameleon
    if (!this->mpChameleonContext) {
        std::cerr << "No initialized context of Chameleon, Please initialize a hardware first" << std::endl;
        exit(1);
    } else {
        CHAMELEON_Finalize()
        this->mpChameleonContext = nullptr;
    }
#ifdef USE_HICMA
    // In case of HiCMA, It may be enabled but without initialize the hardware. aka: dense or dst computation.
    if (this->mpHicmaContext) {
        HICMA_Finalize();
        this->mpHicmaContext = nullptr;
    }
#endif
    exageostat::helpers::CommunicatorMPI::GetInstance()->RemoveHardwareInitialization();
    exageostat::results::Results::GetInstance()->PrintEndSummary();
}

#ifdef USE_HICMA

void *ExaGeoStatHardware::GetHicmaContext() const {
    if (!this->mpHicmaContext) {
        throw std::runtime_error("Hardware is not initialized!");
    }
    return this->mpHicmaContext;
}

#endif

void *ExaGeoStatHardware::GetChameleonContext() const {
    if (!this->mpChameleonContext) {
        throw std::runtime_error("Hardware is not initialized!");
    }
    return this->mpChameleonContext;
}

void *ExaGeoStatHardware::GetContext(exageostat::common::Computation aComputation) const {
    if (aComputation == exageostat::common::EXACT_DENSE || aComputation == exageostat::common::DIAGONAL_APPROX) {
        return GetChameleonContext();
    }
    if (aComputation == exageostat::common::TILE_LOW_RANK) {
#ifdef USE_HICMA
        return GetHicmaContext();
#endif
    }
    return nullptr;
}

