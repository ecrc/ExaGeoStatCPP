
// Copyright (c) 2017-2023 King Abdullah University of Science and Technology,
// All rights reserved.
// ExaGeoStat is a software package, provided by King Abdullah University of Science and Technology (KAUST).

/**
 * @file ExaGeoStatHardware.cpp
 * @brief Contains the implementation of the ExaGeoStatHardware class.
 * @version 1.0.0
 * @author Mahmoud ElKarargy
 * @author Sameh Abdulah
 * @date 2023-08-07
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
    if (!this->mpChameleonContext) {
        CHAMELEON_user_tag_size(tag_width, tag_sep);
        CHAMELEON_Init(aCoreNumber, aGpuNumber)
        this->mpChameleonContext = chameleon_context_self();
    }

    // Init hardware using Hicma
    if (aComputation == common::TILE_LOW_RANK) {
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
    helpers::CommunicatorMPI::GetInstance()->SetHardwareInitialization();
}

ExaGeoStatHardware::~ExaGeoStatHardware() {
    // finalize hardware using Hicma
    // finalize hardware using Chameleon
    if (!this->mpChameleonContext) {
        std::cerr << "No initialized context of Chameleon, Please initialize a hardware first" << std::endl;
        exit(1);
    } else {
        CHAMELEON_Finalize()
        this->mpChameleonContext = nullptr;
    }
    if (this->mComputation == common::TILE_LOW_RANK) {
#ifdef USE_HICMA
        if (!this->mpHicmaContext) {
            std::cout
                    << "No initialized context of HiCMA, Please use 'ExaGeoStatHardware::ExaGeoStatHardware(aComputation, CoreNumber, aGpuNumber);'"
                    << std::endl;
        } else {
            HICMA_Finalize();
            this->mpHicmaContext = nullptr;
        }
#endif
    }
    helpers::CommunicatorMPI::GetInstance()->RemoveHardwareInitialization();
    results::Results::GetInstance()->PrintEndSummary();
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

void *ExaGeoStatHardware::GetContext(common::Computation aComputation) const {
    if (aComputation == common::EXACT_DENSE || aComputation == common::DIAGONAL_APPROX) {
        return GetChameleonContext();
    }
    if (aComputation == common::TILE_LOW_RANK) {
#ifdef USE_HICMA
        return GetHicmaContext();
#endif
    }
    return nullptr;
}
