
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

#include <hardware/ExaGeoStatHardware.hpp>

#if DEFAULT_RUNTIME
#ifdef USE_MPI
#include <starpu_mpi.h>
#endif

#include <linear-algebra-solvers/concrete/ChameleonHeaders.hpp>
#include <linear-algebra-solvers/concrete/HicmaHeaders.hpp>
#else
#include <runtime/parsec/ParsecHeader.h>
#endif

#include <results/Results.hpp>
#include <helpers/CommunicatorMPI.hpp>
#include <utilities/Logger.hpp>
#include <utilities/EnumStringParser.hpp>

using namespace exageostat::common;
using namespace exageostat::results;
using namespace std;

ExaGeoStatHardware::ExaGeoStatHardware(exageostat::configurations::Configurations &aConfigurations){

    // These variables are named according to HiCMA-X inputs
    const int N = aConfigurations.GetProblemSize();
    const int t = aConfigurations.GetDenseTileSize();
    const int e = aConfigurations.GetAccuracy();
    const int a = aConfigurations.GetAdaptiveDecision();
    const int g = aConfigurations.GetGPUsNumbers();
    const int c = aConfigurations.GetCoresNumber();
    const int j = aConfigurations.GetDiagonalAddition();
    const int J = aConfigurations.GetTimeSlot();
    const int K = aConfigurations.GetObjectsNumber();
    const int I = aConfigurations.GetDenseBandDP();
    const int time_slot_per_file = aConfigurations.GetTimeSlotPerFile();
    const int num_file = aConfigurations.GetFileNumber();

    int v = 0;
    if (aConfigurations.GetVerbosity() == Verbose::DETAILED_MODE){
        v = 1;
    }

    // Create a vector to store the arguments as strings
    std::vector<std::string> new_args = {
        "-g", to_string(g),
        "-NB", to_string(t),
        "-K", to_string(t),
        "-N", to_string(N),
        "-v", to_string(v),
        "-I", to_string(I),
        "-a", to_string(a),
        "-J", to_string(J),
        "-c", to_string(c),
        "-K", to_string(K),
        "-j", to_string(j)
    };

    // Convert std::vector<std::string> to char** for the new argv
    int new_argc = new_args.size();
    char **new_argv = new char*[new_argc];

    for (int i = 0; i < new_argc; ++i) {
        new_argv[i] = new char[new_args[i].length() + 1];
        strcpy(new_argv[i], new_args[i].c_str());
    }

#if !DEFAULT_RUNTIME
    int iparam[IPARAM_SIZEOF] = {0};
    double dparam[DPARAM_SIZEOF];
    char *cparam[CPARAM_SIZEOF];
    this->mpHicmaParams = make_unique<hicma_parsec_params_t>();
    this->mpParamsKernel = make_unique<starsh_params_t>();
    this->mpHicmaData = make_unique<hicma_parsec_data_t>();
    this->mpAnalysis = make_unique<hicma_parsec_matrix_analysis_t>();

    this->mpParsecContext = hicma_parsec_init(new_argc, new_argv, iparam, dparam, cparam, this->mpHicmaParams.get(), this->mpParamsKernel.get(), this->mpHicmaData.get());
    SetParsecMPIRank(this->mpHicmaParams->rank);
#endif
    exageostat::helpers::CommunicatorMPI::GetInstance()->SetHardwareInitialization();
}

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

#if DEFAULT_RUNTIME
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
#ifdef USE_MPI
            // Due to a bug in HiCMA if HICMA_user_tag_size is called twice with MPI an error happens due to a static variable in HiCMA that doesn't change.
            if(!mIsMPIInit){
                HICMA_user_tag_size(tag_width, tag_sep);
                mIsMPIInit = true;
            }
#else
            HICMA_user_tag_size(tag_width, tag_sep);
#endif
            HICMA_Init(aCoreNumber, aGpuNumber);
            mpHicmaContext = hicma_context_self();
        }
#else
        throw std::runtime_error("You need to enable HiCMA to use TLR computation!");
#endif
    }
#endif
    exageostat::helpers::CommunicatorMPI::GetInstance()->SetHardwareInitialization();
    LOGGER("** Initialize ExaGeoStat hardware **")
}

void ExaGeoStatHardware::FinalizeHardware() {

#if DEFAULT_RUNTIME
    // finalize hardware using HiCMA
    #ifdef USE_HICMA
    if (mpHicmaContext) {
        HICMA_Finalize();
        mpHicmaContext = nullptr;
    }
    #endif

    // finalize hardware using Chameleon
    if (mpChameleonContext) {
    #if defined(USE_MPI) && defined(USE_HICMA)
        // Since already HiCMA do so, then no need to remove empty cache.
        starpu_mpi_cache_set(0);
    #endif
        CHAMELEON_Finalize()
        mpChameleonContext = nullptr;
    }
#else
    if (mpParsecContext) {

        int iparam[IPARAM_SIZEOF] = {0};
        double dparam[DPARAM_SIZEOF];
        char *cparam[CPARAM_SIZEOF];

        hicma_parsec_fini((parsec_context_t *) mpParsecContext, 0, NULL, iparam, dparam, cparam, this->mpHicmaParams.get(), this->mpParamsKernel.get(), this->mpHicmaData.get(), this->mpAnalysis.get());
        mpParsecContext = nullptr;
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

void *ExaGeoStatHardware::GetParsecContext() {
    if (!mpParsecContext) {
        throw std::runtime_error("PaRSEC Hardware is not initialized!");
    }
    return mpParsecContext;
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

void ExaGeoStatHardware::SetParsecMPIRank(int aRank){
    mParsecMPIRank = aRank;
}

int ExaGeoStatHardware::GetParsecMPIRank() {
    return mParsecMPIRank;
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

#if !DEFAULT_RUNTIME
hicma_parsec_params_t* ExaGeoStatHardware::GetHicmaParams() {
    return mpHicmaParams.get();
}

starsh_params_t* ExaGeoStatHardware::GetParamsKernel() {
    return mpParamsKernel.get();
}

hicma_parsec_data_t* ExaGeoStatHardware::GetHicmaData() {
    return mpHicmaData.get();
}

hicma_parsec_matrix_analysis_t* ExaGeoStatHardware::GetAnalysis() {
    return mpAnalysis.get();
}
#endif

void *ExaGeoStatHardware::mpChameleonContext = nullptr;
void *ExaGeoStatHardware::mpHicmaContext = nullptr;
void *ExaGeoStatHardware::mpParsecContext = nullptr;
int ExaGeoStatHardware::mParsecMPIRank = 0;
int ExaGeoStatHardware::mPGrid = 1;
int ExaGeoStatHardware::mQGrid = 1;
bool ExaGeoStatHardware::mIsMPIInit = false;
#if !DEFAULT_RUNTIME
unique_ptr<hicma_parsec_params_t> ExaGeoStatHardware::mpHicmaParams = nullptr;
unique_ptr<starsh_params_t> ExaGeoStatHardware::mpParamsKernel = nullptr;
unique_ptr<hicma_parsec_data_t> ExaGeoStatHardware::mpHicmaData = nullptr;
unique_ptr<hicma_parsec_matrix_analysis_t> ExaGeoStatHardware::mpAnalysis = nullptr;
#endif