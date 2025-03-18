
// Copyright (c) 2017-2024 King Abdullah University of Science and Technology,
// All rights reserved.
// ExaGeoStat is a software package, provided by King Abdullah University of Science and Technology (KAUST).

/**
 * @file Configurations.cpp
 * @brief This file defines the Configurations class which stores the configuration parameters for ExaGeoStat.
 * @version 1.1.0
 * @author Mahmoud ElKarargy
 * @author Sameh Abdulah
 * @date 2024-02-04
**/

#include <configurations/Configurations.hpp>
#include <configurations/Validator.hpp>
#include <configurations/Parser.hpp>
#include <utilities/Logger.hpp>
#include <kernels/Kernel.hpp>

using namespace std;

using namespace exageostat::configurations;
using namespace exageostat::configurations::parser;
using namespace exageostat::configurations::validator;
using namespace exageostat::common;

Verbose Configurations::mVerbosity = Verbose::STANDARD_MODE;
bool Configurations::mIsThetaInit = false;
bool Configurations::mHeapAllocated = false;
bool Configurations::mFirstInit = false;

Configurations::Configurations() {

    // Set default values from the json file
    Parser::ParseJSON(DEFAULT_CONFIGURATION_PATH, this->mDictionary);
    Validator::Validate(this->mDictionary);

    // initialize thetas with empty vectors
    vector<double> theta;
    SetInitialTheta(theta);
    SetLowerBounds(theta);
    SetUpperBounds(theta);
    SetEstimatedTheta(theta);
    mIsThetaInit = false;

#if !DEFAULT_RUNTIME
    // Set default values for Hicma-Parsec params
    SetTolerance(0);
    //TODO:currently,we support real data only in parsec.In the future,we should support synthetic and real data for both runtimes
    SetIsSynthetic(false);
    SetStageZero(false);
#endif
}

void Configurations::InitializeArguments(const int &aArgC, char **apArgV, const bool &aEnableR) {

    this->mArgC = aArgC;
    this->mpArgV = apArgV;
    mHeapAllocated = aEnableR;

    // the CLI arguments overwrite the arguments of the constructor
    Parser::ParseCLI(aArgC, apArgV, this->mDictionary);
    Validator::Validate(this->mDictionary);
    ValidateConfiguration();
}

void Configurations::ValidateConfiguration() {

    // Throw Errors if any of these arguments aren't given by the user.
    if (GetProblemSize() == 0 && GetIsSynthetic()) {
        throw domain_error("You need to set the problem size, before starting");
    }

    if (GetDenseTileSize() == 0) {
        throw domain_error("You need to set the Dense tile size, before starting");
    }

    if (!GetDataPath().empty()) {
        SetIsSynthetic(false);
    }

    if (GetIsMSPE() || GetIsMLOEMMOM() || GetIsIDW()) {
        if (GetUnknownObservationsNb() <= 1) {
            throw domain_error(
                    "You need to set ZMiss number, as the number of missing values should be bigger than one");
        }
    }

    if (!GetLoggerPath().empty() && !GetLogger()) {
        throw domain_error("To enable logging, please utilize the '--log' option in order to specify a log file.");
    }

    if (GetUnknownObservationsNb() >= GetProblemSize()) {
        throw range_error("Invalid value for ZmissNumber. Please make sure it's smaller than Problem size");
    }

    if (GetComputation() == DIAGONAL_APPROX) {
        if (GetBand() == 0) {
            throw domain_error("You need to set the tile band thickness, before starting");
        }
    }

#if DEFAULT_RUNTIME
    // Throw Errors if any of these arguments aren't given by the user.
    if (GetKernelName().empty()) {
        throw domain_error("You need to set the Kernel, before starting");
    }
    if (GetMaxRank() == -1) {
        SetMaxRank(1);
    }
//#else
    if(GetMaxRank() == -1){
        SetMaxRank(GetDenseTileSize() / 2);
    }
//    if(GetTolerance() >= 0){
        SetTolerance(8);
    }
     if (GetDataPath().empty()) {
        throw domain_error("You need to set the data path, before starting");
    }
#else
    if(GetStageZero() && GetKernelName().empty()){
        throw domain_error("You need to set the Kernel for stage 0, before starting");
    }
#endif

    size_t found = GetKernelName().find("NonGaussian");
    // Check if the substring was found
    if (found != std::string::npos) {
        SetIsNonGaussian(true);
    }

    if (GetDimension() != DimensionST) {
        if (GetTimeSlot() != 1) {
#if DEFAULT_RUNTIME
            throw std::runtime_error("Time Slot can only be greater than 1 if the dimensions are set to SpaceTime.");
#endif
        }
    } else if (GetTimeSlot() < 1) {
        throw std::runtime_error("Time Slot must be at least 1 if the dimensions are set to SpaceTime.");
    }


    if (GetComputation() == TILE_LOW_RANK) {
#ifdef USE_HICMA
        if (GetLowTileSize() == 0) {
            throw domain_error("You need to set the Low tile size, before starting");
        }
#endif
    }
}

void Configurations::InitializeAllTheta() {

    if (!mIsThetaInit) {

        int parameters_number = kernels::KernelsConfigurations::GetParametersNumberKernelMap()[this->GetKernelName()];
        InitTheta(GetInitialTheta(), parameters_number);
        SetInitialTheta(GetInitialTheta());


        if (this->GetIsNonGaussian()) {
            GetInitialTheta()[GetInitialTheta().size() - 1] = 0.2;
            GetInitialTheta()[GetInitialTheta().size() - 2] = 0.2;
        }

        InitTheta(GetLowerBounds(), parameters_number);
        SetLowerBounds(GetLowerBounds());
        InitTheta(GetUpperBounds(), parameters_number);
        SetUpperBounds(GetUpperBounds());
        SetStartingTheta(GetLowerBounds());
        InitTheta(GetEstimatedTheta(), parameters_number);
        SetEstimatedTheta(GetEstimatedTheta());

        for (int i = 0; i < parameters_number; i++) {
            if (GetEstimatedTheta()[i] != -1) {
                GetLowerBounds()[i] = GetEstimatedTheta()[i];
                GetUpperBounds()[i] = GetEstimatedTheta()[i];
                GetStartingTheta()[i] = GetEstimatedTheta()[i];
            }
        }
        mIsThetaInit = true;
    }
}

void Configurations::PrintUsage() {

    LOGGER("\n\t*** Available Arguments For ExaGeoStat Configurations ***")
    LOGGER("--N=value : Problem size.")
    LOGGER("--kernel=value : Used Kernel.")
    LOGGER("--dimension=value : Used Dimension.")
    LOGGER("--p=value : Used P-Grid.")
    LOGGER("--q=value : Used P-Grid.")
    LOGGER("--time_slot=value : Time slot value for ST.")
    LOGGER("--computation=value : Used computation.")
    LOGGER("--precision=value : Used precision.")
    LOGGER("--cores=value : Used to set the number of cores.")
    LOGGER("--gpus=value : Used to set the number of GPUs.")
    LOGGER("--dts=value : Used to set the Dense Tile size.")
    LOGGER("--lts=value : Used to set the Low Tile size.")
    LOGGER("--band=value : Used to set the Tile diagonal thickness.")
    LOGGER("--Zmiss=value : Used to set number of unknown observation to be predicted.")
    LOGGER("--observations_file=PATH/TO/File : Used to pass the observations file path.")
    LOGGER("--max_rank=value : Used to the max rank value.")
    LOGGER("--initial_theta=value : Initial theta parameters for optimization.")
    LOGGER("--estimated_theta=value : Estimated kernel parameters for optimization.")
    LOGGER("--seed=value : Seed value for random number generation.")
    LOGGER("--verbose=value : Run mode whether quiet/standard/detailed.")
    LOGGER("--log_path=value : Path to log file.")
    LOGGER("--distance_metric=value : Used distance metric either eg or gcd.")
    LOGGER("--max_mle_iterations=value : Maximum number of MLE iterations.")
    LOGGER("--tolerance : MLE tolerance between two iterations.")
    LOGGER("--data_path : Used to enter the path to the real data file.")
    LOGGER("--mspe: Used to enable mean square prediction error.")
    LOGGER("--fisher: Used to enable fisher tile prediction function.")
    LOGGER("--idw: Used to IDW prediction auxiliary function.")
    LOGGER("--mloe-mmom: Used to enable MLOE MMOM.")
    LOGGER("--OOC : Used to enable Out of core technology.")
    LOGGER("--approximation_mode : Used to enable Approximation mode.")
    LOGGER("--log : Enable logging.")
    LOGGER("--accuracy : Used to set the accuracy when using tlr.")
    LOGGER("--band_dense=value : Used to set the dense band double precision, Used with PaRSEC runtime only.")
    LOGGER("--objects_number=value : Used to set the number of objects (number of viruses within a population), Used with PaRSEC runtime only.")
    LOGGER("--adaptive_decision=value : Used to set the adaptive decision of each tile's format using norm approach, if enabled, otherwise 0, Used with PaRSEC runtime only.")
    LOGGER("--add_diagonal=value : Used to add this number to diagonal elements to make the matrix positive definite in electrodynamics problem, Used with PaRSEC runtime only.")
    LOGGER("--file_time_slot=value : Used to set time slot per file, Used with PaRSEC runtime only.")
    LOGGER("--file_number=value : Used to set file number, Used with PaRSEC runtime only.")
    LOGGER("--enable-inverse : Used to enable inverse spherical harmonics transform, Used with PaRSEC runtime only.")
    LOGGER("--mpiio : Used to enable MPI IO, Used with PaRSEC runtime only.")
    LOGGER("--log-file-path: Used to set path of file where events and results are logged.")
    LOGGER("\n\n")

    exit(0);
}

Verbose Configurations::GetVerbosity() {
    return Configurations::mVerbosity;
}

void Configurations::SetVerbosity(const Verbose &aVerbose) {
    Configurations::mVerbosity = aVerbose;
}

void Configurations::InitTheta(vector<double> &aTheta, const int &size) {

    // If null, this mean user have not passed the values arguments, Make values equal -1
    if (aTheta.empty()) {
        for (int i = 0; i < size; i++) {
            aTheta.push_back(-1);
        }
    } else if (aTheta.size() < size) {
        // Also allocate new memory as maybe they are not the same size.
        for (size_t i = aTheta.size(); i < size; i++) {
            aTheta.push_back(0);
        }
    }
}

void Configurations::PrintSummary() {

#ifndef USE_R
    Verbose temp = Configurations::GetVerbosity();
    mVerbosity = STANDARD_MODE;

    if (!mFirstInit) {

        LOGGER("********************SUMMARY**********************")
#if DEFAULT_RUNTIME
        if (this->GetIsSynthetic()) {
            LOGGER("#Synthetic Data generation")
        } else {
            LOGGER("#Real Data loader")
        }
        LOGGER("#Number of Locations: " << this->GetProblemSize())
        LOGGER("#Threads per node: " << this->GetCoresNumber())
        LOGGER("#GPUs: " << this->GetGPUsNumbers())
        if (this->GetPrecision() == 1) {
            LOGGER("#Precision: Double")
        } else if (this->GetPrecision() == 0) {
            LOGGER("#Precision: Single")
        } else if (this->GetPrecision() == 2) {
            LOGGER("#Precision: Single/Double")
        }
        LOGGER("#Dense Tile Size: " << this->GetDenseTileSize())
#ifdef USE_HICMA
        LOGGER("#Low Tile Size: " << this->GetLowTileSize())
#endif
        if (this->GetComputation() == TILE_LOW_RANK) {
            LOGGER("#Computation: Tile Low Rank")
        } else if (this->GetComputation() == EXACT_DENSE) {
            LOGGER("#Computation: Exact")
        } else if (this->GetComputation() == DIAGONAL_APPROX) {
            LOGGER("#Computation: Diagonal Approx")
        }

        if (this->GetDimension() == Dimension2D) {
            LOGGER("#Dimension: 2D")
        } else if (this->GetDimension() == Dimension3D) {
            LOGGER("#Dimension: 3D")
        } else if (this->GetDimension() == DimensionST) {
            LOGGER("#Dimension: ST")
        }
        LOGGER("#Kernel: " << this->GetKernelName())
        if (this->GetDistanceMetric() == EUCLIDEAN_DISTANCE) {
            LOGGER("#Distance Metric: Euclidean distance")
        } else {
            LOGGER("#Distance Metric: Great Circle Distance")
        }
        LOGGER("#p: " << this->GetPGrid() << "\t\t #q: " << this->GetQGrid())
        if (this->GetIsOOC()) {
            LOGGER("#Out Of Core (OOC) technology is enabled")
        }
#else
            LOGGER("\t#L: " << this->GetDenseTileSize() << "\t\t\t\t#T: " << this->GetTimeSlot())
            LOGGER("\t#NB: " << this->GetDenseTileSize() << "\t\t\t\t#gpus: " << this->GetGPUsNumbers())
            LOGGER("\t#Nodes: " << this->GetCoresNumber() << "\t\t\t#Time slot per file: " << GetTimeSlotPerFile());
            LOGGER("\t#Number of files: " << this->GetFileNumber() << "\t#File per node: " << ((this->GetFileNumber()%this->GetCoresNumber())? this->GetFileNumber()/this->GetCoresNumber()+1 : this->GetFileNumber()/this->GetCoresNumber()))
#endif
        LOGGER("*************************************************")
        mFirstInit = true;
    }
    mVerbosity = temp;
#endif
}

int Configurations::CalculateZObsNumber() {
    return (this->GetProblemSize()) - this->GetUnknownObservationsNb();
}

Configurations::~Configurations() {

    if (mHeapAllocated) {
        for (size_t i = 0; i < this->mArgC; ++i) {
            delete[] this->mpArgV[i];  // Delete each string
        }
        delete[] this->mpArgV;  // Delete the array of pointers
    }
    this->mpArgV = nullptr;
    mFirstInit = false;
}

void Configurations::SetTolerance(double aTolerance) {
    mDictionary["tolerance"] = pow(10, -1 * aTolerance);
}
