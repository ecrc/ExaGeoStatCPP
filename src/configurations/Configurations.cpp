
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
    SetSeed(static_cast<unsigned int>(time(0)));
    SetLogger(false);
    SetUnknownObservationsNb(0);
    SetApproximationMode(1);
    SetActualObservationsFilePath("");
    SetRecoveryFile("");
    SetPrecision(DOUBLE);
    SetIsMSPE(false);
    SetIsFisher(false);
    SetIsIDW(false);
    SetIsMLOEMMOM(false);
    SetDataPath("");
    SetDistanceMetric(EUCLIDEAN_DISTANCE);
    SetAccuracy(0);
    SetIsNonGaussian(false);
    mIsThetaInit = false;

#if !DEFAULT_RUNTIME
    // Set default values for Hicma-Parsec params
    SetTolerance(0);
    //TODO:currently,we support real data only in parsec.In the future,we should support synthetic and real data for both runtimes
    SetIsSynthetic(false);
    SetMeanTrendRemoval(false);
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

    if (GetMeanTrendRemoval()) {
        if (GetResultsPath().empty()) {
            throw domain_error("You need to set the results path (--resultspath) before starting");
        }

        if (GetLatitudeBand() < 0) {
            throw domain_error("You need to set the latitude band (--lat) for Mean Trend Removal");
        }
        
        if (GetLongitudeCount() <= 0) {
            throw domain_error("You need to set the longitude count (--lon) for Mean Trend Removal");
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
    if (mDictionary.find("tolerance") == mDictionary.end()) {
        SetTolerance(8);
    }
     if (GetDataPath().empty()) {
        throw domain_error("You need to set the data path, before starting");
    }
#else
    if(GetMeanTrendRemoval() && GetKernelName().empty()){
        throw domain_error("You need to set the Kernel for Mean Trend Removal, before starting");
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
    LOGGER("--dimension=value : Used Dimension (2D, 3D, ST).")
    LOGGER("--p=value : Used P-Grid.")
    LOGGER("--q=value : Used Q-Grid.")
    LOGGER("--time_slot=value : Time slot value for ST.")
    LOGGER("--computation=value : Used computation (exact, tlr, diagonal_approx).")
    LOGGER("--precision=value : Used precision (single, double, mixed).")
    LOGGER("--cores=value : Used to set the number of cores.")
    LOGGER("--gpus=value : Used to set the number of GPUs.")
    LOGGER("--dts=value : Used to set the Dense Tile size.")
    LOGGER("--lts=value : Used to set the Low Tile size.")
    LOGGER("--band=value : Used to set the Tile diagonal thickness.")
    LOGGER("--max_rank=value : Used to set the max rank value.")
    LOGGER("--hnb=value : Used to set HNB value.")
    LOGGER("--gen_max_rank=value : Used to set generation max rank.")
    LOGGER("--comp_max_rank=value : Used to set computation max rank.")
    LOGGER("--auto_band=value : Used to set auto band.")
    LOGGER("--band_dense_sp=value : Used to set band dense single precision.")
    LOGGER("--band_low_rank_dp=value : Used to set band low rank double precision.")
    LOGGER("--observations_file=PATH/TO/File : Used to pass the observations file path.")
    LOGGER("--seed=value : Seed value for random number generation.")
    LOGGER("--verbose=value : Run mode (0=quiet, 1=standard, 2=detailed).")
    LOGGER("--log_path=value : Path to log file.")
    LOGGER("--log : Enable logging.")
    LOGGER("--initial_theta=value1,value2,... : Initial theta parameters for optimization.")
    LOGGER("--estimated_theta=value1,value2,... : Estimated kernel parameters for optimization.")
    LOGGER("--lb=value1,value2,... : Lower bounds for optimization.")
    LOGGER("--ub=value1,value2,... : Upper bounds for optimization.")
    LOGGER("--starting_theta=value1,value2,... : Starting theta parameters.")
    LOGGER("--is_non_gaussian : Enable non-Gaussian mode.")
    LOGGER("--OOC : Used to enable Out of Core (OOC) technology.")
    LOGGER("--approximation_mode=value : Used to enable Approximation mode.")
    LOGGER("--accuracy=value : Used to set the accuracy when using TLR.")
    LOGGER("\t=== DATA GENERATION ARGUMENTS ===")
    LOGGER("--data_path=PATH : Used to enter the path to the real data file.")
    LOGGER("--is_synthetic : Use synthetic data generation.")
    LOGGER("--resultspath=PATH : Used to set the output directory path for results.")
    LOGGER("\t=== DATA MODELING ARGUMENTS ===")
    LOGGER("--distance_metric=value : Used distance metric (eg=Euclidean, gcd=Great Circle Distance).")
    LOGGER("--max_mle_iterations=value : Maximum number of MLE iterations.")
    LOGGER("--tolerance=value : MLE tolerance between two iterations.")
    LOGGER("--recovery_file=PATH : Path to recovery file.")
    LOGGER("\t=== DATA PREDICTION ARGUMENTS ===")
    LOGGER("--Zmiss=value : Used to set number of unknown observations to be predicted.")
    LOGGER("--observation_number=value : Used to set the number of observations.")
    LOGGER("--mspe : Used to enable Mean Square Prediction Error.")
    LOGGER("--fisher : Used to enable Fisher tile prediction function.")
    LOGGER("--idw : Used to enable IDW prediction auxiliary function.")
    LOGGER("--mloe-mmom : Used to enable MLOE MMOM.")
    LOGGER("\t=== PARSEC RUNTIME SPECIFIC ARGUMENTS ===")
    LOGGER("--band_dense=value : Used to set the dense band double precision (PaRSEC only).")
    LOGGER("--band_dense_dp=value : Used to set dense band double precision (PaRSEC only).")
    LOGGER("--band_dense_hp=value : Used to set dense band high precision (PaRSEC only).")
    LOGGER("--objects_number=value : Used to set the number of objects (PaRSEC only).")
    LOGGER("--adaptive_decision=value : Used to set adaptive decision for tile format (PaRSEC only).")
    LOGGER("--add_diagonal=value : Add value to diagonal elements (PaRSEC only).")
    LOGGER("--file_time_slot=value : Used to set time slot per file (PaRSEC only).")
    LOGGER("--file_number=value : Used to set file number (PaRSEC only).")
    LOGGER("--enable-inverse : Enable inverse spherical harmonics transform (PaRSEC only).")
    LOGGER("--mpiio : Enable MPI IO (PaRSEC only).")
    LOGGER("--log-file-path=PATH : Path to file where events and results are logged (PaRSEC only).")
    LOGGER("\t=== MEAN TREND REMOVAL / CLIMATE EMULATOR ARGUMENTS ===")
    LOGGER("--mean_trend_removal : Enable Mean Trend Removal.")
    LOGGER("--is_climate_emulator : Enable Climate Emulator mode.")
    LOGGER("--forcing_data_path=PATH : Path to forcing data file.")
    LOGGER("--netcdf_data_path=PATH : Path to NetCDF data file.")
    LOGGER("--start-year=value : Starting year for NetCDF data processing.")
    LOGGER("--end-year=value : Ending year for NetCDF data processing.")
    LOGGER("--lat=value : Latitude band index for climate data processing (required for MeanTrendRemoval).")
    LOGGER("--lon=value : Longitude count for climate data processing (required for MeanTrendRemoval).")
    LOGGER("\n")

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
            LOGGER("#L: " << this->GetDenseTileSize())
            LOGGER("#T: " << this->GetTimeSlot())
            LOGGER("#NB: " << this->GetDenseTileSize())
            LOGGER("#gpus: " << this->GetGPUsNumbers())
            LOGGER("#Nodes: " << this->GetCoresNumber())
            LOGGER("#Time slot per file: " << GetTimeSlotPerFile())
            LOGGER("#Number of files: " << this->GetFileNumber())
            LOGGER("#File per node: " << ((this->GetFileNumber()%this->GetCoresNumber())? this->GetFileNumber()/this->GetCoresNumber()+1 : this->GetFileNumber()/this->GetCoresNumber()))
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
