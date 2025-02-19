
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

#include <cstring>

#include <configurations/Configurations.hpp>
#include <utilities/Logger.hpp>
#include <kernels/Kernel.hpp>

using namespace std;

using namespace exageostat::configurations;
using namespace exageostat::common;

Verbose Configurations::mVerbosity = Verbose::STANDARD_MODE;
bool Configurations::mIsThetaInit = false;
bool Configurations::mHeapAllocated = false;
bool Configurations::mFirstInit = false;

Configurations::Configurations() {

    // Set default values for arguments!
    SetComputation(EXACT_DENSE);
    SetCoresNumber(1);
    SetGPUsNumbers(0);
    SetPGrid(1);
    SetQGrid(1);
    SetMaxRank(-1);
    SetIsOOC(false);
    SetKernelName("");
    SetDimension(Dimension2D);
    SetTimeSlot(1);
    SetProblemSize(0);
    SetDenseTileSize(0);
#ifdef USE_HICMA
    SetLowTileSize(0);
#endif
    SetBand(0);
    SetLoggerPath("");
    SetIsSynthetic(true);
    vector<double> theta;
    SetInitialTheta(theta);
    SetLowerBounds(theta);
    SetUpperBounds(theta);
    SetEstimatedTheta(theta);
    SetSeed(0);
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
    SetDenseBandDP(0);
    SetObjectsNumber(0);
    SetAdaptiveDecision(0);
    SetDiagonalAddition(0);
    SetTimeSlotPerFile(1);
    SetFileNumber(1);
    SetEnableInverse(false);
    SetMPIIO(true);
    SetTolerance(0);
    // TODO: currently, we support real data only in parsec. In the future, we should support synthetic and real data for both runtimes
    SetIsSynthetic(false);
#endif
}

void Configurations::InitializeArguments(const int &aArgC, char **apArgV, const bool &aEnableR) {

    this->mArgC = aArgC;
    this->mpArgV = apArgV;
    mHeapAllocated = aEnableR;

    // Get the example name
    string example_name = apArgV[0];
    // Remove the './'
    example_name.erase(0, 2);
    string argument;
    string argument_name;
    string argument_value;
    int equal_sign_Idx;

    // Loop through the arguments
    for (int i = 1; i < aArgC; ++i) {
        argument = apArgV[i];
        equal_sign_Idx = static_cast<int>(argument.find('='));
        argument_name = argument.substr(0, equal_sign_Idx);

        // Check if argument has an equal sign.
        if (equal_sign_Idx != string::npos) {
            argument_value = argument.substr(equal_sign_Idx + 1);

            // Check the argument name and set the corresponding value
            if (argument_name == "--N" || argument_name == "--n") {
                SetProblemSize(CheckNumericalValue(argument_value));
            } else if (argument_name == "--Kernel" || argument_name == "--kernel") {
                CheckKernelValue(argument_value);
            } else if (argument_name == "--P" || argument_name == "--p") {
                SetPGrid(max(CheckNumericalValue(argument_value), GetPGrid()));
            } else if (argument_name == "--Q" || argument_name == "--q") {
                SetQGrid(max(CheckNumericalValue(argument_value), GetQGrid()));
            } else if (argument_name == "--Dimension" || argument_name == "--dimension" || argument_name == "--dim" ||
                       argument_name == "--Dim") {
                SetDimension(CheckDimensionValue(argument_value));
            } else if (argument_name == "--TimeSlot" || argument_name == "--timeslot" ||
                       argument_name == "--time_slot") {
                SetTimeSlot(CheckNumericalValue(argument_value));
            } else if (argument_name == "--Computation" || argument_name == "--computation") {
                SetComputation(CheckComputationValue(argument_value));
            } else if (argument_name == "--precision" || argument_name == "--Precision") {
                SetPrecision(CheckPrecisionValue(argument_value));
            } else if (argument_name == "--cores" || argument_name == "--coresNumber" ||
                       argument_name == "--cores_number" || argument_name == "--ncores") {
                SetCoresNumber(CheckNumericalValue(argument_value));
            } else if (argument_name == "--gpus" || argument_name == "--GPUsNumbers" ||
                       argument_name == "--gpu_number" || argument_name == "--ngpus") {
                SetGPUsNumbers(CheckNumericalValue(argument_value));
            } else if (argument_name == "--DTS" || argument_name == "--dts" || argument_name == "--Dts") {
                SetDenseTileSize(CheckNumericalValue(argument_value));
            } else if (argument_name == "--LTS" || argument_name == "--lts" || argument_name == "--Lts") {
                SetLowTileSize(CheckNumericalValue(argument_value));
            } else if (argument_name == "--maxRank" || argument_name == "--maxrank" || argument_name == "--max_rank") {
                SetMaxRank(CheckNumericalValue(argument_value));
            } else if (argument_name == "--initial_theta" || argument_name == "--itheta" ||
                       argument_name == "--iTheta") {
                vector<double> theta = ParseTheta(argument_value);
                SetInitialTheta(theta);
            } else if (argument_name == "--lb" || argument_name == "--olb" || argument_name == "--lower_bounds") {
                vector<double> theta = ParseTheta(argument_value);
                SetLowerBounds(theta);
                SetStartingTheta(theta);
            } else if (argument_name == "--ub" || argument_name == "--oub" || argument_name == "--upper_bounds") {
                vector<double> theta = ParseTheta(argument_value);
                SetUpperBounds(theta);
            } else if (argument_name == "--estimated_theta" || argument_name == "--etheta" ||
                       argument_name == "--eTheta") {
                vector<double> theta = ParseTheta(argument_value);
                SetEstimatedTheta(theta);
            } else if (argument_name == "--ObservationsFile" || argument_name == "--observationsfile" ||
                       argument_name == "--observations_file") {
                SetActualObservationsFilePath(argument_value);
            } else if (argument_name == "--Seed" || argument_name == "--seed") {
                SetSeed(CheckNumericalValue(argument_value));
            } else if (argument_name == "--verbose" || argument_name == "--Verbose") {
                ParseVerbose(argument_value);
            } else if (argument_name == "--distance_metric" || argument_name == "--distanceMetric") {
                ParseDistanceMetric(argument_value);
            } else if (argument_name == "--logpath" || argument_name == "--log_path" || argument_name == "--logPath") {
                SetLoggerPath(argument_value);
            } else if(argument_name == "--band_dense" || argument_name == "--bandDense"){
                SetDenseBandDP(CheckNumericalValue(argument_value));
            } else if(argument_name == "--numobj" || argument_name == "--obj-num" || argument_name == "--objnum"){
                SetObjectsNumber(CheckNumericalValue(argument_value));
            } else if(argument_name == "--adaptive_decision"){
                SetAdaptiveDecision(CheckNumericalValue(argument_value));
            } else if(argument_name == "--add_diag" || argument_name == "--add-diag" || argument_name == "--add-diagonal"
            || argument_name == "--add-diag-elements"){
                SetDiagonalAddition(CheckNumericalValue(argument_value));
            } else if(argument_name == "--file_time_slot" || argument_name == "--fileTimeSlot"){
                SetTimeSlotPerFile(CheckNumericalValue(argument_value));
            } else if(argument_name == "--numFiles" || argument_name == "--file-num"){
                SetFileNumber(CheckNumericalValue(argument_value));
            } else {
                if (!(argument_name == "--ZmissNumber" || argument_name == "--Zmiss" ||
                      argument_name == "--ZMiss" || argument_name == "--predict" || argument_name == "--Predict" ||
                      argument_name == "--iterations" ||
                      argument_name == "--Iterations" || argument_name == "--max_mle_iterations" ||
                      argument_name == "--maxMleIterations" || argument_name == "--opt_iters" ||
                      argument_name == "--tolerance" || argument_name == "--opt_tol" ||
                      argument_name == "--distanceMetric" || argument_name == "--distance_metric" ||
                      argument_name == "--log_file_name" || argument_name == "--logFileName" ||
                      argument_name == "--Band" || argument_name == "--band" ||
                      argument_name == "--DataPath" || argument_name == "--dataPath" ||
                      argument_name == "--data_path" ||
                      argument_name == "--acc" || argument_name == "--Acc")) {
                    LOGGER("!! " << argument_name << " !!")
                    throw invalid_argument(
                            "This argument is undefined, Please use --help to print all available arguments");
                }
            }
        } else {
            if (argument_name == "--help") {
                PrintUsage();
            }
            if (argument_name == "--OOC" || argument_name == "--ooc") {
                SetIsOOC(true);
            } else if (argument_name == "--ApproximationMode" || argument_name == "--approximationmode" ||
                       argument_name == "--approximation_mode") {
                SetApproximationMode(true);
            } else if (argument_name == "--log" || argument_name == "--Log") {
                SetLogger(true);
            } else if(argument_name == "--enable-inverse" || argument_name == "--enable_inverse"){
                SetEnableInverse(true);
            } else if(argument_name == "--mpiio"){
                SetMPIIO(true);
            } else {
                if (!(argument_name == "--mspe" || argument_name == "--MSPE" ||
                      argument_name == "--idw" || argument_name == "--IDW" ||
                      argument_name == "--mloe-mmom" || argument_name == "--mloe_mmom" ||
                      argument_name == "--fisher" || argument_name == "--Fisher")) {
                    LOGGER("!! " << argument_name << " !!")
                    throw invalid_argument(
                            "This argument is undefined, Please use --help to print all available arguments");
                }
            }
        }
    }

    // Throw Errors if any of these arguments aren't given by the user.
    if (GetProblemSize() == 0 && GetIsSynthetic()) {
        throw domain_error("You need to set the problem size, before starting");
    }

    if (GetDenseTileSize() == 0) {
        throw domain_error("You need to set the Dense tile size, before starting");
    }

#if DEFAULT_RUNTIME
    // Throw Errors if any of these arguments aren't given by the user.
    if (GetKernelName().empty()) {
        throw domain_error("You need to set the Kernel, before starting");
    }
    if(GetMaxRank() == -1){
        SetMaxRank(1);
    }
#else
    if(GetMaxRank() == -1){
        SetMaxRank(GetDenseTileSize() / 2);
    }
    if(GetTolerance() >= 0){
        SetTolerance(8);
    }
#endif

    size_t found = GetKernelName().find("NonGaussian");
    // Check if the substring was found
    if (found != std::string::npos) {
        SetIsNonGaussian(true);
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

void Configurations::InitializeDataGenerationArguments() {

    this->InitializeAllTheta();
    string argument;
    string argument_name;
    string argument_value;
    int equal_sign_Idx;

    // Loop through the arguments that are specific for data generation.
    for (int i = 1; i < this->mArgC; ++i) {
        argument = this->mpArgV[i];
        equal_sign_Idx = static_cast<int>(argument.find('='));
        argument_name = argument.substr(0, equal_sign_Idx);

        // Check if argument has an equal sign.
        if (equal_sign_Idx != string::npos) {
            argument_value = argument.substr(equal_sign_Idx + 1);

            // Check the argument name and set the corresponding value
            if (argument_name == "--DataPath" || argument_name == "--dataPath" ||
                argument_name == "--data_path") {
                SetDataPath(argument_value);
                SetIsSynthetic(false);
            }
        }
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

#if !DEFAULT_RUNTIME
    if (GetDataPath().empty()) {
        throw domain_error("You need to set the data path, before starting");
    }
#endif
}

void Configurations::InitializeDataModelingArguments() {

    this->InitializeAllTheta();
    string argument;
    string argument_name;
    string argument_value;
    int equal_sign_Idx;

    // Loop through the arguments that are specific for data modeling.
    for (int i = 1; i < this->mArgC; ++i) {
        argument = this->mpArgV[i];
        equal_sign_Idx = static_cast<int>(argument.find('='));
        argument_name = argument.substr(0, equal_sign_Idx);

        // Check if argument has an equal sign.
        if (equal_sign_Idx != string::npos) {
            argument_value = argument.substr(equal_sign_Idx + 1);

            // Check the argument name and set the corresponding value
            if (argument_name == "--max_mle_iterations" || argument_name == "--maxMleIterations" ||
                argument_name == "--opt_iters") {
                SetMaxMleIterations(CheckNumericalValue(argument_value));
            } else if (argument_name == "--tolerance" || argument_name == "--opt_tol") {
                SetTolerance(CheckNumericalValue(argument_value));
            } else if (argument_name == "--Band" || argument_name == "--band") {
                SetBand(CheckNumericalValue(argument_value));
            } else if (argument_name == "--acc" || argument_name == "--Acc") {
                SetAccuracy(CheckNumericalValue(argument_value));
            } else if (argument_name == "--log_file_name" || argument_name == "--logFileName") {
                if (!GetLogger()) {
                    throw domain_error(
                            "To enable logging, please utilize the '--log' option in order to specify a log file.");
                }
                SetFileLogName(argument_value);
            }
        }
    }
    if (GetComputation() == DIAGONAL_APPROX) {
        if (GetBand() == 0) {
            throw domain_error("You need to set the tile band thickness, before starting");
        }
    }
    if (GetComputation() == TILE_LOW_RANK) {
#ifdef USE_HICMA
        if (GetLowTileSize() == 0) {
            throw domain_error("You need to set the Low tile size, before starting");
        }
#endif
    }
}

void Configurations::InitializeDataPredictionArguments() {

    this->InitializeAllTheta();
    string argument;
    string argument_name;
    string argument_value;
    int equal_sign_Idx;

    for (int i = 1; i < this->mArgC; ++i) {
        argument = this->mpArgV[i];
        equal_sign_Idx = static_cast<int>(argument.find('='));
        argument_name = argument.substr(0, equal_sign_Idx);
        if (equal_sign_Idx != string::npos) {
            argument_value = argument.substr(equal_sign_Idx + 1);
            if (argument_name == "--ZmissNumber" || argument_name == "--Zmiss" || argument_name == "--ZMiss" ||
                argument_name == "--predict" || argument_name == "--Predict") {
                SetUnknownObservationsNb(CheckUnknownObservationsValue(argument_value));
            }
        }
    }

    // Loop through the arguments that are specific for Prediction.
    for (int i = 1; i < this->mArgC; ++i) {
        argument = this->mpArgV[i];
        equal_sign_Idx = static_cast<int>(argument.find('='));
        argument_name = argument.substr(0, equal_sign_Idx);

        if (argument_name == "--mspe" || argument_name == "--MSPE") {
            if (GetUnknownObservationsNb() <= 1) {
                throw domain_error(
                        "You need to set ZMiss number, as the number of missing values should be bigger than one");
            }
            SetIsMSPE(true);
        } else if (argument_name == "--idw" || argument_name == "--IDW") {
            if (GetUnknownObservationsNb() <= 1) {
                throw domain_error(
                        "You need to set ZMiss number, as the number of missing values should be bigger than one");
            }
            SetIsIDW(true);
        } else if (argument_name == "--mloe-mmom" || argument_name == "--MLOE_MMOM" || argument_name == "--mloe_mmom") {
            if (GetUnknownObservationsNb() <= 1) {
                throw domain_error(
                        "You need to set ZMiss number, as the number of missing values should be bigger than one");
            }
            SetIsMLOEMMOM(true);
        } else if (argument_name == "--Fisher" || argument_name == "--fisher") {
            //Fisher can be performed without zmiss.
            SetIsFisher(true);
        }
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
    LOGGER("--olb=value : Lower bounds for optimization.")
    LOGGER("--oub=value : Upper bounds for optimization.")
    LOGGER("--itheta=value : Initial theta parameters for optimization.")
    LOGGER("--etheta=value : Estimated kernel parameters for optimization.")
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
    LOGGER("--acc : Used to set the accuracy when using tlr.")
    LOGGER("--band_dense=value : Used to set the dense band double precision, Used with PaRSEC runtime only.")
    LOGGER("--numobj=value : Used to set the number of objects (number of viruses within a population), Used with PaRSEC runtime only.")
    LOGGER("--adaptive_decision=value : Used to set the adaptive decision of each tile's format using norm approach, if enabled, otherwise 0, Used with PaRSEC runtime only.")
    LOGGER("--add_diag=value : Used to add this number to diagonal elements to make the matrix positive definite in electrodynamics problem, Used with PaRSEC runtime only.")
    LOGGER("--file_time_slot=value : Used to set time slot per file, Used with PaRSEC runtime only.")
    LOGGER("--numFiles=value : Used to set file number, Used with PaRSEC runtime only.")
    LOGGER("--enable-inverse : Used to enable inverse spherical harmonics transform, Used with PaRSEC runtime only.")
    LOGGER("--mpiio : Used to enable MPI IO, Used with PaRSEC runtime only.")
    LOGGER("\n\n")

    exit(0);
}

Verbose Configurations::GetVerbosity() {
    return Configurations::mVerbosity;
}

void Configurations::SetVerbosity(const Verbose &aVerbose) {
    Configurations::mVerbosity = aVerbose;
}

int Configurations::CheckNumericalValue(const string &aValue) {

    int numericalValue;
    try {
        numericalValue = stoi(aValue);
    }
    catch (...) {
        throw range_error("Invalid value. Please use Numerical values only.");
    }

    if (numericalValue < 0) {
        throw range_error("Invalid value. Please use positive values");
    }
    return numericalValue;
}

Computation Configurations::CheckComputationValue(const std::string &aValue) {

    if (aValue != "exact" and aValue != "Exact" and aValue != "Dense" and aValue != "dense" and
        aValue != "diag_approx" and aValue != "diagonal_approx" and aValue != "lr_approx" and aValue != "tlr" and
        aValue != "TLR") {
        throw range_error("Invalid value for Computation. Please use Exact, diagonal_approx or TLR.");
    }
    if (aValue == "exact" or aValue == "Exact" or aValue == "Dense" or aValue == "dense") {
        return EXACT_DENSE;
    } else if (aValue == "diag_approx" or aValue == "diagonal_approx") {
        return DIAGONAL_APPROX;
    }
    return TILE_LOW_RANK;
}

Precision Configurations::CheckPrecisionValue(const std::string &aValue) {

    if (aValue != "single" and aValue != "Single" and aValue != "double" and aValue != "Double" and aValue != "mix" and
        aValue != "Mix" and aValue != "Mixed" and aValue != "mixed") {
        throw range_error("Invalid value for Computation. Please use Single, Double or Mixed.");
    }
    if (aValue == "single" or aValue == "Single") {
        return SINGLE;
    } else if (aValue == "double" or aValue == "Double") {
        return DOUBLE;
    }
    return MIXED;
}

void Configurations::ParseVerbose(const std::string &aVerbosity) {
    if (aVerbosity == "quiet" || aVerbosity == "Quiet") {
        mVerbosity = Verbose::QUIET_MODE;
    } else if (aVerbosity == "standard" || aVerbosity == "Standard") {
        mVerbosity = Verbose::STANDARD_MODE;
    } else if (aVerbosity == "detailed" || aVerbosity == "Detailed" || aVerbosity == "detail") {
        mVerbosity = Verbose::DETAILED_MODE;
    } else {
        LOGGER("Error: " << aVerbosity << " is not valid ")
        throw range_error("Invalid value. Please use verbose or standard values only.");
    }
}


void Configurations::CheckKernelValue(const string &aKernel) {

    // Check if the kernel name exists in the availableKernels set.
    if (availableKernels.count(aKernel) <= 0) {
        throw range_error("Invalid value for Kernel. Please check manual.");
    }
    // Check if the string is already in CamelCase format
    if (IsCamelCase(aKernel)) {
        this->SetKernelName(aKernel);
        return;
    }
    string str = aKernel;
    // Replace underscores with spaces and split the string into words
    std::replace(str.begin(), str.end(), '_', ' ');
    std::istringstream iss(str);
    std::string word, result;
    while (iss >> word) {
        // Capitalize the first letter of each word and append it to the result
        word[0] = static_cast<char>(toupper(word[0]));
        result += word;
    }
    this->SetKernelName(result);
}

bool Configurations::IsCamelCase(const std::string &aString) {
    // If the string contains an underscore, it is not in CamelCase format
    if (aString.find('_') != std::string::npos) {
        return false;
    }
    // If the string starts with a lowercase letter, it is not in CamelCase format
    if (islower(aString[0])) {
        return false;
    }
    // If none of the above conditions hold, the string is in CamelCase format
    return true;
}

vector<double> Configurations::ParseTheta(const std::string &aInputValues) {
    // Count the number of values in the string
    int num_values = 1;
    for (char aInputValue: aInputValues) {
        if (aInputValue == ':') {
            num_values++;
        }
    }
    // Allocate memory for the array of doubles
    vector<double> theta;

    // Split the string into tokens using strtok()
    const char *delim = ":";
    char *token = strtok((char *) aInputValues.c_str(), delim);
    int i = 0;
    while (token != nullptr) {
        // Check if the token is a valid double or "?"
        if (!strcmp(token, "?")) {
            theta.push_back(-1);
        } else {
            try {
                theta.push_back(stod(token));
            }
            catch (...) {
                LOGGER("Error: " << token << " is not a valid double or '?' ")
                throw range_error("Invalid value. Please use Numerical values only.");
            }
        }

        // Get the next token
        token = strtok(nullptr, delim);
        i++;
    }

    // Check if the number of values in the array is correct
    if (i != num_values) {
        throw range_error(
                "Error: the number of values in the input string is invalid, please use this example format as a reference 1:?:0.1");
    }

    return theta;
}

Dimension Configurations::CheckDimensionValue(const string &aDimension) {

    if (aDimension != "2D" and aDimension != "2d" and aDimension != "3D" and aDimension != "3d" and
        aDimension != "st" and aDimension != "ST") {
        throw range_error("Invalid value for Dimension. Please use 2D, 3D or ST.");
    }
    if (aDimension == "2D" or aDimension == "2d") {
        return Dimension2D;
    } else if (aDimension == "3D" or aDimension == "3d") {
        return Dimension3D;
    }
    return DimensionST;
}

int Configurations::CheckUnknownObservationsValue(const string &aValue) {
    int value = CheckNumericalValue(aValue);
    if (value >= GetProblemSize()) {
        throw range_error("Invalid value for ZmissNumber. Please make sure it's smaller than Problem size");
    }
    return value;
}

void Configurations::ParseDistanceMetric(const std::string &aDistanceMetric) {
    if (aDistanceMetric == "eg" || aDistanceMetric == "EG" || aDistanceMetric == "euclidean") {
        SetDistanceMetric(EUCLIDEAN_DISTANCE);
    } else if (aDistanceMetric == "gcd" || aDistanceMetric == "GCD" || aDistanceMetric == "great_circle") {
        SetDistanceMetric(GREAT_CIRCLE_DISTANCE);
    } else {
        throw range_error("Invalid value. Please use eg or gcd values only.");
    }
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
    Verbose temp = this->GetVerbosity();
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
    mDictionary["Tolerance"] = pow(10, -1 * aTolerance);
}
