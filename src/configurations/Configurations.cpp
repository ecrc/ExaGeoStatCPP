
// Copyright (c) 2017-2023 King Abdullah University of Science and Technology,
// All rights reserved.
// ExaGeoStat is a software package, provided by King Abdullah University of Science and Technology (KAUST).

/**
 * @file Configurations.cpp
 * @brief This file defines the Configurations class which stores the configuration parameters for ExaGeoStat.
 * @version 1.0.0
 * @author Sameh Abdulah
 * @author Mahmoud ElKarargy
 * @date 2023-01-31
**/

#include <iostream>
#include <algorithm>
#include <cstring>

#include <configurations/Configurations.hpp>

using namespace std;

using namespace exageostat::configurations;
using namespace exageostat::common;

RunMode Configurations::mRunMode = RunMode::STANDARD_MODE;


Configurations::Configurations() {

    // Set default values for arguments!
    SetComputation(common::EXACT_DENSE);
    SetCoresNumber(1);
    SetGPUsNumbers(0);
    SetPGrid(1);
    SetQGrid(1);
    SetMaxRank(1);
    SetIsOOC(false);
    SetKernelName("");
    SetDimension(common::Dimension2D);
    SetTimeSlot(1);
    SetProblemSize(0);
#ifdef EXAGEOSTAT_USE_CHAMELEON
    SetDenseTileSize(0);
#endif
#ifdef EXAGEOSTAT_USE_HiCMA
    SetLowTileSize(0);
#endif
    SetLoggerPath("");
    SetIsSynthetic(true);
    std::vector<double> theta;
    SetInitialTheta(theta);
    SetLowerBounds(theta);
    SetUpperBounds(theta);
    SetTargetTheta(theta);
    SetP(1);
    SetSeed(0);
    SetLogger(false);
    SetUnknownObservationsNb(0);
    SetMeanSquareError(0.0);
    SetApproximationMode(1);
    SetActualObservationsFilePath("");
    SetRecoveryFile("");
    SetPrecision(common::DOUBLE);
}


void Configurations::InitializeArguments(int aArgC, char **apArgV) {

    this->mArgC = aArgC;
    this->mpArgV = apArgV;

    // Get the example name
    string example_name = apArgV[0];
    // Remove the './'
    example_name.erase(0, 2);
    cout << "Running " << example_name << endl;

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
            } else if (argument_name == "--PGrid" || argument_name == "--pGrid" || argument_name == "--pgrid" ||
                       argument_name == "--p_grid") {
                SetPGrid(CheckNumericalValue(argument_value));
            } else if (argument_name == "--QGrid" || argument_name == "--qGrid" || argument_name == "--qgrid" ||
                       argument_name == "--q_grid") {
                SetQGrid(CheckNumericalValue(argument_value));
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
            } else if (argument_name == "--Gpus" || argument_name == "--GPUsNumbers" ||
                       argument_name == "--gpu_number" || argument_name == "--gpus") {
                SetGPUsNumbers(CheckNumericalValue(argument_value));
            } else if (argument_name == "--DTS" || argument_name == "--dts" || argument_name == "--Dts") {
                SetDenseTileSize(CheckNumericalValue(argument_value));
            } else if (argument_name == "--LTS" || argument_name == "--lts" || argument_name == "--Lts") {
                SetLowTileSize(CheckNumericalValue(argument_value));
            } else if (argument_name == "--maxRank" || argument_name == "--maxrank" ||
                       argument_name == "--max_rank") {
                SetMaxRank(CheckNumericalValue(argument_value));
            } else if (argument_name == "--ObservationsFile" || argument_name == "--observationsfile" ||
                       argument_name == "--observations_file") {
                SetActualObservationsFilePath(argument_value);
            } else if (argument_name == "--Seed" || argument_name == "--seed") {
                SetSeed(CheckNumericalValue(argument_value));
            } else if (argument_name == "--runmode" || argument_name == "--runMode" ||
                       argument_name == "--run_mode") {
                ParseRunMode(argument_value);
            } else if (argument_name == "--logpath" || argument_name == "--log_path" ||
                       argument_name == "--logPath") {
                SetLoggerPath(argument_value);
            } else if (argument_name == "--initial_theta" || argument_name == "--itheta" ||
                       argument_name == "--iTheta") {
                std::vector<double> theta = ParseTheta(argument_value);
                SetInitialTheta(theta);
            } else if (argument_name == "--lb" || argument_name == "--olb" || argument_name == "--lower_bounds") {
                std::vector<double> theta = ParseTheta(argument_value);
                SetLowerBounds(theta);
                SetStartingTheta(theta);
            } else if (argument_name == "--ub" || argument_name == "--oub" || argument_name == "--upper_bounds") {
                std::vector<double> theta = ParseTheta(argument_value);
                SetUpperBounds(theta);
            } else {
                if (!(argument_name == "--Dimension" || argument_name == "--dimension" ||
                      argument_name == "--dim" || argument_name == "--Dim" || argument_name == "--ZmissNumber" ||
                      argument_name == "--Zmiss" || argument_name == "--initial_theta" ||
                      argument_name == "--itheta" || argument_name == "--iTheta" ||
                      argument_name == "--target_theta" || argument_name == "--ttheta" ||
                      argument_name == "--tTheta" || argument_name == "--iterations" ||
                      argument_name == "--Iterations" ||
                      argument_name == "--max_mle_iterations" || argument_name == "--maxMleIterations" ||
                      argument_name == "--tolerance" ||
                      argument_name == "--distanceMetric" || argument_name == "--distance_metric" ||
                      argument_name == "--log_file_name" || argument_name == "--logFileName")) {
                    cout << "!! " << argument_name << " !!" << endl;
                    throw invalid_argument(
                            "This argument is undefined, Please use --help to print all available arguments");
                }
            }
        } else {
            if (argument_name == "--help") {
                PrintUsage();
            }
            if (argument_name == "--OOC") {
                SetIsOOC(true);
            } else if (argument_name == "--ApproximationMode" || argument_name == "--approximationmode" ||
                       argument_name == "--approximation_mode") {
                SetApproximationMode(true);
            } else if (argument_name == "--log" || argument_name == "--Log") {
                SetLogger(true);
            } else {
                if (!(argument_name == "--syntheticData" || argument_name == "--SyntheticData" ||
                      argument_name == "--synthetic_data" || argument_name == "--synthetic"))
                    throw invalid_argument(
                            "This argument is undefined, Please use --help to print all available arguments");
            }
        }
    }

    // Throw Errors if any of these arguments aren't given by the user.
    if (GetProblemSize() == 0) {
        throw domain_error("You need to set the problem size, before starting");
    }
#ifdef EXAGEOSTAT_USE_CHAMELEON
    if (GetDenseTileSize() == 0) {
        throw domain_error("You need to set the Dense tile size, before starting");
    }
#endif
#ifdef EXAGEOSTAT_USE_HiCMA
    if (GetLowTileSize() == 0) {
        throw domain_error("You need to set the Low tile size, before starting");
    }
#endif
    // Throw Errors if any of these arguments aren't given by the user.
    if (GetKernelName().empty()) {
        throw domain_error("You need to set the Kernel, before starting");
    }

    if (GetLogger()) {
        //initlog
        InitLog();
    }
}

void Configurations::InitializeDataGenerationArguments() {

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
            if (argument_name == "--Dimension" || argument_name == "--dimension"
                || argument_name == "--dim" || argument_name == "--Dim") {
                SetDimension(CheckDimensionValue(argument_value));
            } else if (argument_name == "--ZmissNumber" || argument_name == "--Zmiss") {
                SetUnknownObservationsNb(CheckUnknownObservationsValue(argument_value));
            } else if (argument_name == "--target_theta" || argument_name == "--ttheta" ||
                       argument_name == "--tTheta") {
                std::vector<double> theta = ParseTheta(argument_value);
                SetTargetTheta(theta);
            }
        } else {
            if (argument_name == "--syntheticData" || argument_name == "--SyntheticData" ||
                argument_name == "--synthetic_data" || argument_name == "--synthetic") {
                SetIsSynthetic(true);
            }
        }
    }
}

void Configurations::InitializeDataModelingArguments() {

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
            if (argument_name == "--distance_metric" || argument_name == "--distanceMetric") {
                ParseDistanceMetric(argument_value);
            } else if (argument_name == "--max_mle_iterations" || argument_name == "--maxMleIterations") {
                SetMaxMleIterations(CheckNumericalValue(argument_value));
            } else if (argument_name == "--tolerance") {
                SetTolerance(CheckNumericalValue(argument_value));
            } else if (argument_name == "--log_file_name" || argument_name == "--logFileName") {
                if (!GetLogger()) {
                    throw domain_error(
                            "To enable logging, please utilize the '--log' option in order to specify a log file.");
                }
                SetFileLogName(argument_value);
            }
        }
    }
}

void Configurations::InitializeDataPredictionArguments() {
}

void Configurations::PrintUsage() {
    cout << "\n\t*** Available Arguments For Synthetic Data Configurations***" << endl;
    cout << "\t\t --N=value : Problem size." << endl;
    cout << "\t\t --kernel=value : Used Kernel." << endl;
    cout << "\t\t --dimension=value : Used Dimension." << endl;
    cout << "\t\t --p_grid=value : Used P-Grid." << endl;
    cout << "\t\t --q_grid=value : Used P-Grid." << endl;
    cout << "\t\t --time_slot=value : Time slot value for ST." << endl;
    cout << "\t\t --computation=value : Used computation." << endl;
    cout << "\t\t --precision=value : Used precision." << endl;
    cout << "\t\t --cores=value : Used to set the number of cores." << endl;
    cout << "\t\t --gpus=value : Used to set the number of GPUs." << endl;
    cout << "\t\t --dts=value : Used to set the Dense Tile size." << endl;
    cout << "\t\t --lts=value : Used to set the Low Tile size." << endl;
    cout << "\t\t --Zmiss=value : Used to set number of unknown observation to be predicted." << endl;
    cout << "\t\t --observations_file=PATH/TO/File : Used to path the observations file path." << endl;
    cout << "\t\t --max_rank=value : Used to the max rank value." << endl;
    cout << "\t\t --olb=value : Lower bounds for optimization." << endl;
    cout << "\t\t --oub=value : Upper bounds for optimization." << endl;
    cout << "\t\t --itheta=value : Initial theta parameters for optimization." << endl;
    cout << "\t\t --ttheta=value : Target kernel parameters for optimization." << endl;
    cout << "\t\t --seed=value : Seed value for random number generation." << endl;
    cout << "\t\t --run_mode=value : Run mode whether verbose/not." << endl;
    cout << "\t\t --log_path=value : Path to log file." << endl;
    cout << "\t\t --distance_metric=value : Used distance metric either eg or gcd.";
    cout << "\t\t --max_mle_iterations=value : Maximum number of MLE iterations.";
    cout << "\t\t --tolerance : MLE tolerance between two iterations.";

    cout << "\t\t --synthetic_data : Used to enable generating synthetic data." << endl;
    cout << "\t\t --OOC : Used to enable Out of core technology." << endl;
    cout << "\t\t --approximation_mode : Used to enable Approximation mode." << endl;
    cout << "\t\t --log : Enable logging." << endl;
    cout << "\n\n";

    exit(0);
}

RunMode Configurations::GetRunMode() {
    return Configurations::mRunMode;
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

    if (aValue != "exact" and aValue != "Exact" and aValue != "Dense" and aValue != "dense"
        and aValue != "diag_approx" and aValue != "diagonal_approx"
        and aValue != "lr_approx" and aValue != "tlr" and aValue != "TLR") {
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

    if (aValue != "single" and aValue != "Single"
        and aValue != "double" and aValue != "Double"
        and aValue != "mix" and aValue != "Mix" and aValue != "Mixed" and aValue != "mixed") {
        throw range_error("Invalid value for Computation. Please use Single, Double or Mixed.");
    }
    if (aValue == "single" or aValue == "Single") {
        return SINGLE;
    } else if (aValue == "double" or aValue == "Double") {
        return DOUBLE;
    }
    return MIXED;
}

void Configurations::ParseRunMode(const std::string &aRunMode) {
    if (aRunMode == "verbose" || aRunMode == "Verbose") {
        mRunMode = RunMode::VERBOSE_MODE;
    } else if (aRunMode == "standard" || aRunMode == "Standard") {
        mRunMode = RunMode::STANDARD_MODE;
    } else {
        cout << "Error: " << aRunMode << " is not valid " << endl;
        throw range_error("Invalid value. Please use verbose or standard values only.");
    }
}


void Configurations::CheckKernelValue(const string &aKernel) {

    // Check if the kernel name exists in the availableKernels set.
    if (availableKernels.count(aKernel) <= 0) {
        throw range_error("Invalid value for Kernel. Please check manual.");
    } else {
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

std::vector<double> Configurations::ParseTheta(const std::string &aInputValues) {
    // Count the number of values in the string
    int num_values = 1;
    for (char aInputValue: aInputValues) {
        if (aInputValue == ':') {
            num_values++;
        }
    }

    // Allocate memory for the array of doubles
    std::vector<double> theta;

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
                cout << "Error: " << token << " is not a valid double or '?' " << endl;
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

    if (aDimension != "2D" and aDimension != "2d"
        and aDimension != "3D" and aDimension != "3d"
        and aDimension != "st" and aDimension != "ST") {
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
    if (aDistanceMetric == "eg" || aDistanceMetric == "EG") {
        SetDistanceMetric(EUCLIDIAN_DISTANCE);
    } else if (aDistanceMetric == "gcd" || aDistanceMetric == "GCD") {
        SetDistanceMetric(GREAT_CIRCLE_DISTANCE);
    } else {
        throw range_error("Invalid value. Please use eg or gcd values only.");
    }
}

void Configurations::InitLog() {
    try {
        SetFileLogPath(fopen(GetFileLogName().c_str(), "w+"));
    }
    catch (std::exception &e) {
        SetFileLogPath(fopen("log_file", "w+"));
    }
    fprintf(GetFileLogPath(), "\t\tlog file is generated by ExaGeoStat application\n");
    fprintf(GetFileLogPath(), "\t\t============================================\n");
}