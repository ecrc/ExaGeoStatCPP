
/*
 * Copyright (c) 2017-2023 King Abdullah University of Science and Technology,
 * Copyright (C) 2023 by Brightskies inc,
 * All rights reserved.
 * ExaGeoStat is a software package, provided by King Abdullah University of Science and Technology (KAUST).
 */

/**
 * @file SyntheticDataConfigurations.cpp
 * @brief Implementation for Synthetic data Configurations.
 * @version 1.0.0
 * @author Sameh Abdulah
 * @date 2023-02-01
**/

#include <configurations/data-generation/concrete/SyntheticDataConfigurations.hpp>
#include <iostream>

using namespace exageostat::configurations::data_configurations;
using namespace exageostat::common;
using namespace std;

//// TODO: Add all arguments info in README.
//// TODO: Add all supported kernels in README.
void SyntheticDataConfigurations::InitializeArguments(int argc, char **argv) {
    // Get the example name
    string example_name = argv[0];
    // Remove the './'
    example_name.erase(0, 2);
    cout << "Running " << example_name << endl;

    string argument;
    string argument_name;
    string argument_value;
    int equal_sign_Idx;

    //Set verbosity level with default value = standard
    SetRunMode(exageostat::common::RunMode::STANDARD_MODE);

    // Loop through the arguments
    for (int i = 1; i < argc; ++i) {
        argument = argv[i];
        equal_sign_Idx = argument.find('=');
        argument_name = argument.substr(0, equal_sign_Idx);

        // Check if argument has an equal sign.
        if (equal_sign_Idx != string::npos) {
            argument_value = argument.substr(equal_sign_Idx + 1);

            // Check the argument name and set the corresponding value
            if (argument_name == "--N" || argument_name == "--n") {
                SetProblemSize(CheckNumericalValue(argument_value));
            } else if (argument_name == "--Kernel" || argument_name == "--kernel") {
                CheckKernelValue(argument_value);
            } else if (argument_name == "--Dimension" || argument_name == "--dimension"
                       || argument_name == "--dim" || argument_name == "--Dim") {
                SetDimension(CheckDimensionValue(argument_value));
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
                       argument_name == "--cores_number") {
                SetCoresNumber(CheckNumericalValue(argument_value));
            } else if (argument_name == "--Gpus" || argument_name == "--GPUsNumbers" ||
                       argument_name == "--gpu_number" || argument_name == "--gpus") {
                SetGPUsNumber(CheckNumericalValue(argument_value));
            } else if (argument_name == "--DTS" || argument_name == "--dts" || argument_name == "--Dts") {
                SetDenseTileSize(CheckNumericalValue(argument_value));
            } else if (argument_name == "--LTS" || argument_name == "--lts" || argument_name == "--Lts") {
                SetLowTileSize(CheckNumericalValue(argument_value));
            } else if (argument_name == "--maxRank" || argument_name == "--maxrank" || argument_name == "--max_rank") {
                SetMaxRank(CheckNumericalValue(argument_value));
            } else if (argument_name == "--ZmissNumber" || argument_name == "--Zmiss") {
                SetUnknownObservationsNb(CheckUnknownObservationsValue(argument_value));
            } else if (argument_name == "--ObservationsFile" || argument_name == "--observationsfile" ||
                       argument_name == "--observations_file") {
                SetActualObservationsFilePath(argument_value);
            } else if (argument_name == "--lb" || argument_name == "--olb" || argument_name == "--lowerBounds") {
                SetLowerBounds(ParseTheta(argument_value));
            } else if (argument_name == "--ub" || argument_name == "--oub" || argument_name == "--upperBounds") {
                SetUpperBounds(ParseTheta(argument_value));
            } else if (argument_name == "--initialTheta" || argument_name == "--itheta" ||
                       argument_name == "--iTheta") {
                SetInitialTheta(ParseTheta(argument_value));
            } else if (argument_name == "--targetTheta" || argument_name == "--ttheta" || argument_name == "--tTheta") {
                SetTargetTheta(ParseTheta(argument_value));
            } else if (argument_name == "--Seed" || argument_name == "--seed") {
                SetSeed(CheckNumericalValue(argument_value));
            } else if (argument_name == "--runmode" || argument_name == "--runMode" || argument_name == "--run_mode") {
                ParseRunMode(argument_value);
            } else if (argument_name == "--logpath" || argument_name == "--log_path" || argument_name == "--logPath") {
                SetLoggerPath(argument_value);
            } else {
                cout << "!! " << argument_name << " !!" << endl;
                throw invalid_argument(
                        "This argument is undefined, Please use --help to print all available arguments");
            }
        }
            // If none then, just set the argument to True.
        else {
            if (argument_name == "--help") {
                PrintUsage();
            }
            if (argument_name == "--syntheticData" || argument_name == "--SyntheticData" ||
                argument_name == "--synthetic_data" || argument_name == "--synthetic") {
                SetIsSynthetic(true);
            } else if (argument_name == "--OOC") {
                SetIsOOC(true);
            } else if (argument_name == "--ApproximationMode" || argument_name == "--approximationmode" ||
                       argument_name == "--approximation_mode") {
                SetApproximationMode(true);
            } else if (argument_name == "--log" || argument_name == "--Log") {
                SetLogger(true);
            } else {
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
    if (GetKernel().empty()) {
        throw domain_error("You need to set the Kernel, before starting");
    }
}

void SyntheticDataConfigurations::PrintUsage() {
    cout << "\n\t*** Available Arguments For Synthetic Data Configurations***" << endl;
    cout << "\t\t --N=value : Problem size." << endl;
    cout << "\t\t --kernel=value : Used Kernel." << endl;
    cout << "\t\t --dimension=value : Used Dimension." << endl;
    cout << "\t\t --p_grid=value : Used P-Grid." << endl;
    cout << "\t\t --q_grid=value : Used P-Grid." << endl;
    cout << "\t\t --time_slot=value : Time slot value for ST." << endl;
    cout << "\t\t --computation=value : Used computation" << endl;
    cout << "\t\t --precision=value : Used precision" << endl;
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

    cout << "\t\t --synthetic_data : Used to enable generating synthetic data." << endl;
    cout << "\t\t --OOC : Used to enable Out of core technology." << endl;
    cout << "\t\t --approximation_mode : Used to enable Approximation mode." << endl;
    cout << "\t\t --log : Enable logging." << endl;

    cout << "\n\n";

    exit(0);

}

Dimension SyntheticDataConfigurations::GetDimension() {
    return this->mDimension;
}

void SyntheticDataConfigurations::SetDimension(Dimension aDimension) {
    this->mDimension = aDimension;
}

Dimension SyntheticDataConfigurations::CheckDimensionValue(const string &aDimension) {

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

int SyntheticDataConfigurations::CheckUnknownObservationsValue(const string &aValue) {
    int value = CheckNumericalValue(aValue);
    if (value >= GetProblemSize()) {
        throw range_error("Invalid value for ZmissNumber. Please make sure it's smaller than Problem size");
    }
    return value;
}
