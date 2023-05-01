
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


SyntheticDataConfigurations::SyntheticDataConfigurations(int argc, char **argv) {
    this->InitializeArguments(argc, argv);
}

SyntheticDataConfigurations::SyntheticDataConfigurations(const string& JSON_path) {
    // Not implemented yet
}

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
            } else if (argument_name == "--PGrid" || argument_name == "--pGrid" || argument_name == "--pgrid") {
                SetPGrid(CheckNumericalValue(argument_value));
            } else if (argument_name == "--QGrid" || argument_name == "--qGrid" || argument_name == "--qgrid") {
                SetQGrid(CheckNumericalValue(argument_value));
            } else if (argument_name == "--TimeSlot" || argument_name == "--timeslot") {
                SetTimeSlot( CheckNumericalValue(argument_value));
            } else if (argument_name == "--Computation" || argument_name == "--computation") {
                SetComputation(CheckComputationValue(argument_value));
            } else if (argument_name == "--precision" || argument_name == "--Precision") {
                SetPrecision(CheckPrecisionValue(argument_value));
            } else if (argument_name == "--cores" || argument_name == "--coresNumber") {
                SetCoresNumber(CheckNumericalValue(argument_value));
            } else if (argument_name == "--Gpus" || argument_name == "--GPUsNumbers") {
                SetGPUsNumber(CheckNumericalValue(argument_value));
            } else if (argument_name == "--DTS" || argument_name == "--dts" || argument_name == "--Dts") {
                SetDenseTileSize(CheckNumericalValue(argument_value));
            } else if (argument_name == "--LTS" || argument_name == "--lts" || argument_name == "--Lts") {
                SetLowTileSize(CheckNumericalValue(argument_value));
            } else if (argument_name == "--maxRank" || argument_name == "--maxrank") {
                SetMaxRank(CheckNumericalValue(argument_value));
            } else if (argument_name == "--ZmissNumber" || argument_name == "--Zmiss") {
                SetUnknownObservationsNb(CheckUnknownObservationsValue(argument_value));
            } else if (argument_name == "--ObservationsFile" || argument_name == "--observationsfile") {
                SetActualObservationsFilePath(argument_value);
            }
            else {
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
            if (argument_name == "--syntheticData" || argument_name == "--SyntheticData" || argument_name == "--synthetic") {
                SetIsSynthetic(true);
            } else if (argument_name == "--OOC") {
                SetIsOOC(true);
            } else if (argument_name == "--ApproximationMode" || argument_name == "--approximationmode") {
                SetApproximationMode(true);
            } else {
                throw invalid_argument(
                        "This argument is undefined, Please use --help to print all available arguments");
            }
        }
    }
    // Throw Errors if any of these arguments aren't given by the user.
    if(GetProblemSize() == 0){
        throw domain_error("You need to set the problem size, before starting");
    }
#ifdef EXAGEOSTAT_USE_CHAMELEON
    if(GetDenseTileSize() == 0){
        throw domain_error("You need to set the Dense tile size, before starting");
    }
#endif
#ifdef EXAGEOSTAT_USE_HiCMA
    if(GetLowTileSize() == 0){
        throw domain_error("You need to set the Low tile size, before starting");
    }
#endif
    if(GetKernel().empty()){
        throw domain_error("You need to set the Kernel, before starting");
    }
}

void SyntheticDataConfigurations::PrintUsage() {
    cout << "\n\t*** Available Arguments For Synthetic Data Configurations***" << endl;
    cout << "\t\t --N=value : Problem size." << endl;
    cout << "\t\t --Kernel=value : Used Kernel." << endl;
    cout << "\t\t --Dimension=value : Used Dimension." << endl;
    cout << "\t\t --PGrid=value : Used P-Grid." << endl;
    cout << "\t\t --QGrid=value : Used P-Grid." << endl;
    cout << "\t\t --TimeSlot=value : Time slot value for ST." << endl;
    cout << "\t\t --Computation=value : Used computation" << endl;
    cout << "\t\t --Precision=value : Used precision" << endl;
    cout << "\t\t --CoresNumber=value : Used to set the number of cores." << endl;
    cout << "\t\t --GPUsNumber=value : Used to set the number of GPUs." << endl;
    cout << "\t\t --Dts=value : Used to set the Dense Tile size." << endl;
    cout << "\t\t --Lts=value : Used to set the Low Tile size." << endl;
    cout << "\t\t --MaxRank=value : Used to the max rank value." << endl;
    cout << "\t\t --ZmissNumber=value : Used to set number of unknown observation to be predicted." << endl;
    cout << "\t\t --ObservationsFile=PATH/TO/File : Used to path the observations file path." << endl;
    cout << "\t\t --SyntheticData : Used to enable generating synthetic data." << endl;
    cout << "\t\t --OOC : Used to enable Out of core technology." << endl;
    cout << "\t\t --ApproximationMode : Used to enable Approximation mode." << endl;
    cout << "\n\n";

    exit(0);

}

Dimension SyntheticDataConfigurations::GetDimension() {
    return this->mDimension;
}

void SyntheticDataConfigurations::SetDimension(Dimension aDimension) {
    this->mDimension = aDimension;
}

Dimension SyntheticDataConfigurations::CheckDimensionValue(const string& aDimension) {

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

int SyntheticDataConfigurations::CheckUnknownObservationsValue(const string& aValue) {
    int value = CheckNumericalValue(aValue);
    if (value >= GetProblemSize()){
        throw range_error("Invalid value for ZmissNumber. Please make sure it's smaller than Problem size");
    }
    return value;
}