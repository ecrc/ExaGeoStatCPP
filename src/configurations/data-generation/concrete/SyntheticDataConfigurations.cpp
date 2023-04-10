
/*
 * Copyright (c) 2017-2023 King Abdullah University of Science and Technology,
 * Copyright (C) 2023 by Brightskies inc,
 * All rights reserved.
 * ExaGeoStat is a software package, provided by King Abdullah University of Science and Technology (KAUST).
 */

/**
 * @file SyntheticDataConfigurations.cpp
 * Implementation for Synthetic data Configurations.
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

SyntheticDataConfigurations::SyntheticDataConfigurations(string JSON_path) {

}

void SyntheticDataConfigurations::InitializeArguments(int argc, char **argv) {

    string example_name = argv[0];
    // Removes the './'
    example_name.erase(0, 2);
    cout << "Running " << example_name << endl;

    string argument;
    string argumentName;
    string argumentValue;
    int equalSignIdx;

    // Skipping first argument as it's the example name.
    for (int i = 1; i < argc; ++i) {

        argument = argv[i];
        equalSignIdx = argument.find('=');
        argumentName = argument.substr(0, equalSignIdx);

        // Check if argument has an equal sign.
        if (equalSignIdx != string::npos) {
            argumentValue = argument.substr(equalSignIdx + 1);

            if (argumentName == "--N" || argumentName == "--n") {
                SetProblemSize(CheckNumericalValue(argumentValue));
            } else if (argumentName == "--Kernel" || argumentName == "--kernel") {
                CheckKernelValue(argumentValue);
                SetKernel(argumentValue);
            } else if (argumentName == "--Dimension" || argumentName == "--dimension"
                       || argumentName == "--dim" || argumentName == "--Dim") {
                SetDimension(CheckDimensionValue(argumentValue));
            } else if (argumentName == "--PGrid" || argumentName == "--pGrid" || argumentName == "--pgrid") {
                SetPGrid(CheckNumericalValue(argumentValue));
            }else if (argumentName == "--QGrid" || argumentName == "--qGrid" || argumentName == "--qgrid") {
                SetQGrid(CheckNumericalValue(argumentValue));
            } else if (argumentName == "--TimeSlot" || argumentName == "--timeslot") {
                SetTimeSlot( CheckNumericalValue(argumentValue));
            } else if (argumentName == "--Computation" || argumentName == "--computation") {
                SetComputation(CheckComputationValue(argumentValue));
            } else if (argumentName == "--precision" || argumentName == "--Precision") {
                SetPrecision(CheckPrecisionValue(argumentValue));
            } else if (argumentName == "--cores" || argumentName == "--coresNumber") {
                SetCoresNumber(CheckNumericalValue(argumentValue));
            } else if (argumentName == "--Gpus" || argumentName == "--GPUsNumbers") {
                SetGPUsNumber(CheckNumericalValue(argumentValue));
            } else if (argumentName == "--DTS" || argumentName == "--dts" || argumentName == "--Dts") {
                SetDenseTileSize(CheckNumericalValue(argumentValue));
            } else if (argumentName == "--LTS" || argumentName == "--lts" || argumentName == "--Lts") {
                SetLowTileSize(CheckNumericalValue(argumentValue));
            } else if (argumentName == "--maxRank" || argumentName == "--maxrank") {
                SetMaxRank(CheckNumericalValue(argumentValue));
            } else if (argumentName == "--ZmissNumber" || argumentName == "--Zmiss") {
                SetUnknownObservationsNb(CheckUnknownObservationsValue(argumentValue));
            } else if (argumentName == "--ObservationsFile" || argumentName == "--observationsfile") {
                SetActualObservationsFilePath(argumentValue);
            }
            else {
                cout << "!! " << argumentName << " !!" << endl;
                throw invalid_argument(
                        "This argument is undefined, Please use --help to print all available arguments");
            }
        }
        // If none then, just set the argument to True.
        else {
            if (argumentName == "--help") {
                PrintUsage();
            }
            if (argumentName == "--syntheticData") {
                SetIsSynthetic(true);
            } else if (argumentName == "--OOC") {
                SetIsOOC(true);
            } else if (argumentName == "--ApproximationMode" || argumentName == "--approximationmode") {
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

Dimension SyntheticDataConfigurations::CheckDimensionValue(string aDimension) {

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

int SyntheticDataConfigurations::CheckUnknownObservationsValue(string aValue) {
    int value = CheckNumericalValue(aValue);
    if (value >= GetProblemSize()){
        throw range_error("Invalid value for ZmissNumber. Please make sure it's smaller than Problem size");
    }
    return value;
}