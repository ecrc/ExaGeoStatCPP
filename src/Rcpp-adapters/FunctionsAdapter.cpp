
// Copyright (c) 2017-2023 King Abdullah University of Science and Technology,
// All rights reserved.
// ExaGeoStat is a software package, provided by King Abdullah University of Science and Technology (KAUST).

/**
 * @file FunctionsAdapter.cpp
 * @brief Header file for function adapters in the ExaGeoStat software.
 * @version 1.1.0
 * @author Mahmoud ElKarargy
 * @date 2024-01-29
**/

#include <cstring>

#include <Rcpp-adapters/FunctionsAdapter.hpp>
#include <data-units/ExaGeoStatData.hpp>
#include <api/ExaGeoStat.hpp>

using namespace std;
using namespace Rcpp;

using namespace exageostat::dataunits;
using namespace exageostat::api;

namespace exageostat::adapters {

    Configurations *
    R_InitializeArguments(const int &aProblemSize, const string &aKernelName, const vector<int> &aTileSize,
                          const vector<int> &aP_QGrid, const int &aTimeSlot, const string &aComputation,
                          const string &aPrecision, const vector<int> &aCoresGPUsNumber, const int &aBand,
                          const int &aMaxRank, const vector<double> &aInitialTheta,
                          const vector <vector<double>> &aLowerUpperBounds, const vector<double> &aEstimatedTheta,
                          const string &aVerbose, const string &aDimension, const int &aMaxMleIterations,
                          const double &aTolerance, const vector<int> &aPrediction, const std::vector<std::string> &aPath, const bool &aSaveData) {

        vector<string> argStrings;

        // Program name or placeholder
        argStrings.emplace_back("--ConfigurationsModuleR");
        // Convert each argument to string and add to argv
        argStrings.push_back("--N=" + to_string(aProblemSize));
        argStrings.push_back("--kernel=" + aKernelName);
        argStrings.push_back("--dts=" + to_string(aTileSize[0]));
        argStrings.push_back("--lts=" + to_string(aTileSize[1]));
        argStrings.push_back("--p=" + to_string(aP_QGrid[0]));
        argStrings.push_back("--q=" + to_string(aP_QGrid[1]));
        argStrings.push_back("--time_slot=" + to_string(aTimeSlot));
        argStrings.push_back("--computation=" + aComputation);
        argStrings.push_back("--precision=" + aPrecision);
        argStrings.push_back("--cores=" + to_string(aCoresGPUsNumber[0]));
        argStrings.push_back("--gpus=" + to_string(aCoresGPUsNumber[1]));
        argStrings.push_back("--band=" + to_string(aBand));
        argStrings.push_back("--max_rank=" + to_string(aMaxRank));

        string initial_theta = "--iTheta=";
        for(int i = 0; i < aInitialTheta.size(); i++){
            if (i != aInitialTheta.size() - 1){
                initial_theta += to_string(aInitialTheta[i]) + ":";
            } else{
                initial_theta += to_string(aInitialTheta[i]);
            }
        }
        argStrings.push_back(initial_theta);
        string lower_bound = "--lb=";
        string upper_bound = "--ub=";
        for(int i = 0; i < aLowerUpperBounds[0].size(); i++){
            if (i != aLowerUpperBounds[0].size() - 1) {
                lower_bound += to_string(aLowerUpperBounds[0][i]) + ":";
            }
            else{
                lower_bound += to_string(aLowerUpperBounds[0][i]);
            }
        }
        for(int i = 0; i < aLowerUpperBounds[1].size(); i++){
            if (i != aLowerUpperBounds[1].size() - 1) {
                upper_bound += to_string(aLowerUpperBounds[1][i]) + ":";
            }
            else{
                upper_bound += to_string(aLowerUpperBounds[1][i]);
            }
        }
        argStrings.push_back(lower_bound);
        argStrings.push_back(upper_bound);

        string estimated_theta = "--eTheta=";
        for(int i = 0; i < aInitialTheta.size(); i++){
            if (i != aInitialTheta.size() - 1) {
                if(i < aEstimatedTheta.size()){
                    estimated_theta += to_string(aEstimatedTheta[i]) + ":";
                }
                else{
                    estimated_theta += "-1:";
                }
            }
            else{
                if(i < aEstimatedTheta.size()) {
                    estimated_theta += to_string(aEstimatedTheta[i]);
                }
                else{
                    estimated_theta += "-1";
                }
            }
        }
        argStrings.push_back(estimated_theta);
        argStrings.push_back("--verbose=" + aVerbose);
        argStrings.push_back("--dimension=" + aDimension);
        argStrings.push_back("--max_mle_iterations=" + to_string(aMaxMleIterations));
        argStrings.push_back("--tolerance=" + to_string(aTolerance));
        if(!aPath[0].empty()){
            argStrings.push_back("--log_path=" + aPath[0]);
        }
        if(!aPath[1].empty()){
            argStrings.push_back("--data_path=" + aPath[1]);
        }
        if(!aPath[2].empty()){
            argStrings.push_back("--observations_file=" + aPath[2]);
        }
        if(!aPath[3].empty()){
            argStrings.push_back("--recovery_file=" + aPath[3]);
        }
        if(aSaveData){
            argStrings.emplace_back("--log");
        }
        // This means that ZMiss > 0, which means prediction is activated
        if (aPrediction[0] > 0) {
            argStrings.push_back("--ZMiss=" + to_string(aPrediction[0]));

            // if mspe is activated
            if (aPrediction[1]) {
                argStrings.emplace_back("--mspe");
            }

            // if idw is activated
            if (aPrediction[2]) {
                argStrings.emplace_back("--idw");
            }

            // if fisher is activated
            if (aPrediction[3]) {
                argStrings.emplace_back("--fisher");
            }

            // if mloe_mmom is activated
            if (aPrediction[4]) {
                argStrings.emplace_back("--mloe_mmom");
            }
        }
        // Now, allocate a char** array and copy the arguments into it.
        // Argv ownership is passed to Configurations, and it's responsibility for freeing it.
        char **argv = new char *[argStrings.size()];

        for (size_t i = 0; i < argStrings.size(); ++i) {
            argv[i] = new char[argStrings[i].size() + 1];  // +1 for the null terminator
            std::strcpy(argv[i], argStrings[i].c_str());
        }

        auto configurations = new Configurations();
        configurations->InitializeArguments(argStrings.size(), argv, true);
        return configurations;
    }

    ExaGeoStatData<double> *R_ExaGeoStatLoadData(ExaGeoStatHardware *apHardware, Configurations *apConfigurations,
                                                 ExaGeoStatData<double> *apData) {

        // Wrap the raw pointer in a unique_ptr
        std::unique_ptr<ExaGeoStatData<double>> apDataPtr(apData);
        exageostat::api::ExaGeoStat<double>::ExaGeoStatLoadData(*apHardware, *apConfigurations, apDataPtr);

        // Since ExaGeoStatLoadData takes a reference, apDataPtr still owns the object.
        // We can safely return the raw pointer, but we need to release it from apDataPtr to avoid deletion.
        return apDataPtr.release();
    }

    ExaGeoStatData<double> *R_ExaGeoStatModelData(ExaGeoStatHardware *apHardware, Configurations *apConfigurations,
                                                  ExaGeoStatData<double> *apData,
                                                  Nullable <NumericMatrix> aMeasurementsMatrix,
                                                  Nullable <NumericMatrix> aLocationsX,
                                                  Nullable <NumericMatrix> aLocationsY,
                                                  Nullable <NumericMatrix> aLocationsZ) {

        std::unique_ptr<ExaGeoStatData<double>> apDataPtr(apData);
        double* pMeasurementsMatrixPtr = nullptr;
        if (aMeasurementsMatrix == nullptr || aMeasurementsMatrix.isNull()){
            // This was the only way to pass C++ tests, if we checked for aMeasurementsMatrix != nullptr a seg fault occurs
        }
        else{

            NumericMatrix mat(aMeasurementsMatrix.get());
            pMeasurementsMatrixPtr = &mat[0]; // Get the raw pointer to the matrix data

            if (!aLocationsX.isNull()) {
                double *pLocationsXPtr;
                NumericMatrix mat_location(aLocationsX.get());
                pLocationsXPtr = &mat_location[0]; // Get the raw pointer to the matrix data
                apData->GetLocations()->SetLocationX(*pLocationsXPtr, mat_location.ncol());
            }

            if (aLocationsY.isNotNull()) {
                double *pLocationsYPtr;
                NumericMatrix mat_location(aLocationsY.get());
                pLocationsYPtr = &mat_location[0]; // Get the raw pointer to the matrix data
                apData->GetLocations()->SetLocationY(*pLocationsYPtr, mat_location.ncol());

            }

            if (aLocationsZ == nullptr || aLocationsZ.isNotNull()) {
                double *pLocationsZPtr;
                NumericMatrix mat_location(aLocationsZ.get());
                pLocationsZPtr = &mat_location[0]; // Get the raw pointer to the matrix data
                apData->GetLocations()->SetLocationZ(*pLocationsZPtr, mat_location.ncol());
            }
        }

        exageostat::api::ExaGeoStat<double>::ExaGeoStatDataModeling(*apHardware, *apConfigurations, apDataPtr, pMeasurementsMatrixPtr);
        return apDataPtr.release();
    }

    ExaGeoStatData<double> *R_ExaGeoStatPredictData(ExaGeoStatHardware *apHardware, Configurations *apConfigurations,
                                                    ExaGeoStatData<double> *apData,
                                                    Nullable <NumericMatrix> aMeasurementsMatrix,
                                                    Nullable <NumericMatrix> aLocationsX,
                                                    Nullable <NumericMatrix> aLocationsY,
                                                    Nullable <NumericMatrix> aLocationsZ) {

        std::unique_ptr<ExaGeoStatData<double>> apDataPtr(apData);
        double* pMeasurementsMatrixPtr = nullptr;
        if (aMeasurementsMatrix == nullptr || aMeasurementsMatrix.isNull()){
            // This was the only way to pass C++ tests, if we checked for aMeasurementsMatrix != nullptr a seg fault occurs
        }
        else{

            NumericMatrix mat(aMeasurementsMatrix.get());
            pMeasurementsMatrixPtr = &mat[0]; // Get the raw pointer to the matrix data

            if (!aLocationsX.isNull()) {
                double *pLocationsXPtr;
                NumericMatrix mat_location(aLocationsX.get());
                pLocationsXPtr = &mat_location[0]; // Get the raw pointer to the matrix data
                apData->GetLocations()->SetLocationX(*pLocationsXPtr, mat_location.ncol());
            }

            if (aLocationsY.isNotNull()) {
                double *pLocationsYPtr;
                NumericMatrix mat_location(aLocationsY.get());
                pLocationsYPtr = &mat_location[0]; // Get the raw pointer to the matrix data
                apData->GetLocations()->SetLocationY(*pLocationsYPtr, mat_location.ncol());

            }

            if (aLocationsZ.isNotNull()) {
                double *pLocationsZPtr;
                NumericMatrix mat_location(aLocationsZ.get());
                pLocationsZPtr = &mat_location[0]; // Get the raw pointer to the matrix data
                apData->GetLocations()->SetLocationZ(*pLocationsZPtr, mat_location.ncol());
            }
        }

        exageostat::api::ExaGeoStat<double>::ExaGeoStatPrediction(*apHardware, *apConfigurations, apDataPtr, pMeasurementsMatrixPtr);
        return apDataPtr.release();
    }
}