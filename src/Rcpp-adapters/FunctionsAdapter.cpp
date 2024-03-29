
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

#include <Rcpp-adapters/FunctionsAdapter.hpp>
#include <api/ExaGeoStat.hpp>
#include <results/Results.hpp>

using namespace std;
using namespace Rcpp;

using namespace exageostat::api;
using namespace exageostat::common;
using namespace exageostat::results;
using namespace exageostat::configurations;

namespace exageostat::adapters {

    NumericVector R_GetLocationX(ExaGeoStatData<double> *apData) {

        // Obtain the pointer to the array of doubles
        double *data = apData->GetLocations()->GetLocationX();
        int length = apData->GetLocations()->GetSize();
        // Create an empty NumericVector of the appropriate length
        NumericVector vec(length);
        // Copy data from the double array to the NumericVector
        copy(data, data + length, vec.begin());
        return vec;
    }

    NumericVector R_GetLocationY(ExaGeoStatData<double> *apData) {

        // Obtain the pointer to the array of doubles
        double *data = apData->GetLocations()->GetLocationY();
        int length = apData->GetLocations()->GetSize();
        // Create an empty NumericVector of the appropriate length
        NumericVector vec(length);
        // Copy data from the double array to the NumericVector
        copy(data, data + length, vec.begin());
        return vec;
    }

    NumericVector R_GetLocationZ(ExaGeoStatData<double> *apData) {

        // Obtain the pointer to the array of doubles
        double *data = apData->GetLocations()->GetLocationZ();
        int length = apData->GetLocations()->GetSize();
        // Create an empty NumericVector of the appropriate length
        NumericVector vec(length);
        // Copy data from the double array to the NumericVector
        copy(data, data + length, vec.begin());
        return vec;
    }

    NumericVector R_GetDescZValues(ExaGeoStatData<double> *apData, const string &aType) {

        DescriptorType descriptorType;
        void *pDescriptor;

        if (aType == "chameleon" || aType == "Chameleon") {
            descriptorType = CHAMELEON_DESCRIPTOR;
            pDescriptor = apData->GetDescriptorData()->GetDescriptor(descriptorType, DESCRIPTOR_Z).chameleon_desc;
        } else if (aType == "hicma" || aType == "Hicma" || aType == "HICMA") {
#ifdef USE_HICMA
            descriptorType = HICMA_DESCRIPTOR;
            pDescriptor = apData->GetDescriptorData()->GetDescriptor(descriptorType, DESCRIPTOR_Z).hicma_desc;
#else
            throw runtime_error("Please enable HiCMA to use HiCMA descriptors.");
#endif
        } else {
            throw domain_error("Invalid type of descriptor, please use chameleon or hicma.");
        }

        // Obtain the pointer to the array of doubles
        double *data = apData->GetDescriptorData()->GetDescriptorMatrix(descriptorType, pDescriptor);
        int length = apData->GetLocations()->GetSize();
        // Create an empty NumericVector of the appropriate length
        NumericVector vec(length);
        // Copy data from the double array to the NumericVector
        copy(data, data + length, vec.begin());
        return vec;
    }

    ExaGeoStatData<double> *
    R_ExaGeoStatLoadData(const string &aKernelName, const vector<double> &aInitialTheta, const string &aDistanceMatrix,
                         const int &aProblemSize, const int &aSeed, const int &aDenseTileSize, const int &aLowTileSize,
                         const string &aDimension, const string &aLogPath, const string &aDataPath,
                         const string &aRecoveryFilePath, const string &aObservationsFilePath) {

        // manually set the configurations values with the user input from R.
        Configurations configurations;
        configurations.SetProblemSize(aProblemSize);
        configurations.CheckKernelValue(aKernelName);
        configurations.ParseDistanceMetric(aDistanceMatrix);
        configurations.SetSeed(aSeed);
        configurations.SetDenseTileSize(aDenseTileSize);
        configurations.SetLowTileSize(aLowTileSize);
        configurations.SetComputation(EXACT_DENSE);
        configurations.SetInitialTheta(aInitialTheta);
        configurations.SetDimension(Configurations::CheckDimensionValue(aDimension));

        if (!aLogPath.empty()) {
            configurations.SetLogger(true);
            configurations.SetLoggerPath(aLogPath);
        }
        if (!aDataPath.empty()) {
            configurations.SetDataPath(aDataPath);
            configurations.SetIsSynthetic(false);
        }
        configurations.SetRecoveryFile(aRecoveryFilePath);
        configurations.SetActualObservationsFilePath(aObservationsFilePath);

        unique_ptr<ExaGeoStatData<double>> data;
        ExaGeoStat<double>::ExaGeoStatLoadData(configurations, data);

        // We can safely return the raw pointer, but we need to release it from data to avoid deletion.
        return data.release();
    }

    vector<double>
    R_ExaGeoStatModelData(const string &aComputation, const string &aKernelName, const string &aDistanceMatrix,
                          const vector<double> &aLowerBound, const vector<double> &aUpperBound, const int &aTolerance,
                          const int &aMleIterations, const int &aDenseTileSize, const int &aLowTileSize,
                          const string &aDimension, const int &aBand, const int &aMaxRank,
                          SEXP apData, Nullable <NumericVector> aMeasurementsVector,
                          Nullable <NumericVector> aLocationsX, Nullable <NumericVector> aLocationsY,
                          Nullable <NumericVector> aLocationsZ) {

        Configurations configurations;
        bool is_initialized = ((SEXP) apData == R_NilValue);
        unique_ptr<ExaGeoStatData<double>> data;
        if (!is_initialized) {
            auto temp_data = (ExaGeoStatData<double> *) Rcpp::internal::as_module_object_internal(apData);
            // Set the unique_ptr to manage the ExaGeoStatData pointer
            data.reset(temp_data);
        } else {
            NumericVector z_values(aMeasurementsVector.get());

            // Create a new ExaGeoStatData object with configurations and manage it with unique_ptr
            data = make_unique<ExaGeoStatData<double>>(z_values.size(), configurations.GetDimension());
        }
        double *pMeasurementsVectorPtr = GetDataFromArguments(aMeasurementsVector, aLocationsX, aLocationsY,
                                                              aLocationsZ, data, configurations,
                                                              aKernelName, aDistanceMatrix, aDenseTileSize,
                                                              aLowTileSize, aDimension,
                                                              Configurations::CheckComputationValue(aComputation));

        configurations.SetLowerBounds(aLowerBound);
        configurations.SetUpperBounds(aUpperBound);
        configurations.SetMaxMleIterations(aMleIterations);
        configurations.SetTolerance(aTolerance);
        configurations.SetBand(aBand);
        configurations.SetMaxRank(aMaxRank);

        ExaGeoStat<double>::ExaGeoStatDataModeling(configurations, data, pMeasurementsVectorPtr);
        // Take back ownership, to avoid deleting apData when the unique_ptr goes out of scope.  
        apData = (SEXP) data.release();
        return configurations.GetStartingTheta();
    }

    vector<double> R_ExaGeoStatPredictData(const string &aKernelName, const string &aDistanceMatrix,
                                 const vector<double> &aEstimatedTheta,
                                 const int &aDenseTileSize, const int &aLowTileSize,
                                 const string &aDimension,
                                 vector<vector<double>> &aTrainData,
                                 vector<vector<double>> &aTestData) {

        Configurations configurations;
        configurations.SetComputation(EXACT_DENSE);

        auto data = make_unique<ExaGeoStatData<double>>(aTrainData[0].size(), configurations.GetDimension());
        auto *train_locations = new dataunits::Locations<double>(aTrainData[0].size(), configurations.GetDimension());
        auto *test_locations = new dataunits::Locations<double>(aTestData[0].size(), configurations.GetDimension());

        auto* z_values = new double[aTrainData.back().size()];
        for(int i = 0; i< aTrainData[0].size(); i++){
            train_locations->SetLocationX(*aTrainData[0].data(), aTrainData[0].size());
            train_locations->SetLocationY(*aTrainData[1].data(), aTrainData[1].size());
            if(aTrainData.size() > 3) {
                train_locations->SetLocationZ(*aTrainData[2].data(), aTrainData[2].size());
            }
        }
        memcpy(z_values, aTrainData.back().data(), aTrainData.back().size() * sizeof(double));
        for(int i = 0; i< aTestData[0].size(); i++){
            test_locations->SetLocationX(*aTestData[0].data(), aTestData[0].size());
            test_locations->SetLocationY(*aTestData[1].data(), aTestData[1].size());
        }

        // Set common configurations.
        configurations.CheckKernelValue(aKernelName);
        configurations.ParseDistanceMetric(aDistanceMatrix);
        configurations.SetDenseTileSize(aDenseTileSize);
        configurations.SetLowTileSize(aLowTileSize);
        configurations.SetDimension(Configurations::CheckDimensionValue(aDimension));
        configurations.SetProblemSize(data->GetLocations()->GetSize());
        configurations.SetEstimatedTheta(aEstimatedTheta);
        configurations.SetIsMSPE(TRUE);

        ExaGeoStat<double>::ExaGeoStatPrediction(configurations, data, z_values, train_locations, test_locations);
        // Take back ownership, to avoid deleting apData when the unique_ptr goes out of scope.
        auto temp_data = data.release();
        delete[] z_values;

        return Results::GetInstance()->GetPredictedMissedValues();
    }

    vector<double> R_ExaGeoStatMLOE_MMOM(const string &aKernelName, const string &aDistanceMatrix,
                                         const vector<double> &aEstimatedTheta, const vector<double> &aTrueTheta,
                                         const int &aDenseTileSize, const int &aLowTileSize, const string &aDimension,
                                         vector<vector<double>> &aTrainData,
                                         vector<vector<double>> &aTestData) {

        Configurations configurations;
        configurations.SetComputation(EXACT_DENSE);

        auto data = make_unique<ExaGeoStatData<double>>(aTrainData[0].size(), configurations.GetDimension());
        auto *train_locations = new dataunits::Locations<double>(aTrainData[0].size(), configurations.GetDimension());
        auto *test_locations = new dataunits::Locations<double>(aTestData[0].size(), configurations.GetDimension());

        auto* z_values = new double[aTrainData.back().size()];
        for(int i = 0; i< aTrainData[0].size(); i++){
            train_locations->SetLocationX(*aTrainData[0].data(), aTrainData[0].size());
            train_locations->SetLocationY(*aTrainData[1].data(), aTrainData[1].size());
            if(aTrainData.size() > 3) {
                train_locations->SetLocationZ(*aTrainData[2].data(), aTrainData[2].size());
            }
        }
        data->SetLocations(*train_locations);
        memcpy(z_values, aTrainData.back().data(), aTrainData.back().size() * sizeof(double));
        for(int i = 0; i< aTestData[0].size(); i++){
            test_locations->SetLocationX(*aTestData[0].data(), aTestData[0].size());
            test_locations->SetLocationY(*aTestData[1].data(), aTestData[1].size());
        }

        // Set common configurations.
        configurations.CheckKernelValue(aKernelName);
        configurations.ParseDistanceMetric(aDistanceMatrix);
        configurations.SetDenseTileSize(aDenseTileSize);
        configurations.SetLowTileSize(aLowTileSize);
        configurations.SetDimension(Configurations::CheckDimensionValue(aDimension));
        configurations.SetProblemSize(data->GetLocations()->GetSize());
        configurations.SetEstimatedTheta(aEstimatedTheta);
        configurations.SetInitialTheta(aTrueTheta);
        configurations.SetIsMLOEMMOM(TRUE);

        ExaGeoStat<double>::ExaGeoStatPrediction(configurations, data, z_values, train_locations, test_locations);

        // Take back ownership, to avoid deleting apData when the unique_ptr goes out of scope.
        auto temp_data = data.release();
        delete[] z_values;

        vector<double> mloe_mmom_values;
        mloe_mmom_values.push_back(Results::GetInstance()->GetMLOE());
        mloe_mmom_values.push_back(Results::GetInstance()->GetMMOM());
        return mloe_mmom_values;
    }


    vector<double>
    R_ExaGeoStatFisher(const string &aKernelName, const string &aDistanceMatrix, const vector<double> &aEstimatedTheta,
                       const int &aDenseTileSize, const int &aLowTileSize, const string &aDimension,
                       vector<vector<double>> &aTrainData,
                       vector<vector<double>> &aTestData) {

        Configurations configurations;
        configurations.SetComputation(EXACT_DENSE);

        auto data = make_unique<ExaGeoStatData<double>>(aTrainData[0].size(), configurations.GetDimension());
        auto *train_locations = new dataunits::Locations<double>(aTrainData[0].size(), configurations.GetDimension());
        auto *test_locations = new dataunits::Locations<double>(aTestData[0].size(), configurations.GetDimension());

        auto* z_values = new double[aTrainData.back().size()];
        for(int i = 0; i< aTrainData[0].size(); i++){
            train_locations->SetLocationX(*aTrainData[0].data(), aTrainData[0].size());
            train_locations->SetLocationY(*aTrainData[1].data(), aTrainData[1].size());
            if(aTrainData.size() > 3) {
                train_locations->SetLocationZ(*aTrainData[2].data(), aTrainData[2].size());
            }
        }
        data->SetLocations(*train_locations);
        memcpy(z_values, aTrainData.back().data(), aTrainData.back().size() * sizeof(double));
        for(int i = 0; i< aTestData[0].size(); i++){
            test_locations->SetLocationX(*aTestData[0].data(), aTestData[0].size());
            test_locations->SetLocationY(*aTestData[1].data(), aTestData[1].size());
        }

        // Set common configurations.
        configurations.CheckKernelValue(aKernelName);
        configurations.ParseDistanceMetric(aDistanceMatrix);
        configurations.SetProblemSize(data->GetLocations()->GetSize());
        configurations.SetDenseTileSize(aDenseTileSize);
        configurations.SetLowTileSize(aLowTileSize);
        configurations.SetEstimatedTheta(aEstimatedTheta);
        configurations.SetIsFisher(TRUE);

        ExaGeoStat<double>::ExaGeoStatPrediction(configurations, data, z_values, train_locations, test_locations);
        // Take back ownership, to avoid deleting apData when the unique_ptr goes out of scope.
        auto temp_data = data.release();
        delete[] z_values;
        return Results::GetInstance()->GetFisherMatrix();
    }

    vector<double> R_ExaGeoStatIDW(const string &aKernelName, const string &aDistanceMatrix,
                                   const vector<double> &aEstimatedTheta,
                                   const int &aDenseTileSize, const int &aLowTileSize,
                                   const string &aDimension,
                                   vector<vector<double>> &aTrainData,
                                   vector<vector<double>> &aTestData, vector<double> &aTestMeasurementsValues) {

        Configurations configurations;
        configurations.SetComputation(EXACT_DENSE);

        auto data = make_unique<ExaGeoStatData<double>>(aTrainData[0].size(), configurations.GetDimension());
        auto *train_locations = new dataunits::Locations<double>(aTrainData[0].size(), configurations.GetDimension());
        auto *test_locations = new dataunits::Locations<double>(aTestData[0].size(), configurations.GetDimension());

        // Allocate memory for z_values to hold elements from both sources
        auto* z_values = new double[aTrainData.back().size() + aTestMeasurementsValues.size()];

        for(int i = 0; i< aTrainData[0].size(); i++){
            train_locations->SetLocationX(*aTrainData[0].data(), aTrainData[0].size());
            train_locations->SetLocationY(*aTrainData[1].data(), aTrainData[1].size());
            if(aTrainData.size() > 3) {
                train_locations->SetLocationZ(*aTrainData[2].data(), aTrainData[2].size());
            }
        }
        // Copy data from the last element of aTrainData to the beginning of z_values
        memcpy(z_values, aTrainData.back().data(), aTrainData.back().size() * sizeof(double));
        // Calculate the starting position for the next part of the data in z_values
        auto* destination = z_values + aTrainData.back().size();
        // Copy data from aTestMeasurementsValues to the next part of z_values, after the previously copied data
        memcpy(destination, aTestMeasurementsValues.data(), aTestMeasurementsValues.size() * sizeof(double));

        for(int i = 0; i< aTestData[0].size(); i++){
            test_locations->SetLocationX(*aTestData[0].data(), aTestData[0].size());
            test_locations->SetLocationY(*aTestData[1].data(), aTestData[1].size());
        }

        // Set common configurations.
        configurations.CheckKernelValue(aKernelName);
        configurations.ParseDistanceMetric(aDistanceMatrix);
        configurations.SetDenseTileSize(aDenseTileSize);
        configurations.SetLowTileSize(aLowTileSize);
        configurations.SetDimension(Configurations::CheckDimensionValue(aDimension));
        configurations.SetProblemSize(data->GetLocations()->GetSize());
        configurations.SetEstimatedTheta(aEstimatedTheta);
        configurations.SetIsIDW(TRUE);

        ExaGeoStat<double>::ExaGeoStatPrediction(configurations, data, z_values, train_locations, test_locations);
        // Take back ownership, to avoid deleting apData when the unique_ptr goes out of scope.
        auto temp_data = data.release();
        delete[] z_values;

        return Results::GetInstance()->GetIDWError();
    }

    double *GetDataFromArguments(Nullable <NumericVector> aMeasurementsVector, Nullable <NumericVector> aLocationsX,
                                 Nullable <NumericVector> aLocationsY, Nullable <NumericVector> aLocationsZ,
                                 unique_ptr <ExaGeoStatData<double>> &aData, Configurations &aConfigurations,
                                 const string &aKernelName, const string &aDistanceMatrix, const int &aDenseTileSize,
                                 const int &aLowTileSize,
                                 const string &aDimension, const Computation &aComputation) {

        double *pMeasurementsVectorPtr = nullptr;
        if (aMeasurementsVector == nullptr || aMeasurementsVector.isNull()) {
            // This was the only way to pass C++ tests, if we checked for aMeasurementsVector != nullptr a seg fault occurs
        } else {
            NumericVector mat(aMeasurementsVector.get());
            pMeasurementsVectorPtr = &mat[0]; // Get the raw pointer to the matrix data

            if (!aLocationsX.isNull()) {
                double *pLocationsXPtr;
                NumericVector mat_location(aLocationsX.get());
                pLocationsXPtr = &mat_location[0]; // Get the raw pointer to the matrix data
                // Check for overflow
                if (mat_location.size() > numeric_limits<int>::max()) {
                    // Handle the overflow by logging an error
                    throw runtime_error("Warning: Size exceeds the limit of int type for Locations X.");
                } else {
                    // Safe to convert
                    int size = static_cast<int>(mat_location.size());
                    aData->GetLocations()->SetLocationX(*pLocationsXPtr, size);
                }
            }

            if (aLocationsY.isNotNull()) {
                double *pLocationsYPtr;
                NumericVector mat_location(aLocationsY.get());
                pLocationsYPtr = &mat_location[0]; // Get the raw pointer to the matrix data
                // Check for overflow
                if (mat_location.size() > numeric_limits<int>::max()) {
                    // Handle the overflow by logging an error
                    throw runtime_error("Warning: Size exceeds the limit of int type for Locations Y.");
                } else {
                    // Safe to convert
                    int size = static_cast<int>(mat_location.size());
                    aData->GetLocations()->SetLocationY(*pLocationsYPtr, size);
                }
            }

            if (aLocationsZ.isNotNull()) {
                double *pLocationsZPtr;
                NumericVector mat_location(aLocationsZ.get());
                pLocationsZPtr = &mat_location[0]; // Get the raw pointer to the matrix data
                // Check for overflow
                if (mat_location.size() > numeric_limits<int>::max()) {
                    // Handle the overflow by logging an error
                    throw runtime_error("Warning: Size exceeds the limit of int type for Locations Z.");
                } else {
                    // Safe to convert
                    int size = static_cast<int>(mat_location.size());
                    aData->GetLocations()->SetLocationZ(*pLocationsZPtr, size);
                }
            }
        }
        // Set common configurations.
        aConfigurations.SetComputation(aComputation);
        aConfigurations.CheckKernelValue(aKernelName);
        aConfigurations.ParseDistanceMetric(aDistanceMatrix);
        aConfigurations.SetDenseTileSize(aDenseTileSize);
        aConfigurations.SetLowTileSize(aLowTileSize);
        aConfigurations.SetDimension(Configurations::CheckDimensionValue(aDimension));
        aConfigurations.SetProblemSize(aData->GetLocations()->GetSize());
        return pMeasurementsVectorPtr;
    }
}