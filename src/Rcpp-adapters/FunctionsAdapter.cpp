
// Copyright (c) 2017-2024 King Abdullah University of Science and Technology,
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

    vector<vector<double>> R_GetLocations(ExaGeoStatData<double> *apData) {

        // Number of points per each dimension
        int length = apData->GetLocations()->GetSize();
        bool is3D = apData->GetLocations()->GetDimension() == Dimension3D;

        double *locationXArray = apData->GetLocations()->GetLocationX();
        double *locationYArray = apData->GetLocations()->GetLocationY();
        double *locationZArray = nullptr;
        if (is3D) {
            locationZArray = apData->GetLocations()->GetLocationZ();
        }

        vector<vector<double>> locations_matrix;
        for (int i = 0; i < length; ++i) {
            vector<double> point;
            point.push_back(locationXArray[i]);
            point.push_back(locationYArray[i]);
            if (is3D) {
                point.push_back(locationZArray[i]);
            }
            locations_matrix.push_back(point);
        }
        return locations_matrix;
    }

    NumericVector R_GetDescZValues(ExaGeoStatData<double> *apData, const string &aType) {

        DescriptorType descriptorType;
        void *pDescriptor;

        if (aType == "chameleon" || aType == "Chameleon") {
            descriptorType = CHAMELEON_DESCRIPTOR;
        } else if (aType == "hicma" || aType == "Hicma" || aType == "HICMA") {
#ifdef USE_HICMA
            descriptorType = HICMA_DESCRIPTOR;
#else
            throw runtime_error("Please enable HiCMA to use HiCMA descriptors.");
#endif
        } else {
            throw domain_error("Invalid type of descriptor, please use chameleon or hicma.");
        }

        // Obtain the pointer to the array of doubles
        double *data = apData->GetDescriptorData()->GetDescriptorMatrix(descriptorType, DESCRIPTOR_Z);
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
                          const string &aDimension, const int &aBand, const int &aMaxRank, SEXP apData,
                          Nullable <NumericVector> aMeasurementsVector, Nullable <NumericVector> aLocationsX,
                          Nullable <NumericVector> aLocationsY, Nullable <NumericVector> aLocationsZ) {

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
                                                              aLocationsZ, data, configurations, aKernelName,
                                                              aDistanceMatrix, aDenseTileSize, aLowTileSize, aDimension,
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
                                           const vector<double> &aEstimatedTheta, const int &aDenseTileSize,
                                           const int &aLowTileSize, const string &aDimension,
                                           vector <vector<double>> &aTrainData, vector <vector<double>> &aTestData) {

        Configurations configurations;
        configurations.SetIsMSPE(TRUE);
        configurations.SetEstimatedTheta(aEstimatedTheta);
        vector<double> empty_vector;
        PredictionSetupHelper(configurations, aKernelName, aDistanceMatrix, aDenseTileSize, aLowTileSize, aDimension,
                              aTrainData, aTestData, aEstimatedTheta, empty_vector);
        return Results::GetInstance()->GetPredictedMissedValues();
    }

    vector<double> R_ExaGeoStatMLOE_MMOM(const string &aKernelName, const string &aDistanceMatrix,
                                         const vector<double> &aEstimatedTheta, const vector<double> &aTrueTheta,
                                         const int &aDenseTileSize, const int &aLowTileSize, const string &aDimension,
                                         vector <vector<double>> &aTrainData, vector <vector<double>> &aTestData) {

        Configurations configurations;
        configurations.SetEstimatedTheta(aEstimatedTheta);
        configurations.SetInitialTheta(aTrueTheta);
        configurations.SetIsMLOEMMOM(TRUE);

        vector<double> empty_vector;
        PredictionSetupHelper(configurations, aKernelName, aDistanceMatrix, aDenseTileSize, aLowTileSize, aDimension,
                              aTrainData, aTestData, aEstimatedTheta, empty_vector);

        vector<double> mloe_mmom_values;
        mloe_mmom_values.push_back(Results::GetInstance()->GetMLOE());
        mloe_mmom_values.push_back(Results::GetInstance()->GetMMOM());
        return mloe_mmom_values;
    }

    vector<double>
    R_ExaGeoStatFisher(const string &aKernelName, const string &aDistanceMatrix, const vector<double> &aEstimatedTheta,
                       const int &aDenseTileSize, const int &aLowTileSize, const string &aDimension,
                       vector <vector<double>> &aTrainData, vector <vector<double>> &aTestData) {

        Configurations configurations;
        configurations.SetIsFisher(TRUE);

        vector<double> empty_vector;
        PredictionSetupHelper(configurations, aKernelName, aDistanceMatrix, aDenseTileSize, aLowTileSize, aDimension,
                              aTrainData, aTestData, aEstimatedTheta, empty_vector);

        return Results::GetInstance()->GetFisherMatrix();
    }

    vector<double>
    R_ExaGeoStatIDW(const string &aKernelName, const string &aDistanceMatrix, const vector<double> &aEstimatedTheta,
                    const int &aDenseTileSize, const int &aLowTileSize, const string &aDimension,
                    vector <vector<double>> &aTrainData, vector <vector<double>> &aTestData,
                    vector<double> &aTestMeasurementsValues) {

        Configurations configurations;
        configurations.SetIsIDW(TRUE);

        PredictionSetupHelper(configurations, aKernelName, aDistanceMatrix, aDenseTileSize, aLowTileSize, aDimension,
                              aTrainData, aTestData, aEstimatedTheta, aTestMeasurementsValues);

        return Results::GetInstance()->GetIDWError();
    }

    double *GetDataFromArguments(Nullable <NumericVector> aMeasurementsVector, Nullable <NumericVector> aLocationsX,
                                 Nullable <NumericVector> aLocationsY, Nullable <NumericVector> aLocationsZ,
                                 unique_ptr <ExaGeoStatData<double>> &aData, Configurations &aConfigurations,
                                 const string &aKernelName, const string &aDistanceMatrix, const int &aDenseTileSize,
                                 const int &aLowTileSize, const string &aDimension, const Computation &aComputation) {

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

    void ValidateDataDimensions(const vector <vector<double>> &aData, const string &aDataType) {
        if (aData.empty() || aData.front().empty()) {
            throw runtime_error("No " + aDataType + " data provided.");
        }

        size_t expectedSize = aData.front().size();
        for (const auto &vec: aData) {
            if (vec.size() != expectedSize) {
                throw runtime_error("Inconsistent dimensions in " + aDataType + " data.");
            }
        }
    }

    void
    PredictionSetupHelper(Configurations &aConfigurations, const string &aKernelName, const string &aDistanceMatrix,
                          const int &aDenseTileSize, const int &aLowTileSize, const string &aDimension,
                          vector <vector<double>> &aTrainData, vector <vector<double>> &aTestData,
                          const vector<double> &aEstimatedTheta, const vector<double> &aTestMeasurementsValues) {

        aConfigurations.SetComputation(EXACT_DENSE);

        ValidateDataDimensions(aTrainData, "train");
        ValidateDataDimensions(aTestData, "test");

        size_t potential_train_data_size = aTrainData[0].size();
        if (potential_train_data_size > static_cast<size_t>(numeric_limits<int>::max())) {
            throw runtime_error("Train data size exceeds the maximum size for int type.");
        }
        int train_data_size = static_cast<int>(potential_train_data_size);

        size_t potential_test_data_size = aTestData[0].size();
        if (potential_test_data_size > static_cast<size_t>(numeric_limits<int>::max())) {
            throw runtime_error("Test data size exceeds the maximum size for int type.");
        }
        int test_data_size = static_cast<int>(potential_test_data_size);

        auto data = make_unique<ExaGeoStatData<double>>(train_data_size, aConfigurations.GetDimension());
        dataunits::Locations<double> train_locations(train_data_size, aConfigurations.GetDimension());
        dataunits::Locations<double> test_locations(test_data_size, aConfigurations.GetDimension());

        // Allocate memory for z_values to hold elements from both sources
        auto *z_values = new double[aTrainData.back().size() + aTestMeasurementsValues.size()];

        for (int i = 0; i < train_data_size; i++) {
            train_locations.SetLocationX(*aTrainData[0].data(), train_data_size);
            train_locations.SetLocationY(*aTrainData[1].data(), train_data_size);
            data->GetLocations()->SetLocationX(*aTrainData[0].data(), train_data_size);
            data->GetLocations()->SetLocationY(*aTrainData[1].data(), train_data_size);
            if (aTrainData.size() > 3) {
                train_locations.SetLocationZ(*aTrainData[2].data(), train_data_size);
                data->GetLocations()->SetLocationZ(*aTrainData[2].data(), train_data_size);
            }
        }
        memcpy(z_values, aTrainData.back().data(), aTrainData.back().size() * sizeof(double));
        // Calculate the starting position for the next part of the data in z_values
        auto *destination = z_values + aTrainData.back().size();
        // Copy data from aTestMeasurementsValues to the next part of z_values, after the previously copied data
        memcpy(destination, aTestMeasurementsValues.data(), aTestMeasurementsValues.size() * sizeof(double));

        for (int i = 0; i < test_data_size; i++) {
            test_locations.SetLocationX(*aTestData[0].data(), test_data_size);
            test_locations.SetLocationY(*aTestData[1].data(), test_data_size);
            if (aTestData.size() > 2) {
                test_locations.SetLocationZ(*aTestData[2].data(), test_data_size);
            }
        }

        // Set common configurations.
        aConfigurations.CheckKernelValue(aKernelName);
        aConfigurations.ParseDistanceMetric(aDistanceMatrix);
        aConfigurations.SetDenseTileSize(aDenseTileSize);
        aConfigurations.SetLowTileSize(aLowTileSize);
        aConfigurations.SetDimension(Configurations::CheckDimensionValue(aDimension));
        aConfigurations.SetProblemSize(data->GetLocations()->GetSize());
        aConfigurations.SetEstimatedTheta(aEstimatedTheta);

        // Temporarily release ownership to pass to the function.
        ExaGeoStatData<double> *tempData = data.release();
        ExaGeoStat<double>::ExaGeoStatPrediction(aConfigurations,
                                                 reinterpret_cast<unique_ptr<ExaGeoStatData<double>> &>(tempData),
                                                 z_values, &train_locations, &test_locations);
        // Take back ownership
        data.reset(tempData);

        delete[] z_values;
    }
}