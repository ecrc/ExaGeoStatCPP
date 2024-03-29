
// Copyright (c) 2017-2023 King Abdullah University of Science and Technology,
// All rights reserved.
// ExaGeoStat is a software package, provided by King Abdullah University of Science and Technology (KAUST).

/**
 * @file FunctionAdapter.hpp
 * @brief Header file for function adapters in the ExaGeoStat software.
 * @details  It provides declarations for functions that adapt and initialize statistical models or algorithms based on user-defined configurations.
 * @version 1.1.0
 * @author Mahmoud ElKarargy
 * @date 2024-01-29
**/

#ifndef EXAGEOSTATCPP_FUNCTIONSADAPTER_HPP
#define EXAGEOSTATCPP_FUNCTIONSADAPTER_HPP

#include <Rcpp.h>
#include <iostream>
#include <string>
#include <vector>

#include <configurations/Configurations.hpp>
#include <data-units/ExaGeoStatData.hpp>

namespace exageostat::adapters {

    /**
     * @brief Initializes and configures the R arguments for ExaGeoStat computations.
     * @details This function prepares the necessary configurations required by ExaGeoStat to perform statistical computations.
     * It includes setting up problem sizes, computational kernels, grid configurations, and other parameters essential
     * for the execution of the ExaGeoStat algorithms.
     * @param[in] aProblemSize The size of the problem to be solved.
     * @param[in] aKernelName The name of the computational kernel to be used.
     * @param[in] aTileSize A vector specifying the size of each tile in the computation.
     * @param[in] aP_QGrid A vector defining the P x Q process grid.
     * @param[in] aTimeSlot The time slot allocated for the computation.
     * @param[in] aComputation The type of computation to be performed (e.g., estimation, prediction).
     * @param[in] aPrecision The precision (e.g., single, double) of the computation.
     * @param[in] aCoresGPUsNumber A vector specifying the number of CPU cores or GPUs to be used.
     * @param[in] aBand The bandwidth of the problem, relevant for band-limited computations.
     * @param[in] aMaxRank The maximum rank for low-rank approximations.
     * @param[in] aInitialTheta A vector of initial values for the model parameters (theta).
     * @param[in] aLowerUpperBounds A 2D vector specifying the lower and upper bounds for model parameters.
     * @param[in] aEstimatedTheta A vector of estimated values for the model parameters after computation.
     * @param[in] aVerbose A string indicating the verbosity level of the output.
     * @param[in] aDimension The dimensionality of the problem (e.g., 2D, 3D).
     * @param[in] aMaxMleIterations The maximum number of iterations for the Maximum Likelihood Estimation (MLE) algorithm.
     * @param[in] aTolerance The tolerance threshold for convergence in iterative algorithms.
     * @param[in] aPrediction A vector indicating prediction locations or settings.
     * @return A pointer to a Configurations object containing the initialized settings.
     *
     */

    /**
     * @brief Retrieves X coordinates of locations from ExaGeoStat data.
     * @details Extracts and returns the X coordinates of geographical or spatial locations stored in an ExaGeoStatData object, facilitating data manipulation and analysis within the ExaGeoStat framework.
     * @param[in] apData Pointer to ExaGeoStatData object containing the spatial data.
     * @return Numeric vector of X coordinates.
     */
    Rcpp::NumericVector R_GetLocationX(ExaGeoStatData<double> *apData);

    /**
     * @brief Retrieves Y coordinates of locations from ExaGeoStat data.
     * @details Extracts and returns the Y coordinates of geographical or spatial locations stored in an ExaGeoStatData object, supporting various spatial data analyses and operations within the ExaGeoStat software.
     * @param[in] apData Pointer to ExaGeoStatData object containing the spatial data.
     * @return Numeric vector of Y coordinates.
     */
    Rcpp::NumericVector R_GetLocationY(ExaGeoStatData<double> *apData);

    /**
     * @brief Retrieves Z coordinates of locations from ExaGeoStat data.
     * @details Extracts and returns the Z coordinates (elevation or depth) of spatial locations stored in an ExaGeoStatData object, enhancing three-dimensional spatial analysis capabilities within the ExaGeoStat framework.
     * @param[in] apData Pointer to ExaGeoStatData object containing the spatial data.
     * @return Numeric vector of Z coordinates.
     */
    Rcpp::NumericVector R_GetLocationZ(ExaGeoStatData<double> *apData);

    /**
     * @brief Retrieves descriptive Z values from ExaGeoStat data based on type.
     * @details Extracts and returns Z values from an ExaGeoStatData object, aiding in targeted spatial data analysis and visualization within ExaGeoStat.
     * @param[in] apData Pointer to ExaGeoStatData object containing the spatial data.
     * @param[in] aType String specifying the type of descriptor value to retrieve (e.g., "Chameleon", "HiCMA").
     * @return Numeric vector of descriptive Z values.
     */
    Rcpp::NumericVector R_GetDescZValues(ExaGeoStatData<double> *apData, const std::string &aType);

    /**
     * @brief Function to load ExaGeoStat data.
     * @details This function loads data into an ExaGeoStatData object using the provided configuration and computational settings.
     * It is designed to initialize the data structure necessary for subsequent statistical model operations within the ExaGeoStat framework.
     * @param[in] apData A pointer to an ExaGeoStatData object where the loaded data will be stored.
     * @return A pointer to an ExaGeoStatData object containing the loaded data.
     *
     */
    ExaGeoStatData<double> *
    R_ExaGeoStatLoadData(const std::string &aKernelName, const std::vector<double> &aInitialTheta,
                         const std::string &aDistanceMatrix, const int &aProblemSize, const int &aSeed,
                         const int &aDenseTileSize, const int &aLowTileSize, const std::string &aDimension,
                         const std::string &aLogPath, const std::string &aDataPath,
                         const std::string &aRecoveryFilePath, const std::string &aObservationsFilePath);

    /**
    * @brief Function to model ExaGeoStat data.
    * @details This function applies the model configurations specified in apConfigurations to the data stored in apData.
    * It prepares the ExaGeoStatData object for statistical analysis and prediction tasks, adjusting the internal data representations and parameters as needed.
    * @param[in] apData A pointer to an ExaGeoStatData object where the loaded data is stored and will be modeled.
    * @param[in] aMeasurementsVector An optional Rcpp::Nullable object containing a vector of measurements. If provided, it is used to enhance the data modeling process.
    * @param[in] aLocationsX An optional Rcpp::Nullable object containing a vector of X coordinates for the locations. If provided, it is used to enhance the data modeling process.
    * @param[in] aLocationsY An optional Rcpp::Nullable object containing a vector of Y coordinates for the locations. If provided, it is used to enhance the data modeling process.
    * @param[in] aLocationsZ An optional Rcpp::Nullable object containing a vector of Z coordinates for the locations. If provided, it is used to enhance the data modeling process.
    * @return A pointer to an ExaGeoStatData object containing the modeled data, ready for statistical analysis and prediction.
    */
    std::vector<double> R_ExaGeoStatModelData(const std::string &aComputation, const std::string &aKernelName,
                                              const std::string &aDistanceMatrix,
                                              const std::vector<double> &aLowerBound,
                                              const std::vector<double> &aUpperBound, const int &aTolerance,
                                              const int &aMleIterations, const int &aDenseTileSize,
                                              const int &aLowTileSize, const std::string &aDimension, const int &aBand, const int &aMaxRank,
                                              SEXP apData,
                                              Rcpp::Nullable<Rcpp::NumericVector> aMeasurementsVector = R_NilValue,
                                              Rcpp::Nullable<Rcpp::NumericVector> aLocationsX = R_NilValue,
                                              Rcpp::Nullable<Rcpp::NumericVector> aLocationsY = R_NilValue,
                                              Rcpp::Nullable<Rcpp::NumericVector> aLocationsZ = R_NilValue);

    /**
     * @brief Function to predict using ExaGeoStat data.
     * @details This function performs predictions based on the modeled ExaGeoStatData object.
     * It utilizes the configuration and computational settings defined by apConfigurations, respectively, to execute prediction algorithms on the data stored in apData.
     * @param[in] apConfigurations A pointer to a Configurations object containing the computational settings.
     * @param[in] apData A pointer to an ExaGeoStatData object where the loaded data is stored and will be modeled.
     * @param[in] aMeasurementsVector An optional Rcpp::Nullable object containing a vector of measurements. If provided, it is used to enhance the data modeling process.
     * @param[in] aLocationsX An optional Rcpp::Nullable object containing a vector of X coordinates for the locations. If provided, it is used to enhance the data modeling process.
     * @param[in] aLocationsY An optional Rcpp::Nullable object containing a vector of Y coordinates for the locations. If provided, it is used to enhance the data modeling process.
     * @param[in] aLocationsZ An optional Rcpp::Nullable object containing a vector of Z coordinates for the locations. If provided, it is used to enhance the data modeling process.
     * @return A pointer to an ExaGeoStatData object containing the modeled data, ready for statistical analysis and prediction.
     */
    std::vector<double> R_ExaGeoStatPredictData(const std::string &aKernelName, const std::string &aDistanceMatrix,
                                 const std::vector<double> &aEstimatedTheta,
                                 const int &aDenseTileSize, const int &aLowTileSize,
                                 const std::string &aDimension, std::vector<std::vector<double>> &aTrainData,
                                 std::vector<std::vector<double>> &aTestData);

    std::vector<double> R_ExaGeoStatMLOE_MMOM(const std::string &aKernelName, const std::string &aDistanceMatrix,
                                                  const std::vector<double> &aEstimatedTheta,
                                                  const std::vector<double> &aTrueTheta,
                                                  const int &aDenseTileSize, const int &aLowTileSize,
                                                  const std::string &aDimension, std::vector<std::vector<double>> &aTrainData,
                                              std::vector<std::vector<double>> &aTestData);

    std::vector<double> R_ExaGeoStatFisher(const std::string &aKernelName, const std::string &aDistanceMatrix,
                                               const std::vector<double> &aEstimatedTheta,
                                               const int &aDenseTileSize,
                                               const int &aLowTileSize, const std::string &aDimension,
                                           std::vector<std::vector<double>> &aTrainData,
                                           std::vector<std::vector<double>> &aTestData);

    std::vector<double> R_ExaGeoStatIDW(const std::string &aKernelName, const std::string &aDistanceMatrix,
                                            const std::vector<double> &aEstimatedTheta,
                                            const int &aDenseTileSize, const int &aLowTileSize,
                                            const std::string &aDimension,
                                        std::vector<std::vector<double>> &aTrainData,
                                        std::vector<std::vector<double>> &aTestData, std::vector<double> &aTestMeasurementsValues);

    double *GetDataFromArguments(Rcpp::Nullable <Rcpp::NumericVector> aMeasurementsVector, Rcpp::Nullable <Rcpp::NumericVector> aLocationsX,Rcpp::Nullable <Rcpp::NumericVector> aLocationsY, Rcpp::Nullable <Rcpp::NumericVector> aLocationsZ, std::unique_ptr<ExaGeoStatData<double>> &aData, configurations::Configurations &aConfigurations,
                                 const std::string &aKernelName, const std::string &aDistanceMatrix,
                                 const int &aDenseTileSize, const int &aLowTileSize, const std::string &aDimension, const common::Computation &aComputation);

}
#endif //EXAGEOSTATCPP_FUNCTIONSADAPTER_HPP
