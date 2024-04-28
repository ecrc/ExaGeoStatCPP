
// Copyright (c) 2017-2024 King Abdullah University of Science and Technology,
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

#include <configurations/Configurations.hpp>
#include <data-units/ExaGeoStatData.hpp>

namespace exageostat::adapters {

    /**
     * @brief Retrieves locations from ExaGeoStat data.
     * @details Extracts and returns the locations stored in an ExaGeoStatData object,
     * @param[in] apData Pointer to ExaGeoStatData object containing the spatial data.
     * @return vector of locations coordinates.
     *
     */
    std::vector<std::vector<double>> R_GetLocations(ExaGeoStatData<double> *apData);

    /**
     * @brief Retrieves descriptive Z values from ExaGeoStat data based on type.
     * @details Extracts and returns Z values from an ExaGeoStatData object, aiding in targeted spatial data analysis and visualization within ExaGeoStat.
     * @param[in] apData Pointer to ExaGeoStatData object containing the spatial data.
     * @param[in] aType String specifying the type of descriptor value to retrieve (e.g., "Chameleon", "HiCMA").
     * @return Numeric vector of descriptive Z values.
     *
     */
    Rcpp::NumericVector R_GetDescZValues(ExaGeoStatData<double> *apData, const std::string &aType);

    /**
     * @brief Function to load ExaGeoStat data.
     * @details This function loads data into an ExaGeoStatData object using the provided configuration and computational settings.
     * It is designed to initialize the data structure necessary for subsequent statistical model operations within the ExaGeoStat framework.
     * @param[in] aKernelName Name of the computational kernel to be utilized.
     * @param[in] aInitialTheta Initial parameter values for the statistical model.
     * @param[in] aDistanceMatrix Type of distance matrix to be used ("euclidean", "manhattan", etc.).
     * @param[in] aProblemSize Size of the problem or dataset.
     * @param[in] aSeed Seed for random number generation, ensuring reproducibility.
     * @param[in] aDenseTileSize Size of the tile for dense computations.
     * @param[in] aLowTileSize Size of the tile for low-rank computations.
     * @param[in] aDimension Dimensionality of the problem ("2D" for two dimensions, "3D" for three dimensions).
     * @param[in] aLogPath Path to the log file where execution details will be stored.
     * @param[in] aDataPath Path to the data file containing spatial observations.
     * @param[in] aRecoveryFilePath Path for saving intermediate computation states, aiding in recovery from interruptions.
     * @param[in] aObservationsFilePath Path to the file containing observation data.
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
     * @brief Models ExaGeoStat data using specified arguments.
     * @details Applies statistical modeling to ExaGeoStatData based on the provided configurations.
     * This function is essential for preparing the data for in-depth statistical analysis and predictions,
     * optimizing internal representations and parameters for the modeling process.
     * @param[in] aComputation Computational method to be used.
     * @param[in] aKernelName Name of the kernel for computations.
     * @param[in] aDistanceMatrix Type of distance matrix ("euclidean", "manhattan", etc.).
     * @param[in] aLowerBound Lower bound for optimization parameters.
     * @param[in] aUpperBound Upper bound for optimization parameters.
     * @param[in] aTolerance Tolerance level for the optimization algorithm.
     * @param[in] aMleIterations Maximum number of iterations for the Maximum Likelihood Estimation (MLE) algorithm.
     * @param[in] aDenseTileSize Tile size for dense matrix computations.
     * @param[in] aLowTileSize Tile size for low-rank approximations.
     * @param[in] aDimension Dimensionality of the problem ("2D" or "3D").
     * @param[in] aBand Bandwidth for band matrices, applicable in certain computational kernels.
     * @param[in] aMaxRank Maximum rank for low-rank approximations.
     * @param[in] apData Pointer to ExaGeoStatData object to be modeled.
     * @param[in] aMeasurementsVector Optional vector of measurements to enhance modeling, can be nullable.
     * @param[in] aLocationsX Optional vector of X coordinates for locations, can be nullable.
     * @param[in] aLocationsY Optional vector of Y coordinates for locations, can be nullable.
     * @param[in] aLocationsZ Optional vector of Z coordinates for locations, can be nullable.
     * @return Vector of doubles represents the modeled theta.
     *
     */
    std::vector<double> R_ExaGeoStatModelData(const std::string &aComputation, const std::string &aKernelName,
                                              const std::string &aDistanceMatrix,
                                              const std::vector<double> &aLowerBound,
                                              const std::vector<double> &aUpperBound, const int &aTolerance,
                                              const int &aMleIterations, const int &aDenseTileSize,
                                              const int &aLowTileSize, const std::string &aDimension, const int &aBand,
                                              const int &aMaxRank, SEXP apData,
                                              Rcpp::Nullable<Rcpp::NumericVector> aMeasurementsVector = R_NilValue,
                                              Rcpp::Nullable<Rcpp::NumericVector> aLocationsX = R_NilValue,
                                              Rcpp::Nullable<Rcpp::NumericVector> aLocationsY = R_NilValue,
                                              Rcpp::Nullable<Rcpp::NumericVector> aLocationsZ = R_NilValue);

    /**
     * @brief Predicts outcomes using ExaGeoStat data and configurations.
     * @details Utilizes a modeled ExaGeoStatData object to perform predictions, leveraging specified computational settings and statistical models.
     * This function is integral for generating spatial predictions based on the data and models within the ExaGeoStat framework.
     * @param[in] aKernelName Name of the kernel used for prediction computations.
     * @param[in] aDistanceMatrix Type of distance matrix used ("euclidean", "manhattan", etc.).
     * @param[in] aEstimatedTheta Vector of estimated parameters from the model.
     * @param[in] aDenseTileSize Tile size for dense matrix operations.
     * @param[in] aLowTileSize Tile size for low-rank matrix operations.
     * @param[in] aDimension Dimensionality of the spatial data ("2D" or "3D").
     * @param[in] aTrainData Training data set used for predictions.
     * @param[in] aTestData Test data set for which predictions are made.
     * @return Vector of predicted values based on the test data.
     *
     */
    std::vector<double> R_ExaGeoStatPredictData(const std::string &aKernelName, const std::string &aDistanceMatrix,
                                                const std::vector<double> &aEstimatedTheta, const int &aDenseTileSize,
                                                const int &aLowTileSize, const std::string &aDimension,
                                                std::vector<std::vector<double>> &aTrainData,
                                                std::vector<std::vector<double>> &aTestData);

    /**
     * @brief Calculates the Mean Logarithmic Error (MLOE) and the Mean Measure of Model Output (MMOM) for ExaGeoStat predictions.
     * @details Assesses the accuracy of spatial predictions made by the ExaGeoStat framework by computing the MLOE and MMOM,
     * which provide insights into the predictive performance and uncertainty of the models.
     * @param[in] aKernelName Kernel used for the prediction computations.
     * @param[in] aDistanceMatrix Type of distance matrix ("euclidean", "manhattan", etc.).
     * @param[in] aEstimatedTheta Vector of estimated parameters from the model.
     * @param[in] aTrueTheta Vector of true parameter values for validation.
     * @param[in] aDenseTileSize Tile size for dense matrix operations.
     * @param[in] aLowTileSize Tile size for low-rank matrix operations.
     * @param[in] aDimension Dimensionality of the spatial data ("2D" or "3D").
     * @param[in] aTrainData Training data set used in the model.
     * @param[in] aTestData Test data set used for validation.
     * @return Vector containing the calculated MLOE and MMOM values.
     *
     */
    std::vector<double> R_ExaGeoStatMLOE_MMOM(const std::string &aKernelName, const std::string &aDistanceMatrix,
                                              const std::vector<double> &aEstimatedTheta,
                                              const std::vector<double> &aTrueTheta, const int &aDenseTileSize,
                                              const int &aLowTileSize, const std::string &aDimension,
                                              std::vector<std::vector<double>> &aTrainData,
                                              std::vector<std::vector<double>> &aTestData);

    /**
     * @brief Computes the Fisher information matrix for ExaGeoStat models.
     * @details Utilizes the estimated parameters and the Fisher information matrix to evaluate the information content and
     * parameter uncertainties within the ExaGeoStat framework, contributing to the understanding of model reliability and sensitivity.
     * @param[in] aKernelName Kernel used for computations.
     * @param[in] aDistanceMatrix Type of distance matrix ("euclidean", "manhattan", etc.).
     * @param[in] aEstimatedTheta Vector of estimated parameters from the model.
     * @param[in] aDenseTileSize Tile size for dense matrix operations.
     * @param[in] aLowTileSize Tile size for low-rank matrix operations.
     * @param[in] aDimension Dimensionality of the spatial data ("2D" or "3D").
     * @param[in] aTrainData Training data set used in the model.
     * @param[in] aTestData Test data set used for validation.
     * @return Vector represents the Fisher information matrix.
     *
     */
    std::vector<double> R_ExaGeoStatFisher(const std::string &aKernelName, const std::string &aDistanceMatrix,
                                           const std::vector<double> &aEstimatedTheta, const int &aDenseTileSize,
                                           const int &aLowTileSize, const std::string &aDimension,
                                           std::vector<std::vector<double>> &aTrainData,
                                           std::vector<std::vector<double>> &aTestData);

    /**
     * @brief Applies Inverse Distance Weighting (IDW) for spatial interpolation using ExaGeoStat data.
     * @details Implements the IDW interpolation method to estimate spatial variables at unsampled locations based on the distances
     * and values of nearby sampled points within the ExaGeoStat framework, enhancing spatial prediction capabilities.
     * @param[in] aKernelName Kernel used for IDW computations.
     * @param[in] aDistanceMatrix Type of distance matrix ("euclidean", "manhattan", etc.).
     * @param[in] aEstimatedTheta Vector of parameters, typically used for weighting in IDW.
     * @param[in] aDenseTileSize Tile size for dense matrix operations.
     * @param[in] aLowTileSize Tile size for low-rank matrix operations.
     * @param[in] aDimension Dimensionality of the spatial data ("2D" or "3D").
     * @param[in] aTrainData Training data set providing sampled locations and values.
     * @param[in] aTestData Test data set providing unsampled locations for which values are interpolated.
     * @param[in] aTestMeasurementsValues Vector of measured values at the test locations, used as reference in some IDW implementations.
     * @return Vector of interpolated values at the test locations.
     *
     */
    std::vector<double> R_ExaGeoStatIDW(const std::string &aKernelName, const std::string &aDistanceMatrix,
                                        const std::vector<double> &aEstimatedTheta, const int &aDenseTileSize,
                                        const int &aLowTileSize, const std::string &aDimension,
                                        std::vector<std::vector<double>> &aTrainData,
                                        std::vector<std::vector<double>> &aTestData,
                                        std::vector<double> &aTestMeasurementsValues);

    /**
     * @brief Extracts and prepares data from given arguments for ExaGeoStat operations.
     * @details This function is designed to parse and prepare spatial and measurement data from provided arguments,
     * making it suitable for processing within the ExaGeoStat framework. It handles optional data vectors for measurements
     * and locations (X, Y, Z coordinates), and configures an ExaGeoStatData object based on these inputs along with other
     * computational and configuration parameters.
     * @param[in] aMeasurementsVector vector of measurements to enhance modeling, can be nullable.
     * @param[in] aLocationsX vector of X coordinates for locations, can be nullable.
     * @param[in] aLocationsY vector of Y coordinates for locations, can be nullable.
     * @param[in] aLocationsZ vector of Z coordinates for locations, can be nullable.
     * @param[in] aData Pointer to ExaGeoStatData object to be modeled.
     * @param[in] aConfigurations Configuration settings specifying computational details such as the kernel type, matrix storage format, etc.
     * @param[in] aKernelName Name of the kernel for computations.
     * @param[in] aDistanceMatrix Type of distance matrix ("euclidean", "manhattan", etc.).
     * @param[in] aDenseTileSize Tile size for dense matrix computations.
     * @param[in] aLowTileSize Tile size for low-rank approximations.
     * @param[in] aDimension Dimensionality of the problem ("2D" or "3D").
     * @param[in] aComputation Computational method to be used.
     * @return Pointer to a double array containing the prepared data, ready for use in ExaGeoStat operations.
     *
     */
    double *GetDataFromArguments(Rcpp::Nullable<Rcpp::NumericVector> aMeasurementsVector,
                                 Rcpp::Nullable<Rcpp::NumericVector> aLocationsX,
                                 Rcpp::Nullable<Rcpp::NumericVector> aLocationsY,
                                 Rcpp::Nullable<Rcpp::NumericVector> aLocationsZ,
                                 std::unique_ptr<ExaGeoStatData<double>> &aData,
                                 configurations::Configurations &aConfigurations, const std::string &aKernelName,
                                 const std::string &aDistanceMatrix, const int &aDenseTileSize, const int &aLowTileSize,
                                 const std::string &aDimension, const common::Computation &aComputation);

    /**
     * @brief Validates the dimensions of input data.
     * @details This function checks the dimensions of the provided data vectors to ensure they meet the expected format and size requirements for a given data type.
     * It's used to verify that the data structures passed into algorithms or processes are correctly formatted, preventing errors or inconsistencies in data processing.
     * @param aData A constant reference to a vector of vectors containing the data to be validated.
     * @param aDataType A string describing the type of data being validated, which influences the expected dimensions and format of the data.
     * @return void
     *
     */
    void ValidateDataDimensions(const std::vector<std::vector<double>> &aData, const std::string &aDataType);

    /**
    * @brief Sets up the prediction environment.
    * @details This function prepares the necessary configurations and data structures for making predictions. It involves setting up various parameters, including kernel names, distance matrices, tile sizes, dimensions, and training/test data. The function is crucial for initializing the prediction process with the appropriate settings and data.
    * @param aConfigurations Reference to a Configurations object containing various prediction and algorithm configurations.
    * @param aKernelName Name of the kernel to be used in predictions.
    * @param aDistanceMatrix String representation of the distance matrix to be used.
    * @param aDenseTileSize Size of the dense tiles in the matrix.
    * @param aLowTileSize Size of the low-resolution tiles in the matrix.
    * @param aDimension String representation of the dimensionality of the data.
    * @param aTrainData Reference to a vector of vectors containing the training data.
    * @param aTestData Reference to a vector of vectors containing the test data.
    * @param aEstimatedTheta Vector containing estimated theta values for the model.
    * @param aTestMeasurementsValues Vector containing the test measurement values.
    * @return Pointer to a double array containing the prepared data, ready for use in ExaGeoStat operations.
    *
    */
    void PredictionSetupHelper(configurations::Configurations &aConfigurations, const std::string &aKernelName,
                               const std::string &aDistanceMatrix, const int &aDenseTileSize, const int &aLowTileSize,
                               const std::string &aDimension, std::vector<std::vector<double>> &aTrainData,
                               std::vector<std::vector<double>> &aTestData,
                               const std::vector<double> &aEstimatedTheta,
                               const std::vector<double> &aTestMeasurementsValues);

}
#endif //EXAGEOSTATCPP_FUNCTIONSADAPTER_HPP
