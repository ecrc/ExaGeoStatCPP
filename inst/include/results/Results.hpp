
// Copyright (c) 2017-2024 King Abdullah University of Science and Technology,
// All rights reserved.
// ExaGeoStat is a software package, provided by King Abdullah University of Science and Technology (KAUST).

/**
 * @file Results.hpp
 * @brief Defines the Results class for storing and accessing result data.
 * @version 1.1.0
 * @author Mahmoud ElKarargy
 * @date 2023-09-14
**/

#ifndef EXAGEOSTATCPP_RESULTS_HPP
#define EXAGEOSTATCPP_RESULTS_HPP

#include <iostream>
#include <vector>
#include <map>

namespace exageostat::results {

    class Results {

    public:

        /**
         * @brief Get a pointer to the singleton instance of the Results class.
         * @return A pointer to the instance of the Results class.
         *
         */
        static Results *GetInstance();

        /**
         * @brief Set whether the dataset is synthetic or not.
         * @param aIsSynthetic Boolean indicating if the dataset is synthetic.
         * @param aKey Custom dictionary key (optional).
         */
        void SetIsSynthetic(bool aIsSynthetic, const std::string &aKey = "");

        /**
         * @brief Set the number of generated locations.
         * @param aNumLocations Integer representing the number of locations.
         * @param aKey Custom dictionary key (optional).
         */
        void SetGeneratedLocationsNumber(int aNumLocations, const std::string &aKey = "");

        /**
         * @brief Enable or disable logging.
         * @param aIsLogger Boolean indicating if logging is enabled.
         * @param aKey Custom dictionary key (optional).
         */
        void SetIsLogger(bool aIsLogger, const std::string &aKey = "");

        /**
         * @brief Set the logger's file path.
         * @param aLoggerPath String representing the logger file path.
         * @param aKey Custom dictionary key (optional).
         */
        void SetLoggerPath(const std::string &aLoggerPath, const std::string &aKey = "");

        /**
         * @brief Set the total data generation execution time.
         * @param aTime Double representing the execution time.
         * @param aKey Custom dictionary key (optional).
         */
        void SetTotalDataGenerationExecutionTime(double aTime, const std::string &aKey = "");

        /**
         * @brief Set the total data generation FLOPS.
         * @param aFlops Double representing the FLOPS.
         * @param aKey Custom dictionary key (optional).
         */
        void SetTotalDataGenerationFlops(double aFlops, const std::string &aKey = "");

        /**
         * @brief Set the log-likelihood value.
         * @param aLogLikValue Double representing the log-likelihood value.
         * @param aKey Custom dictionary key (optional).
         */
        void SetLogLikValue(double aLogLikValue, const std::string &aKey = "");

        /**
         * @brief Set the number of MLE iterations.
         * @param aIterationsNumber Integer representing the number of iterations.
         * @param aKey Custom dictionary key (optional).
         */
        void SetMLEIterations(int aIterationsNumber, const std::string &aKey = "");

        /**
         * @brief Set the maximum theta vector.
         * @param aMaximumTheta Vector of doubles representing the theta values.
         * @param aKey Custom dictionary key (optional).
         */
        void SetMaximumTheta(const std::vector<double> &aMaximumTheta, const std::string &aKey = "");

        /**
         * @brief Set the total modeling execution time.
         * @param aTime Double representing the execution time.
         * @param aKey Custom dictionary key (optional).
         */
        void SetTotalModelingExecutionTime(double aTime, const std::string &aKey = "");

        /**
         * @brief Get the total modeling execution time.
         * @return Double representing the total modeling execution time.
         */
        [[nodiscard]] double GetTotalModelingExecutionTime() const;

        /**
         * @brief Get the MLOE.
         * @return The MLOE.
         *
         */
        [[nodiscard]] double GetMLOE() const;

        /**
        * @brief Get the MSPEError.
        * @return The MSPEError.
        *
        */
        [[nodiscard]] double GetMSPEError() const;

        /**
        * @brief Get the IDW error.
        * @return The the IDW error vector.
        *
        */
        [[nodiscard]] std::vector<double> GetIDWError() const;

        /**
        * @brief Get the MMOM.
        * @return The MMOM.
        *
        */
        [[nodiscard]] double GetMMOM() const;

        /**
         * @brief Get the Fisher matrix elements.
         * @return the Fisher matrix.
         *
         */
        [[nodiscard]] std::vector<double> GetFisherMatrix() const;

        /**
         * @brief Get the Predicted Missed Z matrix elements.
         * @return the Z Predicted matrix.
         *
         */
        [[nodiscard]] std::vector<double> GetPredictedMissedValues() const;

        /**
         * @brief Set the total modeling FLOPS.
         * @param aTime Double representing the FLOPS.
         * @param aKey Custom dictionary key (optional).
         */
        void SetTotalModelingFlops(double aTime, const std::string &aKey = "");

        /**
         * @brief Get the total modeling FLOPS.
         * @return Double representing the total modeling FLOPS.
         */
        [[nodiscard]] double GetTotalModelingFlops() const;

        /**
         * @brief Get the average modeling execution time.
         * @return The average execution time for data modeling.
         *
         */
        [[nodiscard]] double GetAverageModelingExecutionTime() const;

        /**
         * @brief Get the average modeling FLOPs.
         * @return The average number of FLOPs for data modeling.
         *
         */
        [[nodiscard]] double GetAverageModelingFlops() const;

        /**
         * @brief Set the value of ZMiss.
         * @param[in] aZMiss The value of ZMiss.
         *
         */
        void SetZMiss(int aZMiss, const std::string &aKey = "");

        /**
         * @brief Set the value of MSPEError.
         * @param[in] aMSPEError The value of MSPEError.
         *
         */
        void SetMSPEError(double aMSPEError, const std::string &aKey = "");

        /**
         * @brief Set the MSPE execution time.
         * @param[in] aTime The execution time.
         *
         */
        void SetMSPEExecutionTime(double aTime, const std::string &aKey = "");

        /**
         * @brief Set the MSPE number of floating-point operations (FLOPs).
         * @param[in] aFlops The number of FLOPs.
         *
         */
        void SetMSPEFlops(double aFlops, const std::string &aKey = "");

        /**
         * @brief Set the vector of IDW errors.
         * @param[in] aIDWError The vector of IDW errors.
         *
         */
        void SetIDWError(const std::vector<double> &aIDWError, const std::string &aKey = "");

        /**
         * @brief Set the value of MLOE.
         * @param[in] aMLOE The value of MLOE.
         *
         */
        void SetMLOE(double aMLOE, const std::string &aKey = "");

        /**
         * @brief Set the value of MMOM.
         * @param[in] aMMOM The value of MMOM.
         *
         */
        void SetMMOM(double aMMOM, const std::string &aKey = "");

        /**
         * @brief Set the MLOE-MMOM execution time.
         * @param[in] aTime The execution time.
         *
         */
        void SetExecutionTimeMLOEMMOM(double aTime, const std::string &aKey = "");

        /**
         * @brief Set the MLOE-MMOM matrix generation time.
         * @param[in] aTime The execution time.
         *
         */
        void SetMatrixGenerationTimeMLOEMMOM(double aTime, const std::string &aKey = "");

        /**
         * @brief Set the MLOE-MMOM cholesky factorization time.
         * @param[in] aTime The execution time.
         *
         */
        void SetFactoTimeMLOEMMOM(double aTime, const std::string &aKey = "");

        /**
        * @brief Set the MLOE-MMOM loop time.
        * @param[in] aTime The execution time.
        *
        */
        void SetLoopTimeMLOEMMOM(double aTime, const std::string &aKey = "");

        /**
         * @brief Set the MLOE-MMOM number of floating-point operations (FLOPs).
         * @param[in] aFlops The number of FLOPs.
         *
         */
        void SetFlopsMLOEMMOM(double aFlops, const std::string &aKey = "");

        /**
         * @brief Set The total execution time of the fisher tile computation.
         * @param[in] aTime The total execution time for fisher tile computation.
         *
         */
        void SetTotalFisherTime(double aTime, const std::string &aKey = "");

        /**
         * @brief Set the elements of the fisher matrix.
         * @param aFisherMatrix Elements of the fisher matrix.
         *
         */
        void SetFisherMatrix(std::vector<double> aFisherMatrix, const std::string &aKey = "");

        /**
         * @brief Set the elements of the Z missed matrix.
         * @param aPredictedValues Elements of the Predicted Z missed matrix.
         *
         */
        void SetPredictedMissedValues(std::vector<double> aPredictedValues, const std::string &aKey = "");

        /**
         * @brief Print the end summary of the results.
         *
         */
        void PrintEndSummary();

    private:
        /**
         * @brief Add result to dictionary
         * @param[in] aKey string used as title
         * @param[in] aValue string value of result
         *
         */
        void UpdateDictionary(const std::string &aKey, const std::string &aValue);

    private:
        /**
         * @brief Pointer to the singleton instance of the SyntheticGenerator class.
         */
        static Results *mpInstance;
        /// Used is synthetic.
        bool mIsSynthetic = true;
        /// Used number of generated locations.
        int mGeneratedLocationsNumber = 0;
        /// Used is logger
        bool mIsLogger = false;
        /// Used logger path
        std::string mLoggerPath;
        /// Used Data Generation Execution time.
        double mExecutionTimeDataGeneration = 0;
        /// Used Data Generation flops.
        double mFlopsDataGeneration = 0;
        /// Used MLE number of iterations.
        int mMLEIterations = 0;
        /// Used MAX theta.
        std::vector<double> mMaximumTheta;
        /// Used log likelihood value.
        double mLogLikValue = 0;
        /// Used number of Z missed values.
        int mZMiss = 0;
        /// Used MSPE error.
        double mMSPEError = 0;
        /// Used Execution time.
        double mExecutionTimeMSPE;
        /// Used flops.
        double mFlopsMSPE;
        /// Used IDW error.
        std::vector<double> mIDWError;
        /// Used MLOE.
        double mMLOE = 0;
        /// USED MMOM.
        double mMMOM = 0;
        /// Used MLOE-MMOM Execution time.
        double mExecutionTimeMLOEMMOM = 0;
        /// Used MLOE-MMOM Matrix Generation time.
        double mGenerationTimeMLOEMMOM = 0;
        /// Used MLOE-MMOM cholesky factorization time.
        double mFactoTimeMLOEMMOM = 0;
        /// Used MLOE-MMOM loop time.
        double mLoopTimeMLOEMMOM = 0;
        /// Used MLOE-MMOM flops.
        double mFlopsMLOEMMOM = 0;
        /// Used Data Modeling Execution time.
        double mTotalModelingExecutionTime = 0;
        /// Used Data Modeling Number of Flops.
        double mTotalModelingFlops = 0;
        /// Used Total Fisher Time.
        double mTotalFisherTime = 0;
        /// Fisher matrix
        std::vector<double> mFisherMatrix;
        /// Z miss values
        std::vector<double> mPredictedMissedValues;
        /// Map that holds results
        std::map<std::string, std::string> mSummaryDictionary;
    };

}//namespace exageostat
#endif //EXAGEOSTATCPP_RESULTS_HPP