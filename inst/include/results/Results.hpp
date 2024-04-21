
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
         * @brief Set the flag indicating whether the results are synthetic or not.
         * @param[in] aIsSynthetic True if the results are synthetic, false otherwise.
         *
         */
        void SetIsSynthetic(bool aIsSynthetic);

        /**
         * @brief Set the number of generated locations.
         * @param[in] aNumLocations The number of generated locations.
         *
         */
        void SetGeneratedLocationsNumber(int aNumLocations);

        /**
         * @brief Set the flag indicating whether the logger is active or not.
         * @param[in] aIsLogger True if the logger is active, false otherwise.
         *
         */
        void SetIsLogger(bool aIsLogger);

        /**
         * @brief Set the path for the logger.
         * @param[in] aLoggerPath The path for the logger.
         *
         */
        void SetLoggerPath(const std::string &aLoggerPath);

        /**
         * @brief Set the Total Data Generation execution time.
         * @param[in] aTime The execution time.
         *
         */
        void SetTotalDataGenerationExecutionTime(double aTime);

        /**
         * @brief Set the Data Generation floating-point operations (FLOPs).
         * @param[in] aFlops The number of FLOPs.
         *
         */
        void SetTotalDataGenerationFlops(double aFlops);

        /**
         * @brief Set the log-likelihood value.
         * @param[in] aLogLikValue The log-likelihood value.
         *
         */
        void SetLogLikValue(double aLogLikValue);

        /**
         * @brief Set the number of maximum likelihood estimation (MLE) iterations.
         * @param[in] aIterationsNumber The number of MLE iterations.
         *
         */
        void SetMLEIterations(int aIterationsNumber);

        /**
         * @brief Set the vector of maximum theta values.
         * @param[in] aMaximumTheta The vector of maximum theta values.
         *
         */
        void SetMaximumTheta(const std::vector<double> &aMaximumTheta);

        /**
         * @brief Set the total modeling execution time.
         * @param[in] aTime The total execution time for data modeling.
         *
         */
        void SetTotalModelingExecutionTime(double aTime);

        /**
         * @brief Get the total modeling execution time.
         * @return The total execution time for data modeling.
         *
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
         * @brief Set the total modeling FLOPs.
         * @param[in] aTime The total number of FLOPs for data modeling.
         *
         */
        void SetTotalModelingFlops(double aTime);

        /**
         * @brief Get the total modeling FLOPs.
         * @return The total number of FLOPs for data modeling.
         *
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
        void SetZMiss(int aZMiss);

        /**
         * @brief Set the value of MSPEError.
         * @param[in] aMSPEError The value of MSPEError.
         *
         */
        void SetMSPEError(double aMSPEError);

        /**
         * @brief Set the MSPE execution time.
         * @param[in] aTime The execution time.
         *
         */
        void SetMSPEExecutionTime(double aTime);

        /**
         * @brief Set the MSPE number of floating-point operations (FLOPs).
         * @param[in] aFlops The number of FLOPs.
         *
         */
        void SetMSPEFlops(double aFlops);

        /**
         * @brief Set the vector of IDW errors.
         * @param[in] aIDWError The vector of IDW errors.
         *
         */
        void SetIDWError(const std::vector<double> &aIDWError);

        /**
         * @brief Set the value of MLOE.
         * @param[in] aMLOE The value of MLOE.
         *
         */
        void SetMLOE(double aMLOE);

        /**
         * @brief Set the value of MMOM.
         * @param[in] aMMOM The value of MMOM.
         *
         */
        void SetMMOM(double aMMOM);

        /**
         * @brief Set the MLOE-MMOM execution time.
         * @param[in] aTime The execution time.
         *
         */
        void SetExecutionTimeMLOEMMOM(double aTime);

        /**
         * @brief Set the MLOE-MMOM matrix generation time.
         * @param[in] aTime The execution time.
         *
         */
        void SetMatrixGenerationTimeMLOEMMOM(double aTime);

        /**
         * @brief Set the MLOE-MMOM cholesky factorization time.
         * @param[in] aTime The execution time.
         *
         */
        void SetFactoTimeMLOEMMOM(double aTime);

        /**
        * @brief Set the MLOE-MMOM loop time.
        * @param[in] aTime The execution time.
        *
        */
        void SetLoopTimeMLOEMMOM(double aTime);

        /**
         * @brief Set the MLOE-MMOM number of floating-point operations (FLOPs).
         * @param[in] aFlops The number of FLOPs.
         *
         */
        void SetFlopsMLOEMMOM(double aFlops);

        /**
         * @brief Set The total execution time of the fisher tile computation.
         * @param[in] aTime The total execution time for fisher tile computation.
         *
         */
        void SetTotalFisherTime(double aTime);

        /**
         * @brief Set the elements of the fisher matrix.
         * @param aFisherMatrix Elements of the fisher matrix.
         *
         */
        void SetFisherMatrix(std::vector<double> aFisherMatrix);

        /**
         * @brief Set the elements of the Z missed matrix.
         * @param aPredictedValues Elements of the Predicted Z missed matrix.
         *
         */
        void SetPredictedMissedValues(std::vector<double> aPredictedValues);

        /**
         * @brief Print the end summary of the results.
         *
         */
        void PrintEndSummary();

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
    };

}//namespace exageostat
#endif //EXAGEOSTATCPP_RESULTS_HPP