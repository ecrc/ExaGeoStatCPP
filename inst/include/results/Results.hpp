
// Copyright (c) 2017-2023 King Abdullah University of Science and Technology,
// All rights reserved.
// ExaGeoStat is a software package, provided by King Abdullah University of Science and Technology (KAUST).

/**
 * @file Results.hpp
 * @brief Defines the Results class for storing and accessing result data.
 * @version 1.0.0
 * @author Mahmoud ElKarargy
 * @date 2023-09-14
**/

#ifndef EXAGEOSTATCPP_RESULTS_HPP
#define EXAGEOSTATCPP_RESULTS_HPP

#include <iostream>
#include <vector>

namespace exageostat {
    namespace results {

        class Results {

        public:
            /**
             * @brief Get a pointer to the singleton instance of the Results class.
             * @return A pointer to the instance of the Results class.
             */
            static Results *GetInstance();

            /**
             * @brief Set the flag indicating whether the results are synthetic or not.
             * @param[in] aIsSynthetic True if the results are synthetic, false otherwise.
             */
            void SetIsSynthetic(bool aIsSynthetic);

            /**
             * @brief Set the number of generated locations.
             * @param[in] aNumLocations The number of generated locations.
             */
            void SetGeneratedLocationsNumber(int aNumLocations);

            /**
             * @brief Set the flag indicating whether the logger is active or not.
             * @param[in] aIsLogger True if the logger is active, false otherwise.
             */
            void SetIsLogger(bool aIsLogger);

            /**
             * @brief Set the path for the logger.
             * @param[in] aLoggerPath The path for the logger.
             */
            void SetLoggerPath(const std::string &aLoggerPath);

            /**
             * @brief Set the log-likelihood value.
             * @param[in] aLogLikValue The log-likelihood value.
             */
            void SetLogLikValue(double aLogLikValue);

            /**
             * @brief Set the number of maximum likelihood estimation (MLE) iterations.
             * @param[in] aIterationsNumber The number of MLE iterations.
             */
            void SetMLEIterations(int aIterationsNumber);

            /**
             * @brief Set the vector of maximum theta values.
             * @param[in] aMaximumTheta The vector of maximum theta values.
             */
            void SetMaximumTheta(const std::vector<double> &aMaximumTheta);

            /**
             * @brief Set the value of ZMiss.
             * @param[in] aZMiss The value of ZMiss.
             */
            void SetZMiss(int aZMiss);

            /**
             * @brief Set the value of MSPEError.
             * @param[in] aMSPEError The value of MSPEError.
             */
            void SetMSPEError(double aMSPEError);

            /**
             * @brief Set the execution time.
             * @param[in] aExecutionTime The execution time.
             */
            void SetExecutionTime(double aExecutionTime);

            /**
             * @brief Set the number of floating-point operations (FLOPs).
             * @param[in] aFlops The number of FLOPs.
             */
            void SetFlops(double aFlops);

            /**
             * @brief Set the vector of IDW errors.
             * @param[in] aIDWError The vector of IDW errors.
             */
            void SetIDWError(const std::vector<double> &aIDWError);

            /**
             * @brief Set the value of MLOE.
             * @param[in] aMLOE The value of MLOE.
             */
            void SetMLOE(double aMLOE);

            /**
             * @brief Set the value of MMOM.
             * @param[in] aMMOM The value of MMOM.
             */
            void SetMMOM(double aMMOM);

            /**
             * @brief Print the end summary of the results.
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
            double mExecutionTime;
            /// Used flops.
            double mFlops;
            /// Used IDW error.
            std::vector<double> mIDWError;
            /// Used MLOE.
            double mMLOE = 0;
            /// USED MMOM.
            double mMMOM = 0;
        };

    }//namespace results
}//namespace exageostat
#endif //EXAGEOSTATCPP_RESULTS_HPP
