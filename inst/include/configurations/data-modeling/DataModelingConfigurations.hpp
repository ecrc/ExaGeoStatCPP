
/*
 * Copyright (c) 2017-2023 King Abdullah University of Science and Technology,
 * All rights reserved.
 * ExaGeoStat is a software package, provided by King Abdullah University of Science and Technology (KAUST).
 */

/**
 * @file DataModelingConfigurations.hpp
 * @brief Contains the definition of the DataModelingConfigurations class for configuring data modeling in ExaGeoStat.
 * @version 1.0.0
 * @author Sameh Abdulah
 * @date 2023-06-21
**/

#include <configurations/Configurations.hpp>

#ifndef EXAGEOSTATCPP_DATAMODELINGCONFIGURATIONS_HPP
#define EXAGEOSTATCPP_DATAMODELINGCONFIGURATIONS_HPP

namespace exageostat {
    namespace configurations {
        namespace data_modeling {

            /**
             * @class DataModelingConfigurations
             * @brief A class for configuring data modeling.
             *
             */
            class DataModelingConfigurations : public Configurations {
            public:

                /**
                 * @brief Default constructor.
                 *
                 */
                DataModelingConfigurations() = default;

                /**
                 * @brief Virtual destructor to allow calls to the correct concrete destructor.
                 *
                 */
                ~DataModelingConfigurations() override = default;

                /**
                 * @brief Copy constructor.
                 * @param[in] aDataModelingConfigurations Another instance of DataModelingConfigurations to copy from.
                 *
                 */
                DataModelingConfigurations(const DataModelingConfigurations &aDataModelingConfigurations);

                /**
                 * @brief Arguments constructor.
                 * @param[in] aArgC The number of arguments being passed into your program from the command line.
                 * @param[in] apArgV The array of arguments.
                 *
                 */
                DataModelingConfigurations(int aArgC, char **apArgV);

                /**
                 * @brief Set default values for input arguments.
                 * @copydoc Configurations::InitModuleArguments()
                 *
                 */
                void InitModuleArguments(int aArgC, char **apArgV) override;

                /**
                 * Recovery File Setter
                 * @param aRecoveryFile
                 * @return void
                 */
                void
                SetRecoveryFile(std::string &aRecoveryFile);

                /**
                 * Recovery File Getter
                 * @return The Recovery File
                 */
                std::string
                GetRecoveryFile();

                /**
                 * @brief Iteration value setter.
                 * @param[in] aIterationValue
                 * @return void
                 *
                 */
                void
                SetIterationsValue(int aIterationValue);

                /**
                 * @brief Iteration value getter.
                 * @return The Iterations value.
                 *
                 */
                int
                GetIterationsValue() const;

                /**
                 * @brief Determinant value getter.
                 * @return The Determinant value.
                 *
                 */
                double *GetDeterminantValue();

                /**
                 * @brief Variance value getter.
                 * @return The Variance values.
                 *
                 */
                std::vector<double *> &GetDotProduct();

                /**
                 * @brief Variance values getter.
                 * @return The Variance values.
                 *
                 */
                std::vector<double *>GetVariance();


                /**
                 * @brief Distance Metric setter.
                 * @param apDistanceMetric
                 * @return void
                 */
                void
                SetDistanceMetric (const std::string &aDistanceMetric);

                /**
                 * @brief parse user's input to distance metric.
                 * @param aDistanceMetric
                 */
                void
                ParseDistanceMetric (const std::string &aDistanceMetric);

                /**
                 * @brief Distance Metric getter.
                 * @return Distance Metric value.
                 */
                std::string
                GetDistanceMetric ();

                /**
                 * @brief Datalog getter.
                 */
                bool
                GetLog();

                /**
                 * @brief Datalog setter.
                 * @param aLog
                 */
                void
                SetLog(bool aLog);

                /**
                 * @brief parse user's input to Data Log.
                 * @param aLog
                 */
                void
                ParseDataLog(std::string aLog);

                /**
                 * @brief Getter to Log File.
                 * @return Log File.
                 */
                FILE *
                GetFileLog();

                /**
                 * @brief Getter for the average executed time per iteration.
                 * @return The average executed time per iteration as a double.
                 */
                double
                GetAvgExecutedTimePerIteration();

                /**
                 * @brief Setter for the average executed time per iteration.
                 * @param aAvgExecTimePerIter
                 */
                void
                SetAvgExecutedTimePerIteration(double aAvgExecTimePerIter);

                /**
                 * @brief Getter for the average number of floating point operations (FLOPS) per iteration.
                 * @return The average FLOPS per iteration as a double.
                 */
                double
                GetAvgFlopsPerIteration();

                /**
                 * @brief Setter for the average number of floating point operations (FLOPS) per iteration.
                 * @param aAvgFlopsPerIter
                 */
                void
                SetAvgFlopsPerIteration(double aAvgFlopsPerIter);

                /**
                 * @brief Getter for the the final log likelihood.
                 * @return The final log likelihood as a double..
                 */
                double
                GetFinalLogLik();

                /**
                 * @brief Setter for the the final log likelihood.
                 * @param aAvgExecTimePerIter
                 */
                void
                SetFinalLogLik(double aAvgExecTimePerIter);

            private:
                /// Used Recovery File
                std::string mRecoveryFile;
                /// Used Determinant
                double *mpMatrixDeterminant = nullptr;
                /// Used Dot Product Values
                /// where the first, second, and third elements are dotp, dotp1, and dotp2.
                std::vector<double *>mDotProduct;
                /// Used Number of Iterations
                int mIterationsValue = 0;
                /// Used Distance Metric
                std::string mDistanceMetric;
                /// Used Variance values
                /// where the first, second, and third elements are variance, variance1, and variance2.
                std::vector<double *>mVariance;
                /// Used Datalog
                bool mLog = false;
                /// Log File Path.
                FILE *mpFileLog = nullptr;
                /// Average execution time per iteration (only used in verbose mode).
                double mAvgExecutedTimePerIteration = 0;
                /// Average flops per iteration (only used in verbose mode).
                double mAvgFlopsPerIteration = 0;
                /// Final likelihood
                ///TODO: should we create a result class?
                double mFinalLogLik = 0;

            };
        }//namespace data_modeling
    }//namespace configurations
}//namespace exageostat

#endif //EXAGEOSTATCPP_DATAMODELINGCONFIGURATIONS_HPP
