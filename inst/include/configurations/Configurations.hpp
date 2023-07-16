
/*
 * Copyright (c) 2017-2023 King Abdullah University of Science and Technology,
 * All rights reserved.
 * ExaGeoStat is a software package, provided by King Abdullah University of Science and Technology (KAUST).
 */

/**
* @file Configurations.hpp
* @version 1.0.0
* @brief Contains the declaration of the Configurations class and its member functions.
* @author Sameh Abdulah
* @date 2023-01-31
**/

#ifndef EXAGEOSTAT_CPP_CONFIGURATIONS_HPP
#define EXAGEOSTAT_CPP_CONFIGURATIONS_HPP

#include <vector>

#include <common/Definitions.hpp>

namespace exageostat {
    namespace configurations {
        /**
         * @class Configurations
         * @brief Contains methods to set and get.
         *
         */
        class Configurations {
        public:

            /**
             * @brief Virtual destructor to allow calls to the correct concrete destructor.
             *
             */
            virtual ~Configurations() = default;

            /**
             * @brief Initialize the module arguments.
             * @param[in] aArgC The number of arguments being passed into the program from the command line.
             * @param[in] apArgV The array of arguments.
             * @details This method initializes the command line arguments for the module.
             *
             */
            virtual void InitModuleArguments(int aArgC, char **apArgV) = 0;

            /**
             * @brief Set default values for input arguments.
             * @param[in] aArgC The number of arguments being passed into your program from the command line.
             * @param[in] apArgV The array of arguments.
             * @return void
             *
             */
            void InitializeArguments(int aArgC, char **apArgV);

            /**
             * @brief Print the usage and accepted Arguments.
             * @return void
             *
             */
            static void PrintUsage();

            /**
             * @brief Setter for the kernel.
             * @param[in] aKernel The kernel to set.
             * @return void
             *
             */
            void SetKernel(const std::string &aKernel);

            /**
             * @brief Getter for the kernel.
             * @return The kernel.
             *
             */
            std::string GetKernel() const;

            /**
             * @brief Problem size setter.
             * @param[in] aProblemSize
             * @return void
             *
             */
            void
            SetProblemSize(int aProblemSize);

            /**
             * @brief Setter for the number of parameters.
             * @param[in] aParameterNumbers The number of parameters to set.
             * @return void
             *
             */
            void SetParametersNumber(int aParameterNumbers);

            /**
             * @brief Getter for the number of parameters.
             * @return The number of parameters.
             *
             */
            int GetParametersNumber() const;

            /**
             * @brief Problem size getter.
             * @return The problem size.
             *
             */
            int
            GetProblemSize() const;

            /**
             * @brief Time slot setter.
             * @param[in] aTimeSlot The time slot to set.
             * @return void
             *
             */
            void
            SetTimeSlot(int aTimeSlot);

            /**
             * @brief Time slot getter.
             * @return The time slot.
             *
             */
            int
            GetTimeSlot() const;

            /**
             * @brief Computation setter.
             * @param[in] aComputation The computation to set.
             * @return void
             *
             */
            void
            SetComputation(common::Computation aComputation);

            /**
             * @brief
             * Computation getter.
             * @return The computation.
             *
             */
            common::Computation GetComputation() const;

            /**
             * @brief Precision setter.
             * @param[in] aPrecision The precision to set.
             * @return void
             *
             */
            void
            SetPrecision(common::Precision aPrecision);

            /**
             * @brief
             * Precision getter.
             * @return The precision.
             *
             */
            common::Precision GetPrecision() const;

            /**
             * @brief PGrid setter.
             * @param[in] aPGrid The PGrid to set.
             * @return void
             *
             */
            void
            SetPGrid(int aPGrid);

            /**
             * @brief PGrid getter.
             * @return The PGrid.
             *
             */
            int GetPGrid() const;

            /**
             * @brief QGrid setter.
             * @param[in] aQGrid The QGrid to set.
             * @return void
             *
             */
            void
            SetQGrid(int aQGrid);

            /**
             * @brief QGrid getter.
             * @return The QGrid.
             *
             */
            int
            GetQGrid() const;

            /**
             * @brief Cores number setter.
             * @param[in] aCoresNumbers The number of cores to set.
             * @return void
             *
             */
            void
            SetCoresNumber(int aCoresNumbers);

            /**
             * @brief Cores numbers getter.
             * @return The number of cores.
             *
             */
            int
            GetCoresNumber() const;

            /**
             * @brief GPU number setter.
             * @param[in] aGPUsNumber The number of GPUs to set.
             * @return void
             *
             */
            void
            SetGPUsNumber(int aGPUsNumber);

            /**
             * @brief GPUnumber getter.
             * @return The number of GPUs.
             *
             */
            int
            GetGPUsNumber() const;

            /**
             * @brief P setter.
             * @param[in] aP The P value to set.
             * @return void
             *
             */
            void
            SetP(int aP);

            /**
             * @brief P getter.
             * @return The P value.
             *
             */
            int
            GetP() const;

            /**
             * @brief Dense Tile size setter.
             * @param[in] aTileSize The dense Tile size to set.
             * @return void
             *
             */
            void
            SetDenseTileSize(int aTileSize);

            /**
             * @brief Dense Tile size getter.
             * @return The dense Tile size.
             *
             */
            int
            GetDenseTileSize() const;

            /**
             * @brief Low Tile size setter.
             * @param[in] aTileSize The low Tile size to set.
             * @return void
             *
             */
            void
            SetLowTileSize(int aTileSize);

            /**
             * @brief Low tile size getter.
             * @return The low Tile size.
             *
             */
            int
            GetLowTileSize() const;

            /**
             * @brief Out of Core technology setter.
             * @param[in] aIsOOC Flag indicating whether out of core technology should be used.
             * @return void
             *
             */
            void
            SetIsOOC(bool aIsOOC);

            /**
             * @brief Out of Core technology getter.
             * @return Flag indicating whether out of core technology is being used.
             *
             */
            bool
            GetIsOOC() const;

            /**
             * @brief Max rank setter.
             * @param aMaxRank The maximum rank to set.
             * @return void
             *
             */
            void
            SetMaxRank(int aMaxRank);

            /**
             * @brief Getter for the maximum rank.
             * @return The maximum rank.
             */
            int GetMaxRank() const;

            /**
             * @brief Setter for the number of unknown observations to be predicted.
             * @param[in] aUnknownObservationsNumber The number of unknown observations to set.
             * @return void
             *
             */
            void SetUnknownObservationsNb(int aUnknownObservationsNumber);

            /**
             * @brief Getter for the number of unknown observations to be predicted.
             * @return The number of unknown observations.
             */
            int GetUnknownObservationsNb() const;

            /**
             * @brief Setter for the mean square error value.
             * @param[in] aMeanSquareError The mean square error value to set.
             * @return void
             *
             */
            void SetMeanSquareError(double aMeanSquareError);

            /**
             * @brief Getter for the mean square error value.
             * @return The mean square error value.
             *
             */
            double GetMeanSquareError() const;

            /**
             * @brief Setter for the approximation mode.
             * @param[in] aApproximationMode The approximation mode to set.
             * @return void
             */
            void SetApproximationMode(int aApproximationMode);

            /**
             * @brief Getter for the approximation mode.
             * @return The approximation mode.
             *
             */
            int GetApproximationMode() const;

            /**
             * @brief Setter for the number of known observation values.
             * @param[in] aKnownObservationsValues The number of known observation values to set.
             * @return void
             *
             */
            void SetKnownObservationsValues(int aKnownObservationsValues);

            /**
             * @brief Getter for the number of known observation values.
             * @return The number of known observation values.
             *
             */
            int GetKnownObservationsValues() const;

            /**
             * @brief Setter for the determinant value.
             * @param[in] aDeterminantValue The determinant value to set.
             * @return void
             *
             */
            void SetDeterminantValue(double aDeterminantValue);

            /**
             * @brief Getter for the determinant value.
             * @return The determinant value.
             */
            double GetDeterminantValue() const;

            /**
             * @brief Setter for the actual observations file path.
             * @param[in] aActualObservationsFilePath The actual observations file path to set.
             * @return void
             *
             */
            void SetActualObservationsFilePath(const std::string &aActualObservationsFilePath);

            /**
             * @brief Getter for the actual observations file path.
             * @return The actual observations file path.
             *
             */
            std::string GetActualObservationsFilePath() const;

            /**
             * @brief Getter for the vector of C descriptors.
             * @return A reference to the vector of C descriptors.
             *
             */
            std::vector<void *> &GetDescriptorC();

            /**
             * @brief Getter for the vector of Z descriptors.
             * @return A reference to the vector of Z descriptors.
             *
             */
            std::vector<void *> &GetDescriptorZ();

            /**
             * @brief Getter for the Z copy descriptor.
             * @return A reference to the Z copy descriptor.
             *
             */
            void *&GetDescriptorZcpy();

            /**
             * @brief Getter for the vector of Product descriptors.
             * @return A reference to the vector of Product descriptors.
             *
             */
            std::vector<void *> &GetDescriptorProduct();

            /**
             * @brief Getter for the Determinant descriptor.
             * @return A reference to the Determinant descriptor.
             *
             */
            void *&GetDescriptorDeterminant();

            /**
             * @brief Getter for the vector of CD descriptors.
             * @return A reference to the vector of CD descriptors.
             *
             */
            std::vector<void *> &GetDescriptorCD();

            /**
             * @brief Getter for the vector of CUV descriptors.
             * @return A reference to the vector of CUV descriptors.
             *
             */
            std::vector<void *> &GetDescriptorCUV();

            /**
             * @brief Getter for the vector of Crk descriptors.
             * @return A reference to the vector of Crk descriptors.
             *
             */
            std::vector<void *> &GetDescriptorCrk();

            /**
             * @brief Getter for the Unknown Observations Z descriptor.
             * @return A reference to the Unknown Observations Z descriptor.
             *
             */
            void *&GetDescriptorZObservations();

            /**
             * @brief Getter for the Z Actual observations descriptor.
             * @return A reference to the Z Actual observations descriptor.
             *
             */
            void *&GetDescriptorZActual();

            /**
             * @brief Getter for the Mean Square Error descriptor.
             * @returnA reference to the Mean Square Error descriptor.
             *
             */
            void *&GetDescriptorMSE();

            /**
             * @brief Setter for the sequence.
             * @param[in] apSequence Pointer to the sequence to set.
             * @return void
             *
             */
            void SetSequence(void *apSequence);

            /**
             * @brief Getter for the sequence.
             * @return Pointer to the sequence.
             *
             */
            void *GetSequence();

            /**
             * @brief Setter for the request.
             * @param[in] apRequest Pointer to the request to set.
             * @return void
             *
             */
            void SetRequest(void *apRequest);

            /**
             * @brief Getter for the request.
             * @return Pointer to the request.
             *
             */
            void *GetRequest();

            /**
             * @brief Getter for the seed.
             * @return The seed.
             *
             */
            int GetSeed() const;

            /**
             * @brief Setter for the seed.
             * @param[in] aSeed The seed to set.
             * @return void
             *
             */
            void SetSeed(int aSeed);

            /**
             * @brief Getter for the logger.
             * @return The logger.
             *
             */
            bool GetLogger() const;

            /**
             * @brief Setter for the logger.
             * @param[in] aLogger The logger to set.
             * @return void
             *
             */
            void SetLogger(bool aLogger);

            /**
             * @brief Getter for the logger path.
             * @return Pointer to the logger path.
             *
             */
            std::string *GetLoggerPath();

            /**
             * @brief Setter for the logger path.
             * @param[in] aLoggerPath The logger path to set.
             * @return void
             *
             */
            void SetLoggerPath(const std::string &aLoggerPath);

            /**
             * @brief Check if input value is numerical.
             * @param[in] aValue The input from the user side.
             * @return The int casted value.
             *
             */
            static int CheckNumericalValue(const std::string &aValue);

            /**
             * @brief Check input computation value.
             * @param[in] aValue The input from the user side.
             * @return Enum with the selected computation, Error if not exist.
             *
             */
            static common::Computation CheckComputationValue(const std::string &aValue);

            /**
             * @brief Check input precision value.
             * @param[in] aValue The input from the user side.
             * @return Enum with the selected Precision, Error if not exist.
             *
             */
            static common::Precision CheckPrecisionValue(const std::string &aValue);

            /**
             * @brief Getter for the run mode.
             * @return The run mode.
             *
             */
            static exageostat::common::RunMode GetRunMode();

            /**
             * @brief Setter for the run mode.
             * @param[in] aRunMode The run mode to set.
             * @return void
             *
             */
            static void SetRunMode(exageostat::common::RunMode aRunMode);
            /**
             * @brief Setter for the initial theta.
             * @param[in,out] apTheta A pointer to an array of initial theta values to set.
             * @return void
             *
             */
            void SetInitialTheta(std::vector<double> &apTheta);

            /**
             * @brief Getter for the initial theta.
             * @return A pointer to the array of initial theta values.
             *
             */
            std::vector<double> &GetInitialTheta();

            /**
             * @brief Checks the run mode and sets the verbosity level.
             * @param[in] aRunMode A string representing the desired run mode ("verbose" or "standard").
             * @throws std::range_error if the input string isnot "verbose" or "standard".
             * @return void
             *
             */
            static void ParseRunMode(const std::string &aRunMode);

            /**
             * @brief Checks if the kernel value is valid.
             * @param[in] aKernel The kernel to check.
             * @return void
             *
             */
            void CheckKernelValue(const std::string &aKernel);

            /**
             * @brief Checks if a given string is in camel case format.
             * @param[in] aString The string to check.
             * @return true if the string is in camel case format, false otherwise.
             *
             */
            static bool IsCamelCase(const std::string &aString);

            /**
             * @brief Parses a string of theta values and returns an array of doubles.
             * @param[in] aInputValues The input string of theta values.
             * @return A vector of parsed theta values.
             *
             */
            static std::vector<double> ParseTheta(const std::string &aInputValues);

        protected:

            /// Used boolean to identify the first use of initialize arguments
            static bool mIsInitialized;
            /// Used Argument counter
            int mArgC = 0;
            /// Used Argument vectors
            char **mpArgV = nullptr;
            /// The kernel to use.
            std::string mKernel;
            /// Used Problem size.
            int mProblemSize = 0;
            //// The number of parameters to use.
            int mParametersNumber = 0;
            /// Used Time slot.
            int mTimeSlot = 1;
            /// Used PGrid.
            int mPGrid = 1;
            /// Used QGrid.
            int mQGrid = 1;
            /// Used Number of cores.
            int mCoresNumber = 1;
            /// Used Number of GPUs.
            int mGPUsNumber = 0;
            /// Used P.
            int mP = 1;
            //// Used Dense Tile Size.
            int mDenseTileSize = 0;
            //// Used Low Tile Size.
            int mLowTileSize = 0;
            //// Used Out of core technology.
            bool mIsOOC = false;
            //// Used max rank
            int mMaxRank = 1;
            //// Used number of unknown observation to be predicted - NZmiss
            int mUnknownObservationsNumber = 0;
            //// Used number of known observed values. -nZobs
            int mKnownObservationsValues = 0;
            //// Used Approximation mode values.
            int mApproximationMode = 1;
            //// Used Mean Square Error values.
            double mMeanSquareError = 0.0;
            //// Determinant value.
            double mDeterminantValue = 0.0;
            //// Used Actual observations file path in the case of prediction value.
            std::string mActualObservationsFilePath;
            /// Used Computation.
            common::Computation mComputation = common::EXACT_DENSE;
            /// Used Precision.
            common::Precision mPrecision = common::SINGLE;
            /// Used initial theta values.
            std::vector<double> mInitialTheta;
            //// Used vectors of C descriptor.
            std::vector<void *> mpDescriptorC;
            //// Used vectors of Z descriptor.
            std::vector<void *> mpDescriptorZ;
            //// Used copy Z descriptor.
            void *mpDescriptorZcpy = nullptr;
            //// Used Determinant descriptor.
            void *mpDescriptorDeterminant = nullptr;
            //// Used Z observations descriptor.
            void *mpDescriptorZObservations = nullptr;
            //// Used vectors of product descriptor.
            std::vector<void *> mpDescriptorProduct;
            //// Used vectors of CD descriptor.
            std::vector<void *> mpDescriptorCD = {nullptr, nullptr, nullptr};
            //// Used vectors of CUV descriptor.
            std::vector<void *> mpDescriptorCUV = {nullptr, nullptr, nullptr};
            //// Used vectors of Crk descriptor.
            std::vector<void *> mpDescriptorCrk = {nullptr, nullptr, nullptr};
            //// Used MSE descriptor.
            void *mpDescriptorMSE = nullptr;
            //// Used Z Actual observations descriptor.
            void *mpDescriptorZActual = nullptr;
            ////  Used sequence
            void *mpSequence = nullptr;
            //// Used request
            void *mpRequest = nullptr;
            /// The Seed variable, with default value = 0.
            int mSeed = 0;
            //// Used run mode
            static exageostat::common::RunMode mRunMode;
            //// Used logger
            bool mLogger = false;
            //// Used logger path.
            std::string mLoggerPath;
        };

    }//namespace configurations
}//namespace exageostat

#endif //EXAGEOSTAT_CPP_CONFIGURATIONS_HPP
