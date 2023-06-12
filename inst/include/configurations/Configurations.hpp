
/*
 * Copyright (c) 2017-2023 King Abdullah University of Science and Technology,
 * Copyright (C) 2023 by Brightskies inc,
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

#include <common/Definitions.hpp>
#include <vector>

namespace exageostat {
    namespace configurations {
        /**
         * @class Configurations
         * @brief Contains methods to set and get.
         */
        class Configurations {
        public:

            /**
             * @brief Virtual destructor to allow calls to the correct concrete destructor.
             */
            virtual ~Configurations() = default;

            /**
             * @brief
             * Set default values for input arguments
             *
             * @param[in] argc
             * The number of arguments being passed into your program from the command line.
             *
             * @param[in] argv
             * The array of arguments.
             *
             */
            virtual void
            InitializeArguments(int argc, char **argv) = 0;

            /**
             * @brief Print the usage and accepted Arguments.
             */
            virtual void
            PrintUsage() = 0;

            /**
             * @brief Problem size setter.
             * @param aProblemSize
             */
            void
            SetProblemSize(int aProblemSize);

            /**
             * @brief Problem size getter.
             * @return mProblemSize
             */
            int
            GetProblemSize() const;

            /**
            * @brief Time slot setter.
            * @param aTimeSlot
            */
            void
            SetTimeSlot(int aTimeSlot);

            /**
             * @brief Time slot getter.
             * @return mTimeSlot
             */
            int
            GetTimeSlot() const;

            /**
            * @brief
            * Computation setter.
            *
            * @param aComputation
            *
            */
            void
            SetComputation(common::Computation aComputation);

            /**
             * @brief
             * Computation getter.
             *
             * @return mComputation
             *
             */
            common::Computation
            GetComputation();

            /**
             * @brief
             * Precision setter.
             *
             * @param aComputation
             *
             */
            void
            SetPrecision(common::Precision aPrecision);

            /**
             * @brief
             * Operator getter.
             *
             * @return mOperator
             *
             */
            common::Operators
            GetOperator();

            /**
             * @brief
             * Operator setter.
             *
             * @param aOperator
             *
             */
            void
            SetOperator(common::Operators aOperator);

            /**
             * @brief
             * Precision getter.
             *
             * @return mPrecision
             *
             */
            common::Precision
            GetPrecision();

            /**
             * @brief PGrid setter.
             * @param aPGrid
             */
            void
            SetPGrid(int aPGrid);

            /**
             * @brief PGrid getter.
             * @return mPGrid
             */
            int
            GetPGrid() const;

            /**
             * @brief QGrid setter.
             * @param aQGrid
             */
            void
            SetQGrid(int aQGrid);

            /**
             * @brief QGrid getter.
             * @return mQGrid
             */
            int
            GetQGrid() const;

            /**
             * @brief Cores number setter.
             * @param aCoresNumbers
             */
            void
            SetCoresNumber(int aCoresNumbers);

            /**
             * @brief Cores numbers getter.
             * @return mpCoresNumber
             */
            int
            GetCoresNumber() const;

            /**
             * @brief GPU number setter.
             * @param aGPUs
             */
            void
            SetGPUsNumber(int aGPUsNumber);

            /**
             * @brief GPU numbers getter.
             * @return mpGPUsNumber
             */
            int
            GetGPUsNumber() const;

            /**
             * @brief P setter.
             * @param aP
             */
            void
            SetP(int aP);

            /**
             * @brief P getter.
             * @return mP
             */
            int
            GetP() const;

            /**
             * @brief Dense Tile size setter.
             * @param aTileSize
             */
            void
            SetDenseTileSize(int aTileSize);

            /**
             * @brief Dense Tile size getter.
             * @return mTileSize
             */
            int
            GetDenseTileSize() const;

            /**
             * @brief Low Tile size setter.
             * @param aTileSize
             */
            void
            SetLowTileSize(int aTileSize);

            /**
             * @brief Low tile size getter.
             * @return mTileSize
             */
            int
            GetLowTileSize() const;

            /**
             * @brief Out of Core technology setter.
             * @param aIsOOC
             */
            void
            SetIsOOC(bool aIsOOC);

            /**
             * @brief Out of Core technology getter.
             * @return mIsOOC
             */
            bool
            GetIsOOC() const;

            /**
             * @brief Max rank setter.
             * @param aMaxRank
             */
            void
            SetMaxRank(int aMaxRank);

            /**
             * @brief Max Rank getter.
             * @return mMaxRank
             */
            int
            GetMaxRank();

            /**
             * @brief Number of unknown observation to be predicted setter.
             * @param aUnknownObservationsNumber
             */
            void
            SetUnknownObservationsNb(int aUnknownObservationsNumber);

            /**
             * @brief Number of unknown observation to be predicted getter.
             * @return mUnknownObservationsNumber
             */
            int
            GetUnknownObservationsNb() const;

            /**
             * @brief Mean square error value setter.
             * @param aMeanSquareError
             */
            void
            SetMeanSquareError(double aMeanSquareError);

            /**
             * @brief Mean square error value getter.
             * @return mMeanSquareError
             */
            double
            GetMeanSquareError();

            /**
             * @brief Check indicator for approximation mode setter.
             * @param aApproximationMode
             */
            void
            SetApproximationMode(int aApproximationMode);

            /**
             * @brief Mean square error value getter.
             * @return mApproximationMode
             */
            int
            GetApproximationMode() const;

            /**
             * @brief Number of known observation values setter.
             * @param aKnownObservationsValues
             */
            void
            SetKnownObservationsValues(int aKnownObservationsValues);

            /**
             * @brief Number of known observation values getter.
             * @return mKnownObservationsValues
             */
            int
            GetKnownObservationsValues();

            /**
             * @brief Determinant values setter.
             * @param aDeterminantValue
             */
            void
            SetDeterminantValue(double aDeterminantValue);

            /**
             * @brief Determinant values getter.
             * @return mDeterminantValue
             */
            double
            GetDeterminantValue();

            /**
             * @brief Actual Observations File Path setter.
             * @param aKnownObservationsValues
             */
            void
            SetActualObservationsFilePath(std::string aKnownObservationsValues);

            /**
             * @brief Actual Observations File Path getter.
             * @return mActualObservationsFilePath
             */
            std::string
            GetActualObservationsFilePath();

            /**
             * @brief vector of C descriptors getter.
             * @return mpDescriptorC
             */
            std::vector<void *>
            &GetDescriptorC();

            /**
             * @brief vector of Z descriptors getter.
             * @return mpDescriptorZ
             */
            std::vector<void *>
            &GetDescriptorZ();

            /**
             * @brief Z copy descriptors getter.
             * @return mpDescriptorZ
             */
            void *
            &GetDescriptorZcpy();

            /**
             * @brief vector of Product descriptors getter.
             * @return mpDescriptorProduct
             */
            std::vector<void *>
            &GetDescriptorProduct();

            /**
             * @brief Determinant descriptors getter.
             * @return mpDescriptorDeterminant
             */
            void *
            &GetDescriptorDeterminant();

            /**
            * @brief vector of CD descriptors getter.
            * @return mpDescriptorCD
            */
            std::vector<void *>
            &GetDescriptorCD();

            /**
             * @brief CUV descriptors getter.
             * @return mpDescriptorCUV
             */
            std::vector<void *>
            &GetDescriptorCUV();

            /**
             * @brief Crk descriptors getter.
             * @return mpDescriptorCrk
             */
            std::vector<void *>
            &GetDescriptorCrk();

            /**
             * @brief Unknown Observations Z descriptors getter.
             * @return mpDescriptorZObservations
             */
            void *
            &GetDescriptorZObservations();

            /**
             * @brief Z Actual observations descriptors getter.
             * @return mpDescriptorZActual
             */
            void *
            &GetDescriptorZActual();

            /**
             * @brief Mean Square Error descriptors getter.
             * @return mpDescriptorMSE
             */
            void *
            &GetDescriptorMSE();

            /**
             * @brief Sequence setter.
             * @param apSequence
             */
            void
            SetSequence(void *apSequence);

            /**
             * @brief Sequence getter.
             * @return mpSequence
             */
            void *
            GetSequence();

            /**
             * @brief Request setter.
             * @param apRequest
             */
            void
            SetRequest(void *apRequest);

            /**
             * @brief Request getter.
             * @return mpRequest
             */
            void *
            GetRequest();

            /**
             * @brief Getter for the seed.
             * @return mSeed
             */
            int GetSeed();

            /**
             * @brief Getter for the seed.
             * @param aSeed The seed to set.
             */
            void SetSeed(int aSeed);

            /**
             * @brief
             * Check Numerical value.
             *
             * @param aValue
             * The input from the user side
             *
             * @return aValue
             * The int casted value.
             *
             */
            static int CheckNumericalValue(const std::string &aValue);

            /**
             * @brief
             * Check input Computation value.
             *
             * @param aValue
             * The input from the user side
             *
             * @return aComputation
             * Enum with the selected computation, Error if not exist.
             *
             */
            static common::Computation CheckComputationValue(const std::string &aValue);

            /**
             * @brief
             * Check input Precision value.
             *
             * @param aValue
             * The input from the user side
             *
             * @return aComputation
             * Enum with the selected Precision, Error if not exist.
             *
             */
            static common::Precision CheckPrecisionValue(const std::string &aValue);

            /**
             * @brief RunMode getter.
             * @return mRunMode
             */
            static exageostat::common::RunMode GetRunMode();

           /**
            * @brief RunMode setter.
            * @param aRunMode
            */
            static void SetRunMode(exageostat::common::RunMode aRunMode);

            /**
             * @brief Logger getter.
             * @return mLogger
             */
            bool GetLogger();

            /**
             * @brief Logger setter.
             * @param aLogger
             */
            void SetLogger(bool aLogger);

            /**
             * @brief Logger path getter.
             * @return mLoggerPath
             */
            std::string* GetLoggerPath();

            /**
             * @brief Logger path setter.
             * @param aLoggerPath
             */
            void SetLoggerPath(const std::string &aLoggerPath);


        protected:
            /// Used Problem size.
            int mProblemSize = 0;
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
            /// Used Operator.
            common::Operators mOperator = common::MLE;
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
            void *mpSequence;
            //// Used request
            void *mpRequest;
            /// The Seed variable, with default value = 0.
            int mSeed = 0;
            //// Used run mode
            static exageostat::common::RunMode mRunMode;
            //// Used logger
            bool mLogger = false;
            //// Used logger
            std::string mLoggerPath = MY_LOGGER_PATH;
        };

    }//namespace configurations
}//namespace exageostat

#endif //EXAGEOSTAT_CPP_CONFIGURATIONS_HPP
