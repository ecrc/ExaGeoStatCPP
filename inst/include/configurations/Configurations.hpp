
/*
 * Copyright (c) 2017-2023 King Abdullah University of Science and Technology,
 * Copyright (C) 2023 by Brightskies inc,
 * All rights reserved.
 * ExaGeoStat is a software package, provided by King Abdullah University of Science and Technology (KAUST).
 */

/**
* @file Configurations.hpp
* @version 1.0.0
* @author Sameh Abdulah
* @date 2023-01-31
**/

#ifndef EXAGEOSTAT_CPP_CONFIGURATIONS_HPP
#define EXAGEOSTAT_CPP_CONFIGURATIONS_HPP

#include <common/Definitions.hpp>
#include <vector>

namespace exageostat {
    namespace configurations {

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
            GetProblemSize();

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
            GetTimeSlot();

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
            GetPGrid();

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
            GetQGrid();

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
            GetP();

            /**
             * @brief Tile size setter.
             * @param aTileSize
             */
            void
            SetTileSize(int aTileSize);

            /**
             * @brief tile size getter.
             * @return mTileSize
             */
            int
            GetTileSize();

            /**
             * @brief vector of C descriptors getter.
             * @return mpDescriptorC
             */
            std::vector<void *>
            GetDescriptorC();

            /**
             * @brief vector of Z descriptors getter.
             * @return mpDescriptorZ
             */
            std::vector<void *>
            GetDescriptorZ();

            /**
             * @brief Z copy descriptors getter.
             * @return mpDescriptorZ
             */
            void *
            GetDescriptorZcpy();

            /**
             * @brief vector of Product descriptors getter.
             * @return mpDescriptorProduct
             */
            std::vector<void *>
            GetDescriptorProduct();

            /**
             * @brief Determinant descriptors getter.
             * @return mpDescriptorDeterminant
             */
            void *
            GetDescriptorDeterminant();

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
            int CheckNumericalValue(std::string aValue);

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
            common::Computation CheckComputationValue(std::string aValue);

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
            common::Precision CheckPrecisionValue(std::string aValue);

        protected:
            /// Used Problem size.
            int mProblemSize;
            /// Used Time slot.
            int mTimeSlot = 1;
            /// Used PGrid.
            int mPGrid = 1;
            /// Used QGrid.
            int mQGrid = 1;
            /// Used P.
            int mP = 1;
            //// Used Tile Size.
            int mTileSize;
            //// Used vectors of C descriptor.
            std::vector<void *> mpDescriptorC;
            //// Used vectors of Z descriptor.
            std::vector<void *> mpDescriptorZ;
            //// Used copy Z descriptor.
            void * mpDescriptorZcpy;
            //// Used vectors of product descriptor.
            std::vector<void *> mpDescriptorProduct;
            //// Used Determinant descriptor.
            void * mpDescriptorDeterminant;
            /// Used Computation.
            common::Computation mComputation = common::EXACT_DENSE;
            /// Used Precision.
            common::Precision mPrecision = common::SINGLE;

        };

    }//namespace configurations
}//namespace exageostat

#endif //EXAGEOSTAT_CPP_CONFIGURATIONS_HPP
