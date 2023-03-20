
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

#include <data-units/Helpers.hpp>

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
            SetComputation(dataunits::Computation aComputation);

            /**
             * @brief
             * Computation getter.
             *
             * @return mComputation
             *
             */
            dataunits::Computation
            GetComputation();

            /**
             * @brief
             * Precision setter.
             *
             * @param aComputation
             *
             */
            void
            SetPrecision(dataunits::Precision aPrecision);

            /**
             * @brief
             * Precision getter.
             *
             * @return mPrecision
             *
             */
            dataunits::Precision
            GetPrecision();

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
            dataunits::Computation CheckComputationValue(std::string aValue);

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
            dataunits::Precision CheckPrecisionValue(std::string aValue);

        protected:
            /// Used Problem size.
            int mProblemSize;
            /// Used Time slot.
            int mTimeSlot = 1;
            /// Used Computation.
            dataunits::Computation mComputation = dataunits::EXACT_DENSE;
            /// Used Precision.
            dataunits::Precision mPrecision = dataunits::SINGLE;

        };

    }//namespace configurations
}//namespace exageostat

#endif //EXAGEOSTAT_CPP_CONFIGURATIONS_HPP
