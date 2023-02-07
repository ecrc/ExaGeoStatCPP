
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

#include <iostream>
#ifndef EXAGEOSTAT_CPP_CONFIGURATIONS_HPP
#define EXAGEOSTAT_CPP_CONFIGURATIONS_HPP

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
            InitializeArguments(int argc, char** argv) = 0;
            /**
             * @brief Print the usage and accepted Arguments.
             */
            virtual void
            PrintUsage() = 0;
            /**
             * @brief Check Numerical value.
             * @param aValue
             */
            int CheckNumericalValue(std::string aValue);
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

        protected:
            /// Problem size.
            int mProblemSize;
        };

    }//namespace configurations
}//namespace exageostat

#endif //EXAGEOSTAT_CPP_CONFIGURATIONS_HPP
