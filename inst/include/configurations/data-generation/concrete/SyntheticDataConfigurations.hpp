
/*
 * Copyright (c) 2017-2023 King Abdullah University of Science and Technology,
 * Copyright (C) 2023 by Brightskies inc,
 * All rights reserved.
 * ExaGeoStat is a software package, provided by King Abdullah University of Science and Technology (KAUST).
 */

/**
 * @file SyntheticDataConfigurations.hpp
 *
 * @version 1.0.0
 * @author Sameh Abdulah
 * @date 2023-02-01
**/

#ifndef EXAGEOSTAT_CPP_SYNTHETICDATACONFIGURATIONS_HPP
#define EXAGEOSTAT_CPP_SYNTHETICDATACONFIGURATIONS_HPP

#include <iostream>
#include <configurations/data-generation/DataConfigurations.hpp>

namespace exageostat {
    namespace configurations {
        namespace data_configurations {
            class SyntheticDataConfigurations : public DataConfigurations {

            public:
                /**
                 * @brief Default constructor.
                 */
                SyntheticDataConfigurations() = default;

                /**
                * @brief Command line constructor.
                */
                SyntheticDataConfigurations(int argc, char **argv);

                /**
                * @brief JSON file constructor.
                */
                SyntheticDataConfigurations(std::string JSON_path);

                /**
                 * @brief Virtual destructor to allow calls to the correct concrete destructor.
                 */
                virtual ~SyntheticDataConfigurations() = default;

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
                void InitializeArguments(int argc, char **argv) override;

                /**
                 * @brief Print the usage and accepted Arguments.
                 */
                void
                PrintUsage() override;

                /**
                 * @brief Dimension setter.
                 * @param aDimension
                 */
                void
                SetDimension(exageostat::common::Dimension aDimension);

                /**
                 * @brief Dimension getter.
                 * @return mDimension
                 */
                exageostat::common::Dimension
                GetDimension();

                exageostat::common::Dimension CheckDimensionValue(std::string aDimension);
                int CheckUnknownObservationsValue(std::string aValue);

            private:
                /// Used Dimension.
                exageostat::common::Dimension mDimension = common::Dimension2D;
                /// Used Location X.
                double *mpLocationX;
                /// Used Location Y.
                double *mpLocationY;
                /// Used Location Z.
                double *mpLocationZ;

            };

        }//namespace data_configurations
    }//namespace configurations
}//namespace exageostat

#endif //EXAGEOSTAT_CPP_SYNTHETICDATACONFIGURATIONS_HPP
