
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
                SetDimension(exageostat::dataunits::Dimension aDimension);

                /**
                 * @brief Dimension getter.
                 * @return mDimension
                 */
                exageostat::dataunits::Dimension
                GetDimension();

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

                exageostat::dataunits::Dimension CheckDimensionValue(std::string aDimension);

            private:
                /// Used Dimension.
                exageostat::dataunits::Dimension mDimension = dataunits::Dimension2D;
                /// Used PGrid.
                int mPGrid;
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
