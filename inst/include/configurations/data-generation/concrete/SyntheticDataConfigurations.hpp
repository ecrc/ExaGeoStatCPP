
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
            /**
             * @class SyntheticDataConfigurations
             * @brief A class for configuring synthetic data generation.
             */
            class SyntheticDataConfigurations : public DataConfigurations {

            public:
                /**
                 * @brief Default constructor.
                 */
                SyntheticDataConfigurations() = default;

                /**
                * @brief Constructor that takes command line arguments.
                * @param argc The number of arguments being passed into the program from the command line.
                * @param argv The array of arguments.
                */
                SyntheticDataConfigurations(int argc, char **argv);

                /**
                * @brief Constructor that takes a JSON file path.
                * @param JSON_path The path to the JSON file.
                */
                explicit SyntheticDataConfigurations(std::string JSON_path);

                /**
                 * @brief Virtual destructor to allow calls to the correct concrete destructor.
                 */
                ~SyntheticDataConfigurations() override = default;

                /**
                 * @brief Set default values for input arguments.
                 * @param[in] argc The number of arguments being passed into the program from the command line.
                 * @param[in] argv The array of arguments.
                 */
                void InitializeArguments(int argc, char **argv) override;

                /**
                 * @brief Print the usage and accepted arguments.
                 */
                void PrintUsage() override;

                /**
                 * @brief Setter for the dimension.
                 * @param aDimension The dimension to set.
                 */
                void SetDimension(exageostat::common::Dimension aDimension);

                /**
                 * @brief Getter for the dimension.
                 * @return The current dimension.
                 */
                exageostat::common::Dimension GetDimension();

                /**
                 * @brief Checks the value of the dimension parameter.
                 * @param aDimension A string representing the dimension.
                 * @return The corresponding dimension value.
                 */
                static exageostat::common::Dimension CheckDimensionValue(const std::string& aDimension);

                /**
                 * @brief Checks the value of the unknown observations parameter.
                 * @param aValue A string representing the number of unknown observations.
                 * @return The corresponding integer value.
                 */
                int CheckUnknownObservationsValue(std::string aValue);

            private:
                /// The dimension used for data generation.
                exageostat::common::Dimension mDimension = common::Dimension2D;
                /// The X coordinates of the locations.
                double *mpLocationX{};
                /// The Y coordinates of the locations.
                double *mpLocationY{};
                /// The Z coordinates of the locations.
                double *mpLocationZ{};

            };

        }//namespace data_configurations
    }//namespace configurations
}//namespace exageostat

#endif //EXAGEOSTAT_CPP_SYNTHETICDATACONFIGURATIONS_HPP