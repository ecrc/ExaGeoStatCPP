
/*
 * Copyright (c) 2017-2023 King Abdullah University of Science and Technology,
 * All rights reserved.
 * ExaGeoStat is a software package, provided by King Abdullah University of Science and Technology (KAUST).
 */

/**
 * @file SyntheticDataConfigurations.hpp
 * @brief Contains the definition of the SyntheticDataConfigurations class for configuring synthetic data generation in ExaGeoStat.
 * @version 1.0.0
 * @author Sameh Abdulah
 * @date 2023-02-01
**/

#ifndef EXAGEOSTAT_CPP_SYNTHETICDATACONFIGURATIONS_HPP
#define EXAGEOSTAT_CPP_SYNTHETICDATACONFIGURATIONS_HPP

#include <algorithm>
#include <configurations/data-generation/DataConfigurations.hpp>

namespace exageostat {
    namespace configurations {
        namespace data_configurations {

            /**
             * @class SyntheticDataConfigurations
             * @brief A class for configuring synthetic data generation.
             *
             */
            class SyntheticDataConfigurations : public DataConfigurations {

            public:

                /**
                 * @brief Default constructor.
                 *
                 */
                SyntheticDataConfigurations() = default;

                /**
                 * @brief Virtual destructor to allow calls to the correct concrete destructor.
                 *
                 */
                ~SyntheticDataConfigurations() override = default;

                /**
                 * @brief Copy constructor.
                 * @param[in] aSyntheticDataConfigurations Another instance of SyntheticDataConfigurations to copy from.
                 *
                 */
                SyntheticDataConfigurations(const SyntheticDataConfigurations &aSyntheticDataConfigurations) = default;

                /**
                 * @brief Arguments constructor.
                 * @param[in] aArgC The number of arguments being passed into your program from the command line.
                 * @param[in] apArgV The array of arguments.
                 *
                 */
                SyntheticDataConfigurations(int aArgC, char **apArgV);

                /**
                 * @brief Set default values for input arguments.
                 * @copydoc Configurations::InitModuleArguments()
                 *
                 */
                void InitModuleArguments(int aArgC, char **apArgV) override;

                /**
                 * @brief Setter for the dimension.
                 * @param[in] aDimension The dimension to set.
                 * @return void
                 *
                 */
                void SetDimension(exageostat::common::Dimension aDimension);

                /**
                 * @brief Getter for the dimension.
                 * @return The current dimension.
                 *
                 */
                exageostat::common::Dimension GetDimension();

                /**
                 * @brief Checks the value of the dimension parameter.
                 * @param[in] aDimension A string representing the dimension.
                 * @return The corresponding dimension value.
                 *
                 */
                static exageostat::common::Dimension CheckDimensionValue(const std::string &aDimension);

                /**
                 * @brief Checks the value of the unknown observations parameter.
                 * @param[in] aValue A string representing the number of unknown observations.
                 * @return The corresponding integer value.
                 */
                int CheckUnknownObservationsValue(const std::string &aValue);

            private:
                /// The dimension used for data generation.
                exageostat::common::Dimension mDimension = common::Dimension2D;
            };

        }//namespace data_configurations
    }//namespace configurations
}//namespace exageostat

#endif //EXAGEOSTAT_CPP_SYNTHETICDATACONFIGURATIONS_HPP