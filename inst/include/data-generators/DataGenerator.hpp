
/**
 * @file DataGenerator.hpp
 *
 * @version 1.0.0
 * @author Sameh Abdulah
 * @date 2023-02-14
**/

#ifndef EXAGEOSTAT_CPP_DATAGENERATOR_HPP
#define EXAGEOSTAT_CPP_DATAGENERATOR_HPP

#include <data-units/Locations.hpp>
#include <configurations/data-generation/concrete/SyntheticDataConfigurations.hpp>
#include <memory>

namespace exageostat {
    namespace generators {

        class DataGenerator {
        public:

            /**
             * @brief
             * Initialize data locations.
             *
             * @param[in] aLocations
             * X, Y and Z variables.
             *
             * @return aLocations
             * The modified X, Y and Z variables.
             */
            virtual void
            InitializeLocations(int aTimeSlots) = 0;

            /**
             * @brief
             * Factory creation, Whether it's Synthetic or Real data.
             *
             * @param[in] apConfigurations
             *  Pointer to Synthetic data Configurations.
             *
             * @return DataGenerator
             * Unique Pointer to the created type of Data Generators.
             */
            static std::unique_ptr<DataGenerator>
            CreateGenerator(configurations::data_configurations::SyntheticDataConfigurations *apConfigurations);

            /**
             * @brief
             * Configuration map setter.
             *
             * @param apConfigurations
             * Argument pointer to Synthetic Data generation configuration map
             *
             */
            void
            SetConfigurations(configurations::data_configurations::SyntheticDataConfigurations *apConfigurations);

            void
            InitLocationsClass();

            dataunits::Locations *
            GetLocations();

            virtual void Print() = 0;


        protected:
            /// Used Synthetic Configuration.
            configurations::data_configurations::SyntheticDataConfigurations *mpConfigurations{};
            /// Used Locations
            dataunits::Locations * mpLocations;
        };
    }//namespace generators
}//namespace exageostat

#endif //EXAGEOSTAT_CPP_DATAGENERATOR_HPP
