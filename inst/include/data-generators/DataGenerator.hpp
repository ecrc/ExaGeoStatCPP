
/**
 * @file DataGenerator.hpp
 * @brief Contains the declaration of the SyntheticDataConfigurations class.
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

        /**
         * @class DataGenerator
         * @brief Contains methods to set and get.
         */
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
            GenerateLocations() = 0;

            virtual void
            GenerateDescriptors() = 0;

            virtual void
            GenerateObservations() = 0;

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
             * Gets data locations class.
             *
             * @return mpLocations
             * Pointer to locations object.
             */
            dataunits::Locations *
            GetLocations();



            virtual ~DataGenerator() = default;


        protected:
            /// Used Synthetic Configuration.
            configurations::data_configurations::SyntheticDataConfigurations *mpConfigurations{}; // Pointer to SyntheticDataConfigurations object
            /// Used Locations
            dataunits::Locations * mpLocations{}; // Pointer to Locations object
            /// Used Kernel
            exageostat::kernels::Kernel * mpKernel;
        };
    }//namespace generators
}//namespace exageostat

#endif //EXAGEOSTAT_CPP_DATAGENERATOR_HPP