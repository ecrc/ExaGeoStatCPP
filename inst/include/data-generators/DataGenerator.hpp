
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
#include <configurations/data-generation/concrete/SyntheticDataConfigurations.h>

namespace exageostat {
    namespace generators {

        class DataGenerator {
        public:

            /**
             * @brief Virtual destructor to allow calls to the correct concrete destructor.
             */
            virtual ~DataGenerator() = default;
            /**
             * @brief Initialize data locations.
             * @param aLocations
             * @return Locations
             */
            virtual dataunits::Locations
            InitializeLocations(dataunits::Locations aLocations) = 0;

            virtual void Print() = 0;

            DataGenerator *createGenerator(configurations::data_configurations::SyntheticDataConfigurations *aConfigurations);



        protected:
            configurations::data_configurations::SyntheticDataConfigurations mConfigurations;

        };
    }//namespace generators
}//namespace exageostat

#endif //EXAGEOSTAT_CPP_DATAGENERATOR_HPP
