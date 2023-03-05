
/**
 * @file SyntheticGenerator.hpp
 *
 * @version 1.0.0
 * @author Sameh Abdulah
 * @date 2023-02-14
**/

#ifndef EXAGEOSTAT_CPP_SYNTHETICGENERATOR_HPP
#define EXAGEOSTAT_CPP_SYNTHETICGENERATOR_HPP

#include <data-generators/DataGenerator.hpp>

namespace exageostat {
    namespace generators {
        namespace Synthetic {
            class SyntheticGenerator : public DataGenerator{

            public:
//                /**
//                 * @brief
//                 * Constructor for Synthetic Generation.
//                 *
//                 * @param[in] aConfigurations
//                 * The Synthetic data configuration inputs.
//                 */
//                SyntheticGenerator(configurations::data_configurations::SyntheticDataConfigurations aConfigurations);
                /**
                 * @brief Default constructor.
                 */
                SyntheticGenerator() = default;

                /**
                 * @brief
                 * Virtual destructor to allow calls to the correct concrete destructor.
                 */
                virtual ~SyntheticGenerator() = default;

                /**
                 * @brief
                 * Set default values for input arguments
                 *
                 * @param[in] aN
                 * The problem size divided by P-Grid.
                 *
                 * @param[in] aSeed
                 * The input seed.
                 *
                 */
                void GenerateLocations(int aN, int aSeed);

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
                dataunits::Locations
                InitializeLocations(dataunits::Locations aLocations);

                void Print();

            };

        }//namespace Synthetic
    }//namespace generators
}//namespace exageostat

#endif //EXAGEOSTAT_CPP_SYNTHETICGENERATOR_HPP
