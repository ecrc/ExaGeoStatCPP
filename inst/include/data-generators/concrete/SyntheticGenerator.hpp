
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

            class SyntheticGenerator : public DataGenerator {
            public:

                /**
                *  @brief
                *  Class constructor.
                *
                *  @param[in] apConfigurations
                *  Pointer to Synthetic data Configurations.
                 *
                */
                explicit SyntheticGenerator(
                        configurations::data_configurations::SyntheticDataConfigurations *apConfigurations);

                /**
                 * @brief
                 * Virtual destructor to allow calls to the correct concrete destructor.
                 *
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
                 * Generate uniform distribution between rangeLow , rangeHigh.
                 *
                 * @param[in] aRangeLow
                 * The Lower range.
                 *
                 * @param[in] aRangeHigh
                 * The Higher range.
                 *
                 * @return scaled_range
                 * The scaled uniform distribution between the two bounds .
                 *
                 */
                double UniformDistribution(double aRangeLow, double aRangeHigh);

                /**
                 * @brief
                 * Sort in Morton order (input points must be in [0;1]x[0;1] square]).
                 *
                 * @param[in] aN
                 * The problem size.
                 *
                 * @param[in] aLocations
                 * X, Y and Z variables.
                 *
                 */
                void SortLocations(int aN, dataunits::Locations aLocations);

                uint32_t EncodeMorton2(uint32_t x, uint32_t y);
                uint32_t Part1By1(uint32_t x);
                static int compare_uint32(const void *a, const void *b);
                uint32_t Compact1By1(uint32_t x);

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
                InitializeLocations(dataunits::Locations aLocations) override;

                void Print() override;

            private:


            };

        }//namespace Synthetic
    }//namespace generators
}//namespace exageostat

#endif //EXAGEOSTAT_CPP_SYNTHETICGENERATOR_HPP
