
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
                void GenerateLocations(int aN, int aTimeSlots);

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
                void SortLocations(int aN);

                /**
                 * @brief
                 * Spread bits one by one.
                 *
                 * @param[in] aInputByte
                 * The input byte to be spread.
                 *
                 * @return aInputByte
                 * Returns Input byte After being spread.
                 *
                 */
                uint64_t SpreadBits(uint64_t aInputByte);

                /**
                 * @brief
                 * Spread bits one by one.
                 *
                 * @param[in] aInputByte
                 * The input spreaded byte to be compacted.
                 *
                 * @return aInputByte
                 * Returns Input byte After being compacted.
                 *
                 */
                uint64_t ReverseSpreadBits(uint64_t aInputByte);

                /**
                 * @brief
                 * Compares two Unit32 values
                 *
                 * @param[in] apFirstValue
                 * Pointer to the first input Unit32 value
                 *
                 * @param[in] apSecondValue
                 * Pointer to the second input Unit32 value
                 *
                 * @return value
                 *  0  in case of equality,
                 *  1  in case of apFirstValue bigger than apSecondValue
                 *  -1 in case of apSecondValue bigger than apFirstValue
                 *
                 */
                static bool CompareUint64(const uint64_t &aFirstValue, const uint64_t &aSecondValue);

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
                void
                InitializeLocations(int aTimeSlots = 1) override;

                void Print() override;

            private:


            };

        }//namespace Synthetic
    }//namespace generators
}//namespace exageostat

#endif //EXAGEOSTAT_CPP_SYNTHETICGENERATOR_HPP
