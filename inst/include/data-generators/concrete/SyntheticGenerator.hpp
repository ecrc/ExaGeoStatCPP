
/**
 * @file SyntheticGenerator.hpp
 * @brief A class for generating synthetic data.
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

            /**
             * @class SyntheticGenerator
             * @brief A class for generating synthetic data.
             */
            template<typename T>
            class SyntheticGenerator : public DataGenerator<T> {

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
                ~SyntheticGenerator() override = default;

                static std::vector<double> InitTheta(std::vector<double> apTheta, int size);

                /**
                 * @brief
                 * Generates the data locations.
                 * This method generates the X, Y, and Z variables used to define the locations of the data points.
                 *
                 */
                void
                GenerateLocations() override;

                /**
                 * @brief
                 * Generates the data descriptors.
                 * This method generates the descriptors used to define the properties of the data points.
                 *
                 */
                void GenerateDescriptors() override;

                /**
                 * @brief
                 * Generates the data observations.
                 *
                 * This method generates the observations of the data points, which are used to train and test the model.
                 *
                 * @return void
                 */
                void GenerateObservations() override;

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
                static double UniformDistribution(double aRangeLow, double aRangeHigh);

                /**
                 * @brief
                 * Sort in Morton order (input points must be in [0;1]x[0;1] square]).
                 *
                 * @param[in] aN
                 * The problem size divided by P-Grid.
                 *
                 */
                void SortLocations(int aN);

                /**
                 * @brief
                 * Spread bits by three spaces.
                 *
                 * @param[in] aInputByte
                 * The input 64 bit to be spread.
                 *
                 * @returns aInputByte
                 * The byte after being spread.
                 *
                 */
                static uint64_t SpreadBits(uint64_t aInputByte);

                /**
                 * @brief
                 * Reverse Spread bits operation.
                 *
                 * @param[in] aInputByte
                 *  The input spread 64 bit to be compacted.
                 *
                 * @returns aInputByte
                 * The byte after being compacted.
                 *
                 */
                static uint64_t ReverseSpreadBits(uint64_t aInputByte);

                /**
                 * @brief
                 * Compares two Unit64 values
                 *
                 * @param[in] aFirstValue
                 * Constant reference to the first input 64 bit value.
                 *
                 * @param[in] aSecondValue
                 * Constant reference to the second input 64 bit value.
                 *
                 * @return boolean
                 *  True in case of second value bigger than first value.
                 *  False otherwise.
                 *
                 */
                static bool CompareUint64(const uint64_t &aFirstValue, const uint64_t &aSecondValue);

            };

            /**
             * @brief Instantiates the LinearAlgebraMethods class for float and double types.
             */
            EXAGEOSTAT_INSTANTIATE_CLASS(SyntheticGenerator)
        }//namespace Synthetic
    }//namespace generators
}//namespace exageostat

#endif //EXAGEOSTAT_CPP_SYNTHETICGENERATOR_HPP