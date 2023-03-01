
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
                /**
                 * @brief Default constructor.
                 */
                SyntheticGenerator() = default;

                /**
                 * @brief Virtual destructor to allow calls to the correct concrete destructor.
                 */
                virtual ~SyntheticGenerator() = default;

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
                void GenerateLocations(int aN, int aSeed);

            };

        }//namespace Synthetic
    }//namespace generators
}//namespace exageostat

#endif //EXAGEOSTAT_CPP_SYNTHETICGENERATOR_HPP
