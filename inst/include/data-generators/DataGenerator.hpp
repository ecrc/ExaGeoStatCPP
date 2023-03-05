
/**
 * @file DataGenerator.hpp
 *
 * @version 1.0.0
 * @author Sameh Abdulah
 * @date 2023-02-14
**/

#ifndef EXAGEOSTAT_CPP_DATAGENERATOR_HPP
#define EXAGEOSTAT_CPP_DATAGENERATOR_HPP

namespace exageostat {
    namespace generators {

        class DataGenerator {
        public:

            /**
             * @brief Virtual destructor to allow calls to the correct concrete destructor.
             */
            virtual ~DataGenerator() = default;
        };
    }//namespace generators
}//namespace exageostat

#endif //EXAGEOSTAT_CPP_DATAGENERATOR_HPP
