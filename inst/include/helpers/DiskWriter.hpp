// Copyright (c) 2017-2023 King Abdullah University of Science and Technology,
// All rights reserved.
// ExaGeoStat is a software package, provided by King Abdullah University of Science and Technology (KAUST).

/**
 * @file DiskWriter.hpp
 * @brief 
 * @version 1.0.0
 * @author Sameh Abdulah
 * @date 2023-06-08
**/

#ifndef EXAGEOSTATCPP_DISKWRITER_HPP
#define EXAGEOSTATCPP_DISKWRITER_HPP

#include <iostream>
#include <common/Definitions.hpp>
#include <data-units/Locations.hpp>

namespace exageostat {
    namespace helpers {

        /**
         * @class DiskWriter
         * @brief
         * @tparam T The data type of the disk writter.
         */
        template<typename T>
        class DiskWriter {
        public:
            void static WriteVectorsToDisk(T *apMatrixPointer, const int *apProblemSize, const int *apP, std::string *apLoggerPath, exageostat::dataunits::Locations *apLocations);
        };

        /**
         * @brief Instantiates the DiskWriter class for float and double types.
         */
        EXAGEOSTAT_INSTANTIATE_CLASS(DiskWriter)
    }
}
#endif //EXAGEOSTATCPP_DISKWRITER_HPP
