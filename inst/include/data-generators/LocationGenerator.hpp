
// Copyright (c) 2017-2023 King Abdullah University of Science and Technology,
// All rights reserved.
// ExaGeoStat is a software package, provided by King Abdullah University of Science and Technology (KAUST).

/**
 * @file LocationGenerator.hpp
 * @brief Generates and manages spatial locations for ExaGeoStat.
 * @version 1.1.0
 * @author Mahmoud ElKarargy
 * @date 2024-02-04
**/

#ifndef EXAGEOSTATCPP_LOCATIONGENERATOR_HPP
#define EXAGEOSTATCPP_LOCATIONGENERATOR_HPP

#include <data-units/Locations.hpp>

namespace exageostat::generators {

    /**
     * @class LocationGenerator
     * @brief Generates spatial locations based on given parameters.
     * @tparam T Data Type: float or double
     *
     */
    template<typename T>
    class LocationGenerator {

    public:

        /**
         * @brief Generates the data locations.
         * @details This method generates the X, Y, and Z variables used to define the locations of the data points.
         * @param[in] aN The number of data points.
         * @param[in] aTimeSlot The time slot.
         * @param[in] aDimension The dimension of the locations.
         * @param[out] aLocations Reference to the Locations object where the generated data will be stored.
         * @return void
         *
         */
        static void GenerateLocations(const int &aN, const int &aTimeSlot, const common::Dimension &aDimension, dataunits::Locations <T> &aLocations);

        /**
         * @brief Generate uniform distribution between rangeLow , rangeHigh.
         * @param[in] aRangeLow The Lower range.
         * @param[in] aRangeHigh The Higher range.
         * @return The scaled uniform distribution between the two bounds.
         *
         */
        static T UniformDistribution(const T &aRangeLow, const T &aRangeHigh);

        /**
         * @brief Sort locations in Morton order (input points must be in [0;1]x[0;1] square]).
         * @param[in] aN The problem size divided by P-Grid.
         * @param[in] aDimension Dimension of locations.
         * @param[in,out] aLocations Locations to be sorted.
         * @return void
         *
         */
        static void SortLocations(const int &aN, const common::Dimension &aDimension, dataunits::Locations<T> &aLocations);
    };

    /**
     * @brief Instantiates the Data Generator class for float and double types.
     * @tparam T Data Type: float or double
     *
     */
    EXAGEOSTAT_INSTANTIATE_CLASS(LocationGenerator)
}//namespace exageostat

#endif //EXAGEOSTATCPP_LOCATIONGENERATOR_HPP
