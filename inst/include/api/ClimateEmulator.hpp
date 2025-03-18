
// Copyright (c) 2017-2024 King Abdullah University of Science and Technology,
// All rights reserved.
// ExaGeoStat is a software package, provided by King Abdullah University of Science and Technology (KAUST).

/**
 * @file ClimateEmulator.hpp
 * @brief High-Level Wrapper class containing the static API for ClimateEmulator operations.
 * @version 1.1.0
 * @author Mahmoud ElKarargy
 * @date 2024-02-04
**/

#ifndef EXAGEOSTATCPP_CLIMATEEMULATOR_HPP
#define EXAGEOSTATCPP_CLIMATEEMULATOR_HPP

#include <configurations/Configurations.hpp>
#include <data-units/ExaGeoStatData.hpp>

namespace exageostat::api {
    /**
     * @brief High-Level Wrapper class containing the static API for ClimateEmulator operations.
     * @tparam T Data Type: float or double
     *
     */
    template<typename T>
    class ClimateEmulator {
    public:

        /**
         * @brief Generates Data for climate emulator by loading real data.
         * @param[in] aConfigurations Reference to Configurations object containing user input data.
         * @param[out] aData Reference to an ExaGeoStatData<T> object where generated data will be stored.
         * @return void
         *
         */
        static void
        ClimateEmulatorLoadData(configurations::Configurations &aConfigurations,
                                std::unique_ptr <ExaGeoStatData<T>> &aData);

        /**
         * @brief Models loaded data.
         * @param[in] aConfigurations Reference to Configurations object containing user input data.
         * @param[in] aData Reference to an ExaGeoStatData<T> object containing needed descriptors, and locations.
         * @param[in] apMeasurementsMatrix Pointer to the user input measurements matrix.
         * @return the last optimum value of MLE.
         *
         */
        static void
        ClimateEmulatorModeling(configurations::Configurations &aConfigurations,
                                std::unique_ptr <ExaGeoStatData<T>> &aData, T *apMeasurementsMatrix = nullptr);
    };

    /**
     * @brief Instantiates the ClimateEmulator class for float and double types.
     * @tparam T Data Type: float or double
     *
     */
    EXAGEOSTAT_INSTANTIATE_CLASS(ClimateEmulator)
}

#endif //EXAGEOSTATCPP_CLIMATEEMULATOR_HPP