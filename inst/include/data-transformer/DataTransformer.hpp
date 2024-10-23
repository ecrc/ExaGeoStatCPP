
// Copyright (c) 2017-2024 King Abdullah University of Science and Technology,
// All rights reserved.
// ExaGeoStat is a software package, provided by King Abdullah University of Science and Technology (KAUST).

/**
 * @file DataTransformer.hpp
 * @brief Contains the definition of the DataTransformer class.
 * @version 1.1.0
 * @author Mahmoud ElKarargy
 * @author Sameh Abdulah
 * @date 2024-10-15
**/

#ifndef EXAGEOSTATCPP_DATATRANSFORMER_HPP
#define EXAGEOSTATCPP_DATATRANSFORMER_HPP

#include <configurations/Configurations.hpp>
#include <data-units/ExaGeoStatData.hpp>

namespace exageostat::transformer{

    /**
     * @brief Class represents the data transformer for the Climate Emulator.
     * @tparam T Data Type: float or double
     */
    template<typename T>
    class DataTransformer {

        /**
         * @brief Performs the forward spherical harmonics transform (SHT).
         * @param[in] aConfigurations Configurations object containing relevant settings.
         * @param[in,out] aData Descriptor Data object to be populated with descriptors and data.
         */
        static void ForwardSHT(configurations::Configurations &aConfigurations, std::unique_ptr<ExaGeoStatData<T>> &aData);

        /**
         * @brief Reshapes data during the forward phase of the simulation.
         * @param[in] aConfigurations Configurations object containing relevant settings.
         * @param[in,out] aData Descriptor Data object to be populated with descriptors and data.
         */
        static void ForwardReshape(configurations::Configurations &aConfigurations, std::unique_ptr<ExaGeoStatData<T>> &aData);

        /**
         * @brief Performs the inverse spherical harmonics transform (SHT).
         * @param[in] aConfigurations Configurations object containing relevant settings.
         * @param[in,out] aData Descriptor Data object to be populated with descriptors and data.
         */
        static void InverseSHT(configurations::Configurations &aConfigurations, std::unique_ptr<ExaGeoStatData<T>> &aData);
    };
}
#endif // EXAGEOSTATCPP_DATATRANSFORMER_HPP
