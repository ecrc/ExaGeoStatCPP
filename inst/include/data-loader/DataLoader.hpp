
// Copyright (c) 2017-2023 King Abdullah University of Science and Technology,
// All rights reserved.
// ExaGeoStat is a software package, provided by King Abdullah University of Science and Technology (KAUST).

/**
 * @file DataLoader.hpp
 * @brief Manages data loading operations for ExaGeoStat.
 * @version 1.0.1
 * @author Mahmoud ElKarargy
 * @author Sameh Abdulah
 * @date 2024-02-04
**/

#ifndef EXAGEOSTATCPP_DATALOADER_HPP
#define EXAGEOSTATCPP_DATALOADER_HPP

#include <data-generators/DataGenerator.hpp>

namespace exageostat::dataLoader {

    /**
     * @class DataLoader
     * @brief Extends DataGenerator to include data loading functionalities.
     * @tparam T Data Type: float or double
     *
     */
    template<typename T>
    class DataLoader : public generators::DataGenerator<T> {

    public:

        /**
         * @brief Creates the data by synthetically generating it.
         * @copydoc DataGenerator::CreateData()
         *
         */
        std::unique_ptr<dataunits::ExaGeoStatData<T>>
        CreateData(exageostat::configurations::Configurations &aConfigurations,
                   const exageostat::hardware::ExaGeoStatHardware &aHardware,
                   exageostat::kernels::Kernel<T> &aKernel) override;

        /**
         * @brief Reads data from external sources into ExaGeoStat format.
         * @param aConfigurations Configuration settings for data loading.
         * @param aMeasurementsMatrix Vector to store measurement values.
         * @param aXLocations Vector to store X coordinates of locations.
         * @param aYLocations Vector to store Y coordinates of locations.
         * @param aZLocations Vector to store Z coordinates of locations (if applicable).
         * @param aP Partition index for distributed data loading.
         */
        virtual void
        ReadData(exageostat::configurations::Configurations &aConfigurations, std::vector<T> &aMeasurementsMatrix,
                 std::vector<T> &aXLocations, std::vector<T> &aYLocations, std::vector<T> &aZLocations,
                 const int &aP) = 0;

    };

    /**
     * @brief Instantiates the Synthetic Data Generator class for float and double types.
     * @tparam T Data Type: float or double
     *
     */
    EXAGEOSTAT_INSTANTIATE_CLASS(DataLoader)
} // namespace exageostat

#endif //EXAGEOSTATCPP_DATALOADER_HPP
