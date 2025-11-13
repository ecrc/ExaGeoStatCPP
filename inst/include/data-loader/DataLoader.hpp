
// Copyright (c) 2017-2024 King Abdullah University of Science and Technology,
// All rights reserved.
// ExaGeoStat is a software package, provided by King Abdullah University of Science and Technology (KAUST).

/**
 * @file DataLoader.hpp
 * @brief Manages data loading operations for ExaGeoStat.
 * @version 1.1.0
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
        std::unique_ptr<ExaGeoStatData<T>>
        CreateData(configurations::Configurations &aConfigurations, kernels::Kernel<T> &aKernel) override;

        /**
         * @brief Reads data from external sources into ExaGeoStat format.
         * @param aConfigurations Configuration settings for data loading.
         * @param aMeasurementsMatrix Vector to store measurement values.
         * @param aXLocations Vector to store X coordinates of locations.
         * @param aYLocations Vector to store Y coordinates of locations.
         * @param aZLocations Vector to store Z coordinates of locations (if applicable).
         * @param aP Partition index for distributed data loading.
         * @return void
         *
         */
        virtual void
        ReadData(configurations::Configurations &aConfigurations, std::vector<T> &aMeasurementsMatrix,
                 std::vector<T> &aXLocations,
                 std::vector<T> &aYLocations, std::vector<T> &aZLocations, const int &aP) = 0;

        /**
        * @brief Writes a matrix of vectors to disk.
        * @param[in] aMatrixPointer A Reference to the matrix data.
        * @param[in] aProblemSize The size of the problem.
        * @param[in] aP The number of processes.
        * @param[in] aLoggerPath The path to the logger file.
        * @param[in] aLocations A Reference to the Locations object.
        * @return void
        *
        */
        virtual void
        WriteData(const T &aMatrixPointer, const int &aProblemSize, const int &aP, std::string &aLoggerPath,
                  exageostat::dataunits::Locations<T> &aLocations) = 0;

        /**
         * @brief Abstract method for loading data based on provided configurations and kernel.
         * @param[in] aConfigurations Reference to the configurations object that contains parameters for loading data.
         * @param[in] aKernel Reference to the kernel object that defines the operations to be applied while loading the data.
         * @return A unique pointer to the loaded ExaGeoStatData object.
         *
         */
        virtual std::unique_ptr<ExaGeoStatData<T>>
        LoadData(configurations::Configurations &aConfigurations, exageostat::kernels::Kernel<T> &aKernel)  = 0;

        /**
         * @brief Factory method for creating a DataLoader instance based on the given configurations.
         * This method dynamically determines the type of data loader to instantiate based on compile-time conditions.
         * @param[in] aConfigurations Reference to the configurations object that contains parameters for loading data.
         * @return A unique pointer to a DataLoader instance configured as per the specified runtime conditions.
         *
         */
        static std::unique_ptr<DataLoader<T>>
        CreateDataLoader(exageostat::configurations::Configurations &apConfigurations);

        /**
         * @brief Releases the singleton instance of the currently active DataLoader.
         * This method ensures proper deallocation of the singleton instance of the data loader,
         * depending on the selected runtime.
         *
         */
        static void ReleaseDataLoader();

    };

    /**
     * @brief Instantiates the Synthetic Data Generator class for float and double types.
     * @tparam T Data Type: float or double
     *
     */
    EXAGEOSTAT_INSTANTIATE_CLASS(DataLoader)
} // namespace exageostat

#endif //EXAGEOSTATCPP_DATALOADER_HPP
