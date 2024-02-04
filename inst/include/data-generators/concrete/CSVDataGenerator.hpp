
// Copyright (c) 2017-2023 King Abdullah University of Science and Technology,
// All rights reserved.
// ExaGeoStat is a software package, provided by King Abdullah University of Science and Technology (KAUST).

/**
 * @file CSVDataGenerator.hpp
 * @brief A class for generating synthetic data.
 * @version 1.1.0
 * @author Mahmoud ElKarargy
 * @author Sameh Abdulah
 * @date 2023-02-14
**/

#ifndef EXAGEOSTAT_CPP_CSVDATAGENERATOR_HPP
#define EXAGEOSTAT_CPP_CSVDATAGENERATOR_HPP

#include <data-generators/DataGenerator.hpp>

namespace exageostat::generators::csv {

    /**
     * @class CSVDataGenerator
     * @brief A class for creating data by reading CSV files.
     * @tparam T Data Type: float or double
     */
    template<typename T>
    class CSVDataGenerator : public DataGenerator<T> {
    public:

        /**
         * @brief Get a pointer to the singleton instance of the CSVDataGenerator class.
         * @return A pointer to the instance of the CSVDataGenerator class.
         *
         */
        static CSVDataGenerator<T> *GetInstance();

        /**
         * @brief Creates the data by reading a CSV file.
         * @copydoc DataGenerator::CreateData()
         *
         */
        std::unique_ptr<dataunits::ExaGeoStatData<T>> CreateData(Configurations &,
                                                                 const ExaGeoStatHardware &aHardware,
                                                                 exageostat::kernels::Kernel<T> &aKernel) override;

        /**
         * @brief Reads CSV files containing 2D, 3D, or ST locations, and measurements vector.
         * @param aConfigurations Reference to the Configurations object.
         * @param aMeasurementsMatrix Reference to the Measurement matrix to be filled with read data.
         * @param aXLocations Reference to the location's x coordinates matrix to be filled with read data.
         * @param aYLocations Reference to the location's y coordinates matrix to be filled with read data.
         * @param aZLocations Reference to the location's z coordinates matrix to be filled with read data.
         * @param[in] aP the P value of the kernel multiplied by time slot.
         * @return void
         */
        void ReadData(Configurations &aConfigurations, std::vector<T> &aMeasurementsMatrix,
                      std::vector<T> &aXLocations, std::vector<T> &aYLocations, std::vector<T> &aZLocations,
                      const int &aP);

        /**
         * @brief Release the singleton instance of the CSVDataGenerator class.
         * @return void
         *
         */
        static void ReleaseInstance();

    private:
        /**
         * @brief Constructor for the CSVDataGenerator class.
         * @return void
         *
         */
        CSVDataGenerator() = default;

        /**
         * @brief Default destructor.
         *
         */
        ~CSVDataGenerator() override = default;

        /**
         * @brief Pointer to the singleton instance of the CSVDataGenerator class.
         *
         */
        static CSVDataGenerator<T> *mpInstance;

    };

    /**
     * @brief Instantiates the CSV Data Generator class for float and double types.
     * @tparam T Data Type: float or double
     *
     */
    EXAGEOSTAT_INSTANTIATE_CLASS(CSVDataGenerator)
}
#endif //EXAGEOSTAT_CPP_CSVDATAGENERATOR_HPP