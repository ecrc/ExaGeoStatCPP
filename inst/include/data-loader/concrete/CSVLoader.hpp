
// Copyright (c) 2017-2024 King Abdullah University of Science and Technology,
// All rights reserved.
// ExaGeoStat is a software package, provided by King Abdullah University of Science and Technology (KAUST).

/**
 * @file  CSVLoader.hpp
 * @brief A class for loading csv format data.
 * @version 1.1.0
 * @author Mahmoud ElKarargy
 * @author Sameh Abdulah
 * @date 2024-02-04
**/

#ifndef EXAGEOSTAT_CPP_CSVDATALOADER_HPP
#define EXAGEOSTAT_CPP_CSVDATALOADER_HPP

#include <data-loader/DataLoader.hpp>

namespace exageostat::dataLoader::csv {

    /**
     * @class  CSVLoader
     * @brief A class for creating data by reading CSV files.
     * @tparam T Data Type: float or double
     */
    template<typename T>
    class CSVLoader : public DataLoader<T> {
    public:

        /**
         * @brief Get a pointer to the singleton instance of the  CSVLoader class.
         * @return A pointer to the instance of the  CSVLoader class.
         *
         */
        static CSVLoader<T> *GetInstance();

        /**
         * @brief Reads data from external sources into ExaGeoStat format.
         * @copydoc DataLoader::ReadData()
         *
         */
        void ReadData(configurations::Configurations &aConfigurations, std::vector<T> &aMeasurementsMatrix,
                      std::vector<T> &aXLocations, std::vector<T> &aYLocations, std::vector<T> &aZLocations,
                      const int &aP) override;

        /**
        * @brief Writes a matrix of vectors to disk.
        * @copydoc DataLoader::WriteData()
        *
        */
        void
        WriteData(const T &aMatrixPointer, const int &aProblemSize, const int &aP, std::string &aLoggerPath,
                  exageostat::dataunits::Locations<T> &aLocations) override;

        /**
        * @brief Loads data based on given configuration.
        * @copydoc DataLoader::LoadData()
        *
        */
       std::unique_ptr<ExaGeoStatData<T>>
       LoadData(configurations::Configurations &aConfigurations, exageostat::kernels::Kernel<T> &aKernel) override;

       /**
         * @brief Release the singleton instance of the  CSVLoader class.
         * @return void
         *
         */
        static void ReleaseInstance();

    private:
        /**
         * @brief Constructor for the  CSVLoader class.
         * @return void
         *
         */
        CSVLoader() = default;

        /**
         * @brief Default destructor.
         *
         */
        ~CSVLoader() override = default;

        /**
         * @brief Pointer to the singleton instance of the  CSVLoader class.
         *
         */
        static CSVLoader<T> *mpInstance;

    };

    /**
     * @brief Instantiates the CSV Data Generator class for float and double types.
     * @tparam T Data Type: float or double
     *
     */
    EXAGEOSTAT_INSTANTIATE_CLASS(CSVLoader)
}
#endif //EXAGEOSTAT_CPP_CSVDATALOADER_HPP