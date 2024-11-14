
// Copyright (c) 2017-2024 King Abdullah University of Science and Technology,
// All rights reserved.
// ExaGeoStat is a software package, provided by King Abdullah University of Science and Technology (KAUST).

/**
 * @file  ParsecLoader.hpp
 * @brief A class for loading PaRSEC format data.
 * @version 1.1.0
 * @author Mahmoud ElKarargy
 * @author Sameh Abdulah
 * @date 2024-02-04
**/

#ifndef EXAGEOSTAT_CPP_PARSECDATALOADER_HPP
#define EXAGEOSTAT_CPP_PARSECDATALOADER_HPP

#include <data-loader/DataLoader.hpp>

namespace exageostat::dataLoader::parsec {

    /**
     * @class  ParsecLoader
     * @brief A class for creating data by reading PaRSEC files.
     * @tparam T Data Type: float or double
     */
    template<typename T>
    class ParsecLoader : public DataLoader<T> {
    public:

        /**
         * @brief Get a pointer to the singleton instance of the  ParsecLoader class.
         * @return A pointer to the instance of the  ParsecLoader class.
         *
         */
        static ParsecLoader<T> *GetInstance();

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
         * @brief Creates the data by synthetically generating it.
         * @copydoc DataGenerator::LoadData()
         *
         */
        std::unique_ptr<ExaGeoStatData<T>>
        LoadData(configurations::Configurations &aConfigurations, kernels::Kernel<T> &aKernel) override;

        /**
         * @brief Release the singleton instance of the  ParsecLoader class.
         * @return void
         *
         */
        static void ReleaseInstance();

        /**
         * @brief Reads data from a CSV file into a matrix.
         * @param[in] apFilename Name of the CSV file.
         * @param[out] apFileContent Pointer to an array where file contents will be stored.
         * @param[in] aM Number of rows in the matrix.
         * @param[in] aN Number of columns in the matrix.
         * @return 0 on success, or a non-zero error code on failure.
         *
         */
        int ReadCSVFileHelper(const char* apFilename, double *apFileContent, int aM, int aN);

        /**
         * @brief Adapter for the matrix compress operation
         * @param[in] aConfigurations Configurations object containing relevant settings.
         * @param[in,out] aData Descriptor Data object to be populated with descriptors and data.
         * @return void.
         *
         */
        void CompressMatrixHelper(configurations::Configurations &aConfigurations, std::unique_ptr<ExaGeoStatData<T>> &aData);


    private:
        /**
         * @brief Constructor for the  ParsecLoader class.
         * @return void
         *
         */
        ParsecLoader() = default;

        /**
         * @brief Default destructor.
         *
         */
        ~ParsecLoader() override = default;

        /**
         * @brief Pointer to the singleton instance of the  ParsecLoader class.
         *
         */
        static ParsecLoader<T> *mpInstance;

    };

    /**
     * @brief Instantiates the PaRSEC Data Generator class for float and double types.
     * @tparam T Data Type: float or double
     *
     */
    EXAGEOSTAT_INSTANTIATE_CLASS(ParsecLoader)
}
#endif //EXAGEOSTAT_CPP_PARSECDATALOADER_HPP