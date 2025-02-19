
// Copyright (c) 2017-2024 King Abdullah University of Science and Technology,
// All rights reserved.
// ExaGeoStat is a software package, provided by King Abdullah University of Science and Technology (KAUST).

/**
 * @file DataAnalyzer.hpp
 * @brief Contains the definition of the DataAnalyzer class.
 * @version 1.1.0
 * @author Mahmoud ElKarargy
 * @author Sameh Abdulah
 * @date 2024-10-15
**/

#ifndef EXAGEOSTATCPP_DATAANALYZER_HPP
#define EXAGEOSTATCPP_DATAANALYZER_HPP

#include <configurations/Configurations.hpp>
#include <data-units/ExaGeoStatData.hpp>

namespace exageostat::analyzer{
    /**
     * @brief Class represents the data analyzer for the Climate Emulator.
     * @tparam T Data Type: float or double
     */
    template<typename T>
    class DataAnalyzer {

    public:

        /**
         * @brief Analyzes the given matrix data pre computation.
         * @param[in, out] aData Reference to an ExaGeoStatData object that contains matrix to be analyzed.
         * @return void
         *
         */
        static void PreAnalyzeMatrix(std::unique_ptr<ExaGeoStatData<T>> &aData);

        /**
         * @brief Analyzes the given matrix data post computation.
         * @param[in, out] aData Reference to an ExaGeoStatData object that contains matrix to be analyzed.
         * @return void
         *
         */
        static void PostAnalyzeMatrix(std::unique_ptr<ExaGeoStatData<T>> &aData);

        /**
         * @brief Compares between two matrices by getting the difference.
         * @param[in, out] aData Reference to an ExaGeoStatData object that contains matrix to be analyzed.
         * @return the calculated MSE.
         *
         */
        static double CompareMatDifference(std::unique_ptr<ExaGeoStatData<T>> &aData);
    };

    /**
     * @brief Instantiates the ExaGeoStat class for float and double types.
     * @tparam T Data Type: float or double
     *
     */
     EXAGEOSTAT_INSTANTIATE_CLASS(DataAnalyzer)

}//namespace exageostat

#endif // EXAGEOSTATCPP_DATAANALYZER_HPP