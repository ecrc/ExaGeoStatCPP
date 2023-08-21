
// Copyright (c) 2017-2023 King Abdullah University of Science and Technology,
// All rights reserved.
// ExaGeoStat is a software package, provided by King Abdullah University of Science and Technology (KAUST).

/**
 * @file ExaGeoStatData.hpp
 * @brief Contains the definition of the ExaGeoStatData class.
 * @version 1.0.0
 * @author Sameh Abdulah
 * @date 2023-07-19
**/

#ifndef EXAGEOSTATCPP_EXAGEOSTATDATA_HPP
#define EXAGEOSTATCPP_EXAGEOSTATDATA_HPP

#include <data-units/DescriptorData.hpp>
#include <data-units/Locations.hpp>
#include <hardware/ExaGeoStatHardware.hpp>

namespace exageostat {
    namespace dataunits {

        /**
         * @Class ExaGeoStatData
         * @brief  Manages geo-statistical data with functions for location and descriptor manipulation
         * @tparam T Data Type: float or double
         */
        template<typename T>
        class ExaGeoStatData {

        public:
            /**
             * @brief Constructor for ExaGeoStatData.
             * @param[in] aSize The size of the data.
             * @param[in] aDimension The dimension of the data.
             * @param[in] apHardware Reference to the ExaGeoStatHardware object.
             */
            ExaGeoStatData(int aSize, exageostat::common::Dimension aDimension,
                           hardware::ExaGeoStatHardware &apHardware);

            /**
             * @brief Destructor for ExaGeoStatData.
             */
            ~ExaGeoStatData();

            /**
             * @brief Get the locations.
             * @return Pointer to the Locations object.
             */
            Locations<T> *GetLocations();

            /**
             * @brief Set the locations.
             * @param[in] aLocation Pointer to the Locations object.
             */
            void SetLocations(Locations<T> &aLocation);

            /**
             * @brief Get the descriptor data.
             * @return Pointer to the DescriptorData object.
             */
            DescriptorData<T> *GetDescriptorData();

            /**
             * @brief Setter for the number of performed MLE iterations.
             * @param[in] aMleIterations number of performed MLE iterations.
             * @return void
             */
            void SetMleIterations(int aMleIterations);

            /**
             * @brief Get the number of performed MLE iterations.
             * @return Pointer to the DescriptorData object.
             */
            int GetMleIterations();

            /**
             * @brief Calculates Median Locations.
             * @param[in] aKernelName Name of the Kernel used.
             * @param[out] apLocations Location object to save medianLocations in.
             * @return void
             */
            void CalculateMedianLocations(std::string &aKernelName, dataunits::Locations<T> &apLocations);

        private:
            //// Used descriptor data.
            DescriptorData<T> *mpDescriptorData = nullptr;
            //// Used locations data.
            Locations<T> *mpLocations = nullptr;
            //// Current number of performed MLE iterations.
            int mMleIterations = 0;
        };

        /**
         * @brief Instantiates the ExaGeoStatData class for float and double types.
         * @tparam T Data Type: float or double
         */
        EXAGEOSTAT_INSTANTIATE_CLASS(ExaGeoStatData)
    } // namespace dataunits
} // namespace exageostat

#endif //EXAGEOSTATCPP_EXAGEOSTATDATA_HPP