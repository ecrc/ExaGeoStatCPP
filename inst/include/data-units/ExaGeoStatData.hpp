/**
 * @file ExaGeoStatData.hpp
 * @brief 
 * @version 1.0.0
 * @author Sameh Abdulah
 * @date 2023-07-19
**/

#ifndef EXAGEOSTATCPP_EXAGEOSTATDATA_HPP
#define EXAGEOSTATCPP_EXAGEOSTATDATA_HPP

#include <data-units/DescriptorData.hpp>
#include <data-units/Locations.hpp>

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
             * @param aSize Locations Size.
             * @param aDimension Dimensions used.
             */
            ExaGeoStatData(int aSize, exageostat::common::Dimension aDimension);

            /**
             * Destructor for ExaGeoStatData.
             */
            ~ExaGeoStatData();

            /**
             * @brief Getter for Locations.
             * @return Locations
             */
            Locations<T> *GetLocations();

            /**
             * @brief Setter for mpLocations.
             * @param apLocation
             * @return void
             */
            void SetLocations(Locations<T> *apLocation);

            /**
             * @brief Getter for DescriptorData.
             * @return DescriptorData
             */
            DescriptorData<T> *GetDescriptorData();

            /**
             * @brief Calculates Median Locations.
             * @param[in] aKernelName Name of the Kernel used.
             * @return median locations
             */
            Locations<T> *CalculateMedianLocations(std::string &aKernelName);

        private:
            /// Used Descriptors' Data.
            DescriptorData<T> *mpDescriptorData = nullptr;
            /// Used Locations.
            Locations<T> *mpLocations;
        };


        /**
        * @brief Instantiates the chameleon dense class for float and double types.
        * @tparam T Data Type: float or double
        *
        */
        EXAGEOSTAT_INSTANTIATE_CLASS(ExaGeoStatData)
    }//namespace configurations
}//namespace exageostat

#endif //EXAGEOSTATCPP_EXAGEOSTATDATA_HPP
