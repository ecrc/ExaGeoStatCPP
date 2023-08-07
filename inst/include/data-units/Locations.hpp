
// Copyright (c) 2017-2023 King Abdullah University of Science and Technology,
// All rights reserved.
// ExaGeoStat is a software package, provided by King Abdullah University of Science and Technology (KAUST).

/**
 * @file Locations.hpp
 * @brief Header file for the Locations class, which contains methods to set and get location data.
 * @version 1.0.0
 * @author Sameh Abdulah
 * @date 2023-02-27
**/

#ifndef EXAGEOSTAT_CPP_LOCATIONS_HPP
#define EXAGEOSTAT_CPP_LOCATIONS_HPP

#include <common/Definitions.hpp>

namespace exageostat {
    namespace dataunits {

        /**
         * @class Locations
         * @brief A class containing methods to set and get location data.
         *
         */
        template<typename T>
        class Locations {
        public:
            /**
             * @brief Constructor.
             * @param[in] aSize The number of data points.
             * @param[in] aDimension The dimensionality of the data points.
             * @return void
             *
             */
            Locations(int aSize, exageostat::common::Dimension aDimension);

            Locations(const Locations<T> &lo) = default;
//            // define move assignment operator
//            Locations& operator=(Locations&& other) noexcept {
//                // move all fields from other to this object
//                mpLocationX = std::move(other.mpLocationX);
//                mpLocationY = std::move(other.mpLocationY);
//                // ...
//                return *this;
//            }

            /**
             * @brief Virtual destructor to allow calls to the correct concrete destructor.
             *
             */
            virtual ~Locations();

            /**
             * @brief Setter for LocationX.
             * @param[in] apLocationX Pointer to X data.
             * @return void
             *
             */
            void SetLocationX(T *apLocationX);

            /**
             * @brief Getter for LocationX.
             * @return Pointer to X data.
             *
             */
            T *GetLocationX();

            /**
             * @brief Setter for LocationY.
             * @param[in] apLocationY Pointer to Y data.
             * @return void
             *
             */
            void SetLocationY(T *apLocationY);

            /**
             * @brief Getter for LocationY.
             * @return Pointer to Y data.
             *
             */
            T *GetLocationY();

            /**
             * @brief Setter for LocationZ.
             * @param[in] apLocationZ Pointer to Z data.
             * @return void
             *
             */
            void SetLocationZ(T *apLocationZ);

            /**
             * @brief Getter for LocationZ.
             * @return Pointer to Z data.
             *
             */
            T *GetLocationZ();

            /**
             * @brief Setter for mSize.
             * @param[in] aSize.
             * @return void
             *
             */
            void SetSize(int aSize);

            /**
             * @brief Getter for mSize.
             * @return Locations size.
             *
             */
            int GetSize();

            /**
             * @brief Setter for Dimensions.
             * @param[in] aDimension.
             * @return void
             *
             */
            void SetDimension(common::Dimension aDimension);

            /**
             * @brief Getter for Dimension.
             * @return Locations dimension.
             *
             */
            common::Dimension GetDimension();

        public:
        private:
            /// Pointer to X data.
            T *mpLocationX = nullptr;
            /// Pointer to Y data.
            T *mpLocationY = nullptr;
            /// Pointer to Z data.
            T *mpLocationZ = nullptr;
            /// Size of each dimension
            int mSize;
            /// Data dimensions
            common::Dimension mDimension;
        };

        /**
        * @brief Instantiates the Linear Algebra methods class for float and double types.
        * @tparam T Data Type: float or double
        *
        */
        EXAGEOSTAT_INSTANTIATE_CLASS(Locations)
    }//namespace configurations
}//namespace exageostat

#endif //EXAGEOSTAT_CPP_LOCATIONS_HPP