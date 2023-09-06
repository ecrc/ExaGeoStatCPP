
// Copyright (c) 2017-2023 King Abdullah University of Science and Technology,
// All rights reserved.
// ExaGeoStat is a software package, provided by King Abdullah University of Science and Technology (KAUST).

/**
 * @file Locations.hpp
 * @brief Header file for the Locations class, which contains methods to set and get location data.
 * @version 1.0.0
 * @author Mahmoud ElKarargy
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
         * @tparam T Data Type: float or double
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
            Locations(const int &aSize, const exageostat::common::Dimension &aDimension);

            /**
             * @brief Default copy constructor.
             * @param[in] aLocations Locations to be copied.
             */
            Locations(const Locations<T> &aLocations) = default;

            /**
             * @brief Virtual destructor to allow calls to the correct concrete destructor.
             *
             */
            virtual ~Locations();

            /**
             * @brief Setter for LocationX.
             * @param[in] aLocationX Reference to X data.
             * @return void
             *
             */
            void SetLocationX(T &aLocationX);

            /**
             * @brief Getter for LocationX.
             * @return Pointer to X data.
             *
             */
            T *GetLocationX();

            /**
             * @brief Setter for LocationY.
             * @param[in] aLocationY Reference to Y data.
             * @return void
             *
             */
            void SetLocationY(T &aLocationY);

            /**
             * @brief Getter for LocationY.
             * @return Pointer to Y data.
             *
             */
            T *GetLocationY();

            /**
             * @brief Setter for LocationZ.
             * @param[in] aLocationZ Reference to Z data.
             * @return void
             *
             */
            void SetLocationZ(T &aLocationZ);

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
            void SetSize(const int &aSize);

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
            void SetDimension(const common::Dimension &aDimension);

            /**
             * @brief Getter for Dimension.
             * @return Locations dimension.
             *
             */
            common::Dimension GetDimension();

        private:
            /// Pointer to X data.
            T *mpLocationX = nullptr;
            /// Pointer to Y data.
            T *mpLocationY = nullptr;
            /// Pointer to Z data.
            T *mpLocationZ = nullptr;
            /// Size of each dimension
            int mSize = 1;
            /// Data dimensions
            common::Dimension mDimension = common::Dimension2D;
        };

        /**
        * @brief Instantiates the Linear Algebra methods class for float and double types.
        * @tparam T Data Type: float or double
        *
        */
        EXAGEOSTAT_INSTANTIATE_CLASS(Locations)
    }//namespace dataunits
}//namespace exageostat

#endif //EXAGEOSTAT_CPP_LOCATIONS_HPP