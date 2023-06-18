
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
            void SetLocationX(double *apLocationX);

            /**
             * @brief Getter for LocationX.
             * @return Pointer to X data.
             *
             */
            double *GetLocationX();

            /**
             * @brief Setter for LocationY.
             * @param[in] apLocationY Pointer to Y data.
             * @return void
             *
             */
            void SetLocationY(double *apLocationY);

            /**
             * @brief Getter for LocationY.
             * @return Pointer to Y data.
             *
             */
            double *GetLocationY();

            /**
             * @brief Setter for LocationZ.
             * @param[in] apLocationZ Pointer to Z data.
             * @return void
             *
             */
            void SetLocationZ(double *apLocationZ);

            /**
             * @brief Getter for LocationZ.
             * @return Pointer to Z data.
             *
             */
            double *GetLocationZ();

        private:
            /// Pointer to X data.
            double *mpLocationX = nullptr;
            /// Pointer to Y data.
            double *mpLocationY = nullptr;
            /// Pointer to Z data.
            double *mpLocationZ = nullptr;
        };

    }//namespace configurations
}//namespace exageostat

#endif //EXAGEOSTAT_CPP_LOCATIONS_HPP