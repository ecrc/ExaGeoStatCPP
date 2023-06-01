
// Copyright (c) 2017-2023 King Abdullah University of Science and Technology,
// Copyright (C) 2023 by Brightskies inc,
// All rights reserved.
// ExaGeoStat is a software package, provided by King Abdullah University of Science and Technology (KAUST).

/**
 * @file Locations.hpp
 * @brief Header file for the Locations class, which contains methods to set and get location data.
 *
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
         */
        class Locations {
        public:
            /**
             * @brief Default constructor.
             */
            Locations(int aSize, exageostat::common::Dimension aDimension);

            /**
             * @brief Virtual destructor to allow calls to the correct concrete destructor.
             */
            virtual ~Locations();

            /**
             * @brief Setter for LocationX.
             *
             * @param[in] apLocationX Pointer to X data.
             */
            void SetLocationX(double *apLocationX);

            /**
             * @brief Getter for LocationX.
             *
             * @return mpLocationX Pointer to X data.
             */
            double * GetLocationX();

            /**
             * @brief Setter for LocationY.
             *
             * @param[in] apLocationY Pointer to Y data.
             */
            void SetLocationY(double *apLocationY);

            /**
             * @brief Getter for LocationY.
             *
             * @return mpLocationY Pointer to Y data.
             */
            double * GetLocationY();

            /**
             * @brief Setter for LocationZ.
             *
             * @param[in] apLocationZ Pointer to Z data.
             */
            void SetLocationZ(double *apLocationZ);

            /**
             * @brief Getter for LocationZ.
             *
             * @return mpLocationZ Pointer to Z data.
             */
            double * GetLocationZ();

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