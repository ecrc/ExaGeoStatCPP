
// Copyright (c) 2017-2023 King Abdullah University of Science and Technology,
// Copyright (C) 2023 by Brightskies inc,
// All rights reserved.
// ExaGeoStat is a software package, provided by King Abdullah University of Science and Technology (KAUST).

/**
 * @file Locations.hpp
 *
 * @version 1.0.0
 * @author Sameh Abdulah
 * @date 2023-02-27
**/

#ifndef EXAGEOSTAT_CPP_LOCATIONS_HPP
#define EXAGEOSTAT_CPP_LOCATIONS_HPP
namespace exageostat {
    namespace dataunits {

        class Locations {
        public:
            /**
             * @brief Default constructor.
             */
            Locations() = default;

            /**
             * @brief Virtual destructor to allow calls to the correct concrete destructor.
             */
            virtual ~Locations() = default;

            /**
             * @brief
             * LocationX setter.
             *
             * @param[in] apLocationX
             * Pointer to X data.
             *
             */
            void
            SetLocationX(double *apLocationX);

            /**
             * @brief
             * LocationX getter.
             *
             * @return mpLocationX
             * Pointer to X data.
             *
             */
            double *
            GetLocationX();

            /**
             * @brief
             * LocationY setter.
             *
             * @param[in] apLocationY
             * Pointer to Y data.
             *
             */
            void
            SetLocationY(double *apLocationY);

            /**
             * @brief
             * LocationY getter.
             *
             * @return mpLocationY
             * Pointer to Y data.
             *
             */
            double *
            GetLocationY();

            /**
             * @brief
             * LocationZ setter.
             *
             * @param[in] apLocationZ
             * Pointer to Z data.
             *
             */
            void
            SetLocationZ(double *apLocationZ);

            /**
             * @brief
             * LocationZ getter.
             *
             * @return mpLocationZ
             * Pointer to Z data.
             *
             */
            double *
            GetLocationZ();

        private:
            /// Used Location X.
            double *mpLocationX;
            /// Used Location Y.
            double *mpLocationY;
            /// Used Location Z.
            double *mpLocationZ;
        };

    }//namespace configurations
}//namespace exageostat

#endif //EXAGEOSTAT_CPP_LOCATIONS_HPP
