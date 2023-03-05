
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
             * @brief LocationX setter.
             * @param apLocationX
             */
            void
            SetLocationX(double *apLocationX);
            /**
             * @brief LocationX getter.
             * @return mpLocationX
             */
            double*
            GetLocationX();
            /**
             * @brief LocationY setter.
             * @param apLocationY
             */
            void
            SetLocationY(double *apLocationY);
            /**
             * @brief LocationY getter.
             * @return mpLocationY
             */
            double*
            GetLocationY();
            /**
             * @brief LocationZ setter.
             * @param apLocationZ
             */
            void
            SetLocationZ(double *apLocationZ);
            /**
             * @brief LocationZ getter.
             * @return mpLocationZ
             */
            double*
            GetLocationZ();
            /**
             * @brief Check Dimension value.
             * @param aDimensions
             */

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
