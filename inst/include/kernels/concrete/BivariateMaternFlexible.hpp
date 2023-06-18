
// Copyright (c) 2017-2023 King Abdullah University of Science and Technology,
// All rights reserved.
// ExaGeoStat is a software package, provided by King Abdullah University of Science and Technology (KAUST).

/**
 * @file BivariateMaternFlexible.hpp
 * @brief Defines the BivariateMaternFlexible class, a Bivariate Matern Flexible kernel.
 * @version 1.0.0
 * @author Sameh Abdulah
 * @date 2023-04-14
**/

#ifndef EXAGEOSTATCPP_BIVARIATEMATERNFLEXIBLE_HPP
#define EXAGEOSTATCPP_BIVARIATEMATERNFLEXIBLE_HPP

#include <kernels/Kernel.hpp>

namespace exageostat {
    namespace kernels {

        /**
         * @class BivariateMaternFlexible
         * @brief A class representing a Bivariate Matern Flexible kernel.
         * @details This class represents a Bivariate Matern Flexible, which is a subclass of the Kernel class.
         * It provides a method for generating a covariance matrix using a set of input locations and kernel parameters.
         *
         */
        class BivariateMaternFlexible : public Kernel {

        public:

            /**
             * @brief Constructs a new BivariateMaternFlexible object.
             * @details Initializes a new BivariateMaternFlexible object with default values.
             */
            BivariateMaternFlexible();

            /**
             * @brief Virtual destructor to allow calls to the correct concrete destructor.
             *
             */
            ~BivariateMaternFlexible() = default;

            /**
             * @brief Generates a covariance matrix using a set of locations and kernel parameters.
             * @copydoc Kernel::GenerateCovarianceMatrix()
             */
            void GenerateCovarianceMatrix(double *apMatrixA, int &aRowsNumber, int &aColumnsNumber, int &aRowOffset,
                                          int &aColumnOffset, dataunits::Locations *apLocation1,
                                          dataunits::Locations *apLocation2, dataunits::Locations *apLocation3,
                                          double *aLocalTheta, int &aDistanceMetric) override ;

            /**
             * @brief Creates a new BivariateMaternFlexible object.
             * @details This method creates a new BivariateMaternFlexible object and returns a pointer to it.
             * @return A pointer to the new BivariateMaternFlexible object.
             *
             */
            static Kernel *Create();

        private:
            //// Used plugin name for static registration
            static bool plugin_name;
        };
    }//namespace Kernels
}//namespace exageostat
#endif //EXAGEOSTATCPP_BIVARIATEMATERNFLEXIBLE_HPP
