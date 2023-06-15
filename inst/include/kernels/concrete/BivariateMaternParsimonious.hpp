
// Copyright (c) 2017-2023 King Abdullah University of Science and Technology,
// All rights reserved.
// ExaGeoStat is a software package, provided by King Abdullah University of Science and Technology (KAUST).

/**
 * @file BivariateMaternParsimonious.hpp
 * @brief Defines the BivariateMaternParsimonious class, a Bivariate Matern Parsimonious kernel.
 * @version 1.0.0
 * @author Sameh Abdulah
 * @date 2023-04-14
**/

#ifndef EXAGEOSTATCPP_BIVARIATEMATERNPARSIMONIOUS_HPP
#define EXAGEOSTATCPP_BIVARIATEMATERNPARSIMONIOUS_HPP

#include <kernels/Kernel.hpp>

namespace exageostat {
    namespace kernels {

        /**
        * @class BivariateMaternParsimonious
        * @brief A class representing a Bivariate Matern Parsimonious kernel.
        *
        * This class represents a Bivariate Matern Parsimonious kernel, which is a subclass of the Kernel class. It provides
        * a method for generating a covariance matrix using a set of input locations and kernel parameters.
        */
        class BivariateMaternParsimonious : public Kernel {
        public:

            /**
             * @brief Constructs a new BivariateMaternParsimonious object.
             *
             * Initializes a new BivariateMaternParsimonious object with default values.
             */
            BivariateMaternParsimonious();

            /**
             * @brief Virtual destructor to allow calls to the correct concrete destructor.
             */
            ~BivariateMaternParsimonious() = default;

            /**
             * @brief Generates a covariance matrix using a set of locations and kernel parameters.
             * @param[in] apMatrixA The output covariance matrix.
             * @param[in] aRowsNumber The number of rows in the output matrix.
             * @param[in] aColumnsNumber The number of columns in the output matrix.
             * @param[in] aRowOffset The row offset for the input locations.
             * @param[in] aColumnOffset The column offset for the input locations.
             * @param[in] apLocation1 The set of input locations 1.
             * @param[in] apLocation2 The set of input locations 2.
             * @param[in] apLocation3 The set of input locations 3.
             * @param[in] aLocalTheta An array of kernel parameters.
             * @param [in] aDistanceMetric Distance metric to be used (1 = Euclidean, 2 = Manhattan, 3 = Minkowski).
             */
            void GenerateCovarianceMatrix(double *apMatrixA, int &aRowsNumber, int &aColumnsNumber, int &aRowOffset,
                                          int &aColumnOffset, dataunits::Locations *apLocation1,
                                          dataunits::Locations *apLocation2, dataunits::Locations *apLocation3,
                                          double *aLocalTheta, int &aDistanceMetric) override;

            /**
             * @brief Creates a new BivariateMaternParsimonious object.
             * @return A pointer to the new BivariateMaternParsimonious object.
             *
             * This method creates a new BivariateMaternParsimonious object and returns a pointer to it.
             */
            static Kernel *Create();

        private:
            //// Used plugin name for static registration
            static bool plugin_name;
        };
    }//namespace Kernels
}//namespace exageostat

#endif //EXAGEOSTATCPP_BIVARIATEMATERNPARSIMONIOUS_HPP
