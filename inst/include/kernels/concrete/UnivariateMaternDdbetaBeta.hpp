
// Copyright (c) 2017-2023 King Abdullah University of Science and Technology,
// All rights reserved.
// ExaGeoStat is a software package, provided by King Abdullah University of Science and Technology (KAUST).

/**
 * @file UnivariateMaternDdbetaBeta.hpp
 * @brief Defines the UnivariateMaternDdbetaBeta class, a Univariate Matern Ddbeta Beta kernel.
 * @version 1.0.0
 * @author Sameh Abdulah
 * @date 2023-04-14
**/

#ifndef EXAGEOSTATCPP_UNIVARIATEMATERNDDBETABETA_HPP
#define EXAGEOSTATCPP_UNIVARIATEMATERNDDBETABETA_HPP

#include <kernels/Kernel.hpp>

namespace exageostat {
    namespace kernels {

        /**
         * @class UnivariateMaternDdbetaBeta
         * @brief A class representing a Bivariate Matern Flexible kernel.
         * @details This class represents a Bivariate Matern Flexible, which is a subclass of the Kernel class. 
         * It provides a method for generating a covariance matrix using a set of input locations and kernel parameters.
         * 
         */
        class UnivariateMaternDdbetaBeta : public Kernel {

        public:

            /**
             * @brief Constructs a new UnivariateMaternDdbetaBeta object.
             * @details Initializes a new UnivariateMaternDdbetaBeta object with default values.
             */
            UnivariateMaternDdbetaBeta();

            /**
             * @brief Virtual destructor to allow calls to the correct concrete destructor.
             * 
             */
            ~UnivariateMaternDdbetaBeta() = default;

            /**
             * @brief Generates a covariance matrix using a set of locations and kernel parameters.
             * @copydoc Kernel::GenerateCovarianceMatrix()
             */
            void GenerateCovarianceMatrix(double *apMatrixA, int &aRowsNumber, int &aColumnsNumber, int &aRowOffset,
                                          int &aColumnOffset, dataunits::Locations *apLocation1,
                                          dataunits::Locations *apLocation2, dataunits::Locations *apLocation3,
                                          double *aLocalTheta, int &aDistanceMetric) override ;

            /**
             * @brief Creates a new UnivariateMaternDdbetaBeta object.
             * @details This method creates a new UnivariateMaternDdbetaBeta object and returns a pointer to it.
             * @return A pointer to the new UnivariateMaternDdbetaBeta object.
             * 
             */
            static Kernel *Create();

        private:
            //// Used plugin name for static registration
            static bool plugin_name;
        };
    }//namespace Kernels
}//namespace exageostat

#endif //EXAGEOSTATCPP_UNIVARIATEMATERNDDBETABETA_HPP
