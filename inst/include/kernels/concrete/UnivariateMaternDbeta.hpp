
// Copyright (c) 2017-2023 King Abdullah University of Science and Technology,
// All rights reserved.
// ExaGeoStat is a software package, provided by King Abdullah University of Science and Technology (KAUST).

/**
 * @file UnivariateMaternDbeta.hpp
 * @brief Defines the UnivariateMaternDbeta class, a Univariate Matern Dbeta kernel.
 * @version 1.0.0
 * @author Sameh Abdulah
 * @date 2023-04-14
**/

#ifndef EXAGEOSTATCPP_UNIVARIATEMATERNDBETA_HPP
#define EXAGEOSTATCPP_UNIVARIATEMATERNDBETA_HPP

#include <kernels/Kernel.hpp>

namespace exageostat {
    namespace kernels {

        /**
         * @class UnivariateMaternDbeta
         * @brief A class representing a Bivariate Matern Flexible kernel.
         * @details This class represents a Bivariate Matern Flexible, which is a subclass of the Kernel class. 
         * It provides a method for generating a covariance matrix using a set of input locations and kernel parameters.
         * 
         */
        template<typename T>
        class UnivariateMaternDbeta : public Kernel<T> {

        public:

            /**
             * @brief Constructs a new UnivariateMaternDbeta object.
             * @details Initializes a new UnivariateMaternDbeta object with default values.
             */
            UnivariateMaternDbeta();

            /**
             * @brief Virtual destructor to allow calls to the correct concrete destructor.
             * 
             */
            ~UnivariateMaternDbeta() = default;

            /**
             * @brief Generates a covariance matrix using a set of locations and kernel parameters.
             * @copydoc Kernel::GenerateCovarianceMatrix()
             */
            void GenerateCovarianceMatrix(T *apMatrixA, int &aRowsNumber, int &aColumnsNumber, int &aRowOffset,
                                          int &aColumnOffset, dataunits::Locations<T> *apLocation1,
                                          dataunits::Locations<T> *apLocation2, dataunits::Locations<T> *apLocation3,
                                          T *aLocalTheta, int &aDistanceMetric) override ;

            /**
             * @brief Creates a new UnivariateMaternDbeta object.
             * @details This method creates a new UnivariateMaternDbeta object and returns a pointer to it.
             * @return A pointer to the new UnivariateMaternDbeta object.
             * 
             */
            static Kernel<T> *Create();

        private:
            //// Used plugin name for static registration
            static bool plugin_name;
        };
    }//namespace Kernels
}//namespace exageostat

#endif //EXAGEOSTATCPP_UNIVARIATEMATERNDBETA_HPP
