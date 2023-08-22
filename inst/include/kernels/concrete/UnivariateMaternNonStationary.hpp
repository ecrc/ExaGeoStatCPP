
// Copyright (c) 2017-2023 King Abdullah University of Science and Technology,
// All rights reserved.
// ExaGeoStat is a software package, provided by King Abdullah University of Science and Technology (KAUST).

/**
 * @file UnivariateMaternNonStationary.hpp
 * @brief Defines the UnivariateMaternNonStationary class, a Univariate Matern Non Stationary kernel.
 * @version 1.0.0
 * @author Suhas Shankar
 * @author Mary Lai Salvana
 * @author Sameh Abdulah
 * @author Mahmoud ElKarargy
 * @date 2023-04-13
**/

#ifndef EXAGEOSTATCPP_UNIVARIATEMATERNNONSTATIONARY_HPP
#define EXAGEOSTATCPP_UNIVARIATEMATERNNONSTATIONARY_HPP

#include <kernels/Kernel.hpp>

namespace exageostat {
    namespace kernels {

        /**
         * @class UnivariateMaternNonStationary
         * @brief A class representing a Univariate Matern Non Stationary kernel.
         * @details This class represents a Univariate Matern Non Stationary, which is a subclass of the Kernel class.
         * It provides a method for generating a covariance matrix using a set of input locations and kernel parameters.
         *
         */
        template<typename T>
        class UnivariateMaternNonStationary : public Kernel<T> {

        public:

            /**
             * @brief Constructs a new UnivariateMaternNonStationary object.
             * @details Initializes a new UnivariateMaternNonStationary object with default values.
             */
            UnivariateMaternNonStationary();

            /**
             * @brief Virtual destructor to allow calls to the correct concrete destructor.
             *
             */
            ~UnivariateMaternNonStationary() override = default;

            /**
             * @brief Generates a covariance matrix using a set of locations and kernel parameters.
             * @copydoc Kernel::GenerateCovarianceMatrix()
             */
            void GenerateCovarianceMatrix(T *apMatrixA, int &aRowsNumber, int &aColumnsNumber, int &aRowOffset,
                                          int &aColumnOffset, dataunits::Locations<T> *apLocation1,
                                          dataunits::Locations<T> *apLocation2, dataunits::Locations<T> *apLocation3,
                                          T *aLocalTheta, int &aDistanceMetric) override;

            /**
             * @brief Creates a new UnivariateMaternNonStationary object.
             * @details This method creates a new UnivariateMaternNonStationary object and returns a pointer to it.
             * @return A pointer to the new UnivariateMaternNonStationary object.
             *
             */
            static Kernel<T> *Create();

        private:
            //// Used plugin name for static registration
            static bool plugin_name;
        };
    }//namespace Kernels
}//namespace exageostat

#endif //EXAGEOSTATCPP_UNIVARIATEMATERNNONSTATIONARY_HPP
