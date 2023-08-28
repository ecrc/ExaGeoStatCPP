
// Copyright (c) 2017-2023 King Abdullah University of Science and Technology,
// All rights reserved.
// ExaGeoStat is a software package, provided by King Abdullah University of Science and Technology (KAUST).

/**
 * @file UnivariateExpNonGaussian.hpp
 * @brief Defines the UnivariateExpNonGaussian class, a Univariate Exp Non Gaussian kernel.
 * @version 1.0.0
 * @author Suhas Shankar
 * @author Mary Lai Salvana
 * @author Sameh Abdulah
 * @author Mahmoud ElKarargy
 * @date 2023-04-14
**/

#ifndef EXAGEOSTATCPP_UNIVARIATEEXPNONGAUSSIAN_HPP
#define EXAGEOSTATCPP_UNIVARIATEEXPNONGAUSSIAN_HPP

#include <kernels/Kernel.hpp>

namespace exageostat {
    namespace kernels {

        /**
         * @class UnivariateExpNonGaussian
         * @brief A class representing a Univariate Exp Non Gaussian kernel.
         * @details This class represents a Univariate Exp Non Gaussian, which is a subclass of the Kernel class.
         * It provides a method for generating a covariance matrix using a set of input locations and kernel parameters.
         * 
         */
        template<typename T>
        class UnivariateExpNonGaussian : public Kernel<T> {

        public:

            /**
             * @brief Constructs a new UnivariateExpNonGaussian object.
             * @details Initializes a new UnivariateExpNonGaussian object with default values.
             */
            UnivariateExpNonGaussian();

            /**
             * @brief Virtual destructor to allow calls to the correct concrete destructor.
             * 
             */
            ~UnivariateExpNonGaussian() override = default;

            /**
             * @brief Generates a covariance matrix using a set of locations and kernel parameters.
             * @copydoc Kernel::GenerateCovarianceMatrix()
             */
            void GenerateCovarianceMatrix(T *apMatrixA, const int &aRowsNumber, const int &aColumnsNumber,
                                          const int &aRowOffset, const int &aColumnOffset,
                                          dataunits::Locations<T> &aLocation1, dataunits::Locations<T> &aLocation2,
                                          dataunits::Locations<T> &aLocation3, T *apLocalTheta,
                                          const int &aDistanceMetric) override;

            /**
             * @brief Creates a new UnivariateExpNonGaussian object.
             * @details This method creates a new UnivariateExpNonGaussian object and returns a pointer to it.
             * @return A pointer to the new UnivariateExpNonGaussian object.
             * 
             */
            static Kernel<T> *Create();

        private:
            //// Used plugin name for static registration
            static bool plugin_name;
        };

        /**
         * @brief Instantiates the Data Generator class for float and double types.
         * @tparam T Data Type: float or double
         *
         */
        EXAGEOSTAT_INSTANTIATE_CLASS(UnivariateExpNonGaussian)
    }//namespace Kernels
}//namespace exageostat

#endif //EXAGEOSTATCPP_UNIVARIATEEXPNONGAUSSIAN_HPP
