
// Copyright (c) 2017-2024 King Abdullah University of Science and Technology,
// All rights reserved.
// ExaGeoStat is a software package, provided by King Abdullah University of Science and Technology (KAUST).

/**
 * @file UnivariateMaternDdsigmaSquare.hpp
 * @brief Defines the UnivariateMaternDdsigmaSquare class, a univariate stationary Matern kernel.
 * @version 1.1.0
 * @author Mahmoud ElKarargy
 * @author Sameh Abdulah
 * @author Suhas Shankar
 * @author Mary Lai Salvana
 * @date 2023-04-12
 *
 * This file provides the declaration of the UnivariateMaternDdsigmaSquare class, which is a subclass of the Kernel class
 * and represents a univariate stationary Matern kernel. It provides a method for generating a covariance matrix
 * using a set of input locations and kernel parameters.
 *
**/

#ifndef EXAGEOSTATCPP_UNIVARIATEMATERNDDSIGMASQUARE_HPP
#define EXAGEOSTATCPP_UNIVARIATEMATERNDDSIGMASQUARE_HPP

#include <kernels/Kernel.hpp>

namespace exageostat::kernels {

    /**
     * @class UnivariateMaternDdsigmaSquare
     * @brief A class represents a Univariate Matern Ddsigma Square kernel.
     * @details This class represents a Univariate Matern Ddsigma Square, which is a subclass of the Kernel class.
     * It provides a method for generating a covariance matrix using a set of input locations and kernel parameters.
     *
     */
    template<typename T>
    class UnivariateMaternDdsigmaSquare : public Kernel<T> {

    public:

        /**
         * @brief Constructs a new UnivariateMaternDdsigmaSquare object.
         * @details Initializes a new UnivariateMaternDdsigmaSquare object with default values.
         *
         */
        UnivariateMaternDdsigmaSquare();

        /**
         * @brief Virtual destructor to allow calls to the correct concrete destructor.
         *
         */
        ~UnivariateMaternDdsigmaSquare() override = default;

        /**
         * @brief Generates a covariance matrix using a set of locations and kernel parameters.
         * @copydoc Kernel::GenerateCovarianceMatrix()
         *
         */
        void
        GenerateCovarianceMatrix(T *apMatrixA, const int &aRowsNumber, const int &aColumnsNumber, const int &aRowOffset,
                                 const int &aColumnOffset, dataunits::Locations<T> &aLocation1,
                                 dataunits::Locations<T> &aLocation2, dataunits::Locations<T> &aLocation3,
                                 T *apLocalTheta, const int &aDistanceMetric) override;

        /**
         * @brief Creates a new UnivariateMaternDdsigmaSquare object.
         * @details This method creates a new UnivariateMaternDdsigmaSquare object and returns a pointer to it.
         * @return A pointer to the new UnivariateMaternDdsigmaSquare object.
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
    EXAGEOSTAT_INSTANTIATE_CLASS(UnivariateMaternDdsigmaSquare)

}//namespace exageostat

#endif //EXAGEOSTATCPP_UNIVARIATEMATERNDDSIGMASQUARE_HPP
