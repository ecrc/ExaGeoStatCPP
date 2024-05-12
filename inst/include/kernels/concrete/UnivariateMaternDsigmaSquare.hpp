
// Copyright (c) 2017-2024 King Abdullah University of Science and Technology,
// All rights reserved.
// ExaGeoStat is a software package, provided by King Abdullah University of Science and Technology (KAUST).

/**
 * @file UnivariateMaternDsigmaSquare.hpp
 * @brief Defines the UnivariateMaternDsigmaSquare class, a Univariate Matern Dsigma Square kernel.
 * @version 1.1.0
 * @author Mahmoud ElKarargy
 * @author Sameh Abdulah
 * @author Suhas Shankar
 * @author Mary Lai Salvana
 * @date 2023-04-14
**/

#ifndef EXAGEOSTATCPP_UNIVARIATEMATERNDSIGMASQUARE_HPP
#define EXAGEOSTATCPP_UNIVARIATEMATERNDSIGMASQUARE_HPP

#include <kernels/Kernel.hpp>

namespace exageostat::kernels {

    /**
     * @class UnivariateMaternDsigmaSquare
     * @brief A class represents a Univariate Matern Dsigma Square kernel.
     * @details This class represents a Univariate Matern Dsigma Square, which is a subclass of the Kernel class.
     * It provides a method for generating a covariance matrix using a set of input locations and kernel parameters.
     *
     */
    template<typename T>
    class UnivariateMaternDsigmaSquare : public Kernel<T> {

    public:

        /**
         * @brief Constructs a new UnivariateMaternDsigmaSquare object.
         * @details Initializes a new UnivariateMaternDsigmaSquare object with default values.
         *
         */
        UnivariateMaternDsigmaSquare();

        /**
         * @brief Virtual destructor to allow calls to the correct concrete destructor.
         *
         */
        ~UnivariateMaternDsigmaSquare() override = default;

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
         * @brief Creates a new UnivariateMaternDsigmaSquare object.
         * @details This method creates a new UnivariateMaternDsigmaSquare object and returns a pointer to it.
         * @return A pointer to the new UnivariateMaternDsigmaSquare object.
         *
         */
        static Kernel<T> *Create();

    private:
        //// Used plugin name for static registration
        static bool plugin_name;
    };
}//namespace exageostat

#endif //EXAGEOSTATCPP_UNIVARIATEMATERNDSIGMASQUARE_HPP
