
// Copyright (c) 2017-2024 King Abdullah University of Science and Technology,
// All rights reserved.
// ExaGeoStat is a software package, provided by King Abdullah University of Science and Technology (KAUST).

/**
 * @file UnivariateMaternDdsigmaSquareBeta.hpp
 * @brief Defines the UnivariateMaternDdsigmaSquareBeta class, a Univariate Matern Ddsigma Square Beta kernel.
 * @version 1.1.0
 * @author Mahmoud ElKarargy
 * @author Sameh Abdulah
 * @author Suhas Shankar
 * @author Mary Lai Salvana
 * @date 2023-04-14
**/

#ifndef EXAGEOSTATCPP_UNIVARIATEMATERNDDSIGMASQUAREBETA_HPP
#define EXAGEOSTATCPP_UNIVARIATEMATERNDDSIGMASQUAREBETA_HPP

#include <kernels/Kernel.hpp>

namespace exageostat::kernels {

    /**
     * @class UnivariateMaternDdsigmaSquareBeta
     * @brief A class represents a Univariate Matern Ddsigma Square Beta kernel.
     * @details This class represents a Univariate Matern Ddsigma Square Beta, which is a subclass of the Kernel class.
     * It provides a method for generating a covariance matrix using a set of input locations and kernel parameters.
     *
     */
    template<typename T>
    class UnivariateMaternDdsigmaSquareBeta : public Kernel<T> {

    public:

        /**
         * @brief Constructs a new UnivariateMaternDdsigmaSquareBeta object.
         * @details Initializes a new UnivariateMaternDdsigmaSquareBeta object with default values.
         *
         */
        UnivariateMaternDdsigmaSquareBeta();

        /**
         * @brief Virtual destructor to allow calls to the correct concrete destructor.
         *
         */
        ~UnivariateMaternDdsigmaSquareBeta() override = default;

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
         * @brief Creates a new UnivariateMaternDdsigmaSquareBeta object.
         * @details This method creates a new UnivariateMaternDdsigmaSquareBeta object and returns a pointer to it.
         * @return A pointer to the new UnivariateMaternDdsigmaSquareBeta object.
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
    EXAGEOSTAT_INSTANTIATE_CLASS(UnivariateMaternDdsigmaSquareBeta)

}//namespace exageostat

#endif //EXAGEOSTATCPP_UNIVARIATEMATERNDDSIGMASQUAREBETA_HPP
