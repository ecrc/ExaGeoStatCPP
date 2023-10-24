
// Copyright (c) 2017-2023 King Abdullah University of Science and Technology,
// All rights reserved.
// ExaGeoStat is a software package, provided by King Abdullah University of Science and Technology (KAUST).

/**
 * @file UnivariateSpacetimeMaternStationary.hpp
 * @brief Defines the UnivariateSpacetimeMaternStationary class, a Univariate Spacetime Matern Stationary kernel.
 * @version 1.0.0
 * @author Mahmoud ElKarargy
 * @author Sameh Abdulah
 * @author Suhas Shankar
 * @author Mary Lai Salvana
 * @date 2023-04-14
**/

#ifndef EXAGEOSTATCPP_UNIVARIATESPACETIMEMATERNSTATIONARY_HPP
#define EXAGEOSTATCPP_UNIVARIATESPACETIMEMATERNSTATIONARY_HPP

#include <kernels/Kernel.hpp>

namespace exageostat::kernels {

    /**
     * @class UnivariateSpacetimeMaternStationary
     * @brief A class representing a Univariate Spacetime Matern Stationary kernel.
     * @details This class represents a Univariate Spacetime Matern Stationary, which is a subclass of the Kernel class.
     * It provides a method for generating a covariance matrix using a set of input locations and kernel parameters.
     *
     */
    template<typename T>
    class UnivariateSpacetimeMaternStationary : public Kernel<T> {

    public:

        /**
         * @brief Constructs a new UnivariateSpacetimeMaternStationary object.
         * @details Initializes a new UnivariateSpacetimeMaternStationary object with default values.
         */
        UnivariateSpacetimeMaternStationary();

        /**
         * @brief Virtual destructor to allow calls to the correct concrete destructor.
         *
         */
        ~UnivariateSpacetimeMaternStationary() override = default;

        /**
         * @brief Generates a covariance matrix using a set of locations and kernel parameters.
         * @copydoc Kernel::GenerateCovarianceMatrix()
         */
        void
        GenerateCovarianceMatrix(T *apMatrixA, const int &aRowsNumber, const int &aColumnsNumber, const int &aRowOffset,
                                 const int &aColumnOffset, dataunits::Locations<T> &aLocation1,
                                 dataunits::Locations<T> &aLocation2, dataunits::Locations<T> &aLocation3,
                                 T *apLocalTheta, const int &aDistanceMetric) override;

        /**
         * @brief Creates a new UnivariateSpacetimeMaternStationary object.
         * @details This method creates a new UnivariateSpacetimeMaternStationary object and returns a pointer to it.
         * @return A pointer to the new UnivariateSpacetimeMaternStationary object.
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
    EXAGEOSTAT_INSTANTIATE_CLASS(UnivariateSpacetimeMaternStationary)
}//namespace exageostat

#endif //EXAGEOSTATCPP_UNIVARIATESPACETIMEMATERNSTATIONARY_HPP
