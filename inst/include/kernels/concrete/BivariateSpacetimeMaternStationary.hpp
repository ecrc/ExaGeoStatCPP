
// Copyright (c) 2017-2024 King Abdullah University of Science and Technology,
// All rights reserved.
// ExaGeoStat is a software package, provided by King Abdullah University of Science and Technology (KAUST).

/**
 * @file BivariateSpacetimeMaternStationary.hpp
 * @brief Defines the BivariateSpacetimeMaternStationary class, a Bivariate Spacetime Matern Stationary kernel.
 * @version 1.1.0
 * @author Mahmoud ElKarargy
 * @author Sameh Abdulah
 * @author Suhas Shankar
 * @author Mary Lai Salvana
 * @date 2023-04-14
**/

#ifndef EXAGEOSTATCPP_BIVARIATESPACETIMEMATERNSTATIONARY_HPP
#define EXAGEOSTATCPP_BIVARIATESPACETIMEMATERNSTATIONARY_HPP

#include <kernels/Kernel.hpp>

namespace exageostat::kernels {

    /**
     * @class BivariateSpacetimeMaternStationary
     * @brief A class represents a Bivariate Spacetime Matern Stationary kernel.
     * @details This class represents a Bivariate Spacetime Matern Stationary, which is a subclass of the Kernel class.
     * It provides a method for generating a covariance matrix using a set of input locations and kernel parameters.
     *
     */
    template<typename T>
    class BivariateSpacetimeMaternStationary : public Kernel<T> {

    public:

        /**
         * @brief Constructs a new BivariateSpacetimeMaternStationary object.
         * @details Initializes a new BivariateSpacetimeMaternStationary object with default values.
         *
         */
        BivariateSpacetimeMaternStationary();

        /**
         * @brief Virtual destructor to allow calls to the correct concrete destructor.
         *
         */
        ~BivariateSpacetimeMaternStationary() override = default;

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
         * @brief Creates a new BivariateSpacetimeMaternStationary object.
         * @details This method creates a new BivariateSpacetimeMaternStationary object and returns a pointer to it.
         * @return A pointer to the new BivariateSpacetimeMaternStationary object.
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
    EXAGEOSTAT_INSTANTIATE_CLASS(BivariateSpacetimeMaternStationary)
}//namespace exageostat

#endif //EXAGEOSTATCPP_BIVARIATESPACETIMEMATERNSTATIONARY_HPP
