
// Copyright (c) 2017-2024 King Abdullah University of Science and Technology,
// All rights reserved.
// ExaGeoStat is a software package, provided by King Abdullah University of Science and Technology (KAUST).

/**
 * @file UnivariatePowExpStationary.hpp
 * @brief Defines the UnivariatePowExpStationary class, a univariate stationary PowExp kernel.
 * @version 1.1.0
 * @author Mahmoud ElKarargy
 * @author Sameh Abdulah
 * @author Suhas Shankar
 * @author Mary Lai Salvana
 * @date 2024-11-22
 *
 * This file provides the declaration of the UnivariatePowExpStationary class, which is a subclass of the Kernel class
 * and represents a univariate stationary PowExp kernel. It provides a method for generating a covariance matrix
 * using a set of input locations and kernel parameters.
 *
**/

#ifndef EXAGEOSTATCPP_UNIVARIATEMATERNSTATIONARY_HPP
#define EXAGEOSTATCPP_UNIVARIATEMATERNSTATIONARY_HPP

#include <kernels/Kernel.hpp>

namespace exageostat::kernels {

    /**
     * @class UnivariatePowExpStationary
     * @brief A class represents a Univariate PowExp Stationary kernel.
     * @details This class represents a Univariate PowExp Stationary, which is a subclass of the Kernel class.
     * It provides a method for generating a covariance matrix using a set of input locations and kernel parameters.
     *
     */
    template<typename T>
    class UnivariatePowExpStationary : public Kernel<T> {

    public:

        /**
         * @brief Constructs a new UnivariatePowExpStationary object.
         * @details Initializes a new UnivariatePowExpStationary object with default values.
         *
         */
        UnivariatePowExpStationary();

        /**
         * @brief Virtual destructor to allow calls to the correct concrete destructor.
         *
         */
        ~UnivariatePowExpStationary() override = default;

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
         * @brief Creates a new UnivariatePowExpStationary object.
         * @details This method creates a new UnivariatePowExpStationary object and returns a pointer to it.
         * @return A pointer to the new UnivariatePowExpStationary object.
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
    EXAGEOSTAT_INSTANTIATE_CLASS(UnivariatePowExpStationary)

}//namespace exageostat

#endif //EXAGEOSTATCPP_UNIVARIATEMATERNSTATIONARY_HPP
