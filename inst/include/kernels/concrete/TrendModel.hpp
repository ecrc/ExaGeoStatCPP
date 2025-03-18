
// Copyright (c) 2017-2024 King Abdullah University of Science and Technology,
// All rights reserved.
// ExaGeoStat is a software package, provided by King Abdullah University of Science and Technology (KAUST).

/**
 * @file TrendModel.hpp
 * @brief Defines the TrendModel kernel class.
 * @version 2.0.0
 * @author Mahmoud ElKarargy
 * @author Sameh Abdulah
 * @date 2024-11-11
 *
 * This file provides the declaration of the TrendModel class, which is a subclass of the Kernel class.
 * It provides a method for generating a covariance matrix
 * using a set of input locations and kernel parameters.
 *
**/

#ifndef EXAGEOSTATCPP_TRENDMODEL_HPP
#define EXAGEOSTATCPP_TRENDMODEL_HPP

#include <kernels/Kernel.hpp>

namespace exageostat::kernels {

    /**
     * @class TrendModel
     * @brief A class represents a Trend Model kernel.
     * @details This class represents a Trend Model Stationary, which is a subclass of the Kernel class.
     * It provides a method for generating a covariance matrix using a set of input locations and kernel parameters.
     *
     */
    template<typename T>
    class TrendModel : public Kernel<T> {

    public:

        /**
         * @brief Constructs a new TrendModel object.
         * @details Initializes a new TrendModel object with default values.
         *
         */
        TrendModel();

        /**
         * @brief Virtual destructor to allow calls to the correct concrete destructor.
         *
         */
        ~TrendModel() override = default;

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
         * @brief Creates a new TrendModel object.
         * @details This method creates a new TrendModel object and returns a pointer to it.
         * @return A pointer to the new TrendModel object.
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
    EXAGEOSTAT_INSTANTIATE_CLASS(TrendModel)

}//namespace exageostat

#endif //EXAGEOSTATCPP_TRENDMODEL_HPP
