
// Copyright (c) 2017-2024 King Abdullah University of Science and Technology,
// All rights reserved.
// ExaGeoStat is a software package, provided by King Abdullah University of Science and Technology (KAUST).

/**
 * @file TrendModel.hpp
 * @brief Defines the TrendModel class – a time–trend kernel that can be used as a mean/covariate model.
 * @version 2.0.0
 * @author Mahmoud ElKarargy
 * @author Sameh Abdulah
 * @date 2024-11-11
**/

#ifndef EXAGEOSTATCPP_TRENDMODEL_HPP
#define EXAGEOSTATCPP_TRENDMODEL_HPP

#include <kernels/Kernel.hpp>

namespace exageostat::kernels {

    /**
     * @class TrendModel
     * @brief A simple temporal trend model that can be plugged into the unified kernel interface.
     *
     * The implementation replicates the time-series trend model that is originally implemented in the
     * C version of ExaGeoStat.  The class inherits from the abstract `Kernel<T>` and therefore must
     * provide an implementation for `GenerateCovarianceMatrix` in which the design matrix (or the
     * corresponding covariance when used jointly with other kernels) is produced.
     */
    template<typename T>
    class TrendModel : public Kernel<T> {
    public:
        /**
         * @brief Default constructor – sets the internal bookkeeping variables.
         */
        TrendModel();

        /**
         * @brief Virtual destructor to allow calls to the correct concrete destructor.
         */
        ~TrendModel() override = default;

        /**
         * @brief Generates a covariance (or design) matrix block for the trend model.
         * @copydoc Kernel::GenerateCovarianceMatrix()
         */
        void GenerateCovarianceMatrix(T *apMatrixA, const int &aRowsNumber, const int &aColumnsNumber,
                                      const int &aRowOffset, const int &aColumnOffset, dataunits::Locations<T> &aLocation1,
                                      dataunits::Locations<T> &aLocation2, dataunits::Locations<T> &aLocation3,
                                      T *apLocalTheta, const int &aDistanceMetric) override;

        /**
         * @brief Factory function required by the plugin registry mechanism.
         * @return A newly allocated `TrendModel` object.
         */
        static Kernel<T> *Create();

    private:
        //// Used plugin name for static registration
        static bool plugin_name;
    };

    /**
     * @brief Instantiate the TrendModel class for the supported precision types.
     */
    EXAGEOSTAT_INSTANTIATE_CLASS(TrendModel)

}// namespace exageostat::kernels

#endif // EXAGEOSTATCPP_TRENDMODEL_HPP
