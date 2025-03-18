
// Copyright (c) 2017-2024 King Abdullah University of Science and Technology,
// All rights reserved.
// ExaGeoStat is a software package, provided by King Abdullah University of Science and Technology (KAUST).


/**
 * @file Kernels.hpp
 * @brief Header file for the Kernels class, which contains the main kernel functions.
 * @version 1.1.0
 * @author Mahmoud ElKarargy
 * @author Sameh Abdulah
 * @author Suhas Shankar
 * @author Mary Lai Salvana
 * @date 2023-03-05
**/

#ifndef EXAGEOSTAT_CPP_KERNELS_HPP
#define EXAGEOSTAT_CPP_KERNELS_HPP

#include<cmath>

#if DEFAULT_RUNTIME
#include <starpu.h>
#endif

extern "C" {
#include <gsl/gsl_sf_bessel.h>
#include <gsl/gsl_sf_psi.h>
}

#include <common/PluginRegistry.hpp>
#include <data-units/Locations.hpp>
#include <helpers/DistanceCalculationHelpers.hpp>
#include <helpers/BasselFunction.hpp>

/**
 * @def EARTH_RADIUS
 * @brief The radius of the Earth in kilometers.
 * @details This macro defines the radius of the Earth in kilometers, which is used by the Great Circle Distance (GCD) function.
 *
 */
#define EARTH_RADIUS 6371.0

namespace exageostat::kernels {

    struct KernelsConfigurations {
        /**
         * @brief Returns the static map containing kernel parameter numbers.
         * @return Reference to the static map.
         */
        static std::unordered_map<std::string, int> &GetParametersNumberKernelMap() {
            /**
             * @brief Static map containing kernel parameter numbers.
             * @details The map is initialized only once and retains its value across multiple function invocations.
             */
            static std::unordered_map<std::string, int> mKernelParametersNumbers;
            return mKernelParametersNumbers;
        }
    };

    /**
     * @class Kernels
     * @brief A base class for kernel functions.
     * @details This class provides a base class for kernel functions and contains several utility functions for computing distance metrics and Bessel functions.
     *
     */
    template<typename T>
    class Kernel {
    public:

        /**
         * Default virtual destructor to be overridden by the the suitable concrete kernel destructor.
         * 
         */
        virtual ~Kernel() = default;

        /**
         * @brief Generates a covariance matrix using a set of locations and kernel parameters.
         * @param[out] apMatrixA The output covariance matrix.
         * @param[in] aRowsNumber The number of rows in the output matrix.
         * @param[in] aColumnsNumber The number of columns in the output matrix.
         * @param[in] aRowOffset The row offset for the input locations.
         * @param[in] aColumnOffset The column offset for the input locations.
         * @param[in] apLocation1 The set of input locations 1.
         * @param[in] apLocation2 The set of input locations 2.
         * @param[in] apLocation3 The set of input locations 3.
         * @param[in] aLocalTheta An array of kernel parameters.
         * @param [in] aDistanceMetric Distance metric to be used (1 = Euclidean, 2 = Manhattan, 3 = Minkowski).
         * @return void
         *
         */
        virtual void
        GenerateCovarianceMatrix(T *apMatrixA, const int &aRowsNumber, const int &aColumnsNumber, const int &aRowOffset,
                                 const int &aColumnOffset, dataunits::Locations<T> &aLocation1,
                                 dataunits::Locations<T> &aLocation2, dataunits::Locations<T> &aLocation3,
                                 T *apLocalTheta, const int &aDistanceMetric) = 0;

        /**
         * @brief Returns the value of the parameter P used by the kernel function.
         * @return The value of P (Variables Number).
         *
         */
        [[nodiscard]] int GetVariablesNumber() const;

        /**
         * @brief Sets the value of the parameter P used by the kernel function.
         * @param[in] aTimeSlot Value to set `mP` with.
         * @return void
         *
         */
        void SetPValue(int aTimeSlot);

        /**
         * @brief Returns the number of the parameters used by the kernel function.
         * @return The value of ParametersNumber.
         *
         */
        [[nodiscard]] int GetParametersNumbers() const;

    protected:
        //// Used P.
        int mP = 1;
        //// Used Variable number which is P multiplied by timeslot
        int mVariablesNumber = 1;
        //// Used number of parameters.
        int mParametersNumber = 3;
    };

    /**
     * @brief Instantiates the Data Generator class for float and double types.
     * @tparam T Data Type: float or double
     *
     */
    EXAGEOSTAT_INSTANTIATE_CLASS(Kernel)

}//namespace exageostat

#endif //EXAGEOSTAT_CPP_KERNELS_HPP