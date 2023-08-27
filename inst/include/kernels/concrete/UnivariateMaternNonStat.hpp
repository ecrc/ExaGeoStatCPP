
// Copyright (c) 2017-2023 King Abdullah University of Science and Technology,
// All rights reserved.
// ExaGeoStat is a software package, provided by King Abdullah University of Science and Technology (KAUST).

/**
 * @file UnivariateMaternNonStat.hpp
 * @brief Defines the UnivariateMaternNonStat class, a Univariate Matern Non Stat kernel.
 * @version 1.0.0
 * @author Suhas Shankar
 * @author Mary Lai Salvana
 * @author Sameh Abdulah
 * @author Mahmoud ElKarargy
 * @date 2023-04-14
**/

#ifndef EXAGEOSTATCPP_UNIVARIATEMATERNNONSTAT_HPP
#define EXAGEOSTATCPP_UNIVARIATEMATERNNONSTAT_HPP

#include <kernels/Kernel.hpp>

namespace exageostat {
    namespace kernels {

        /**
         * @class UnivariateMaternNonStat
         * @brief A class representing a Univariate Matern Non Stat kernel.
         * @details This class represents a Univariate Matern Non Stat,which is a subclass of the Kernel class.
         * It provides a method for generating a covariance matrix using a set of input locations and kernel parameters.
         *
         */
        template<typename T>
        class UnivariateMaternNonStat : public Kernel<T> {

        public:

            /**
             * @brief Constructs a new UnivariateMaternNonStat object.
             * @details Initializes a new UnivariateMaternNonStat object with default values.
             */
            UnivariateMaternNonStat();

            /**
             * @brief Virtual destructor to allow calls to the correct concrete destructor.
             *
             */
            ~UnivariateMaternNonStat() override = default;

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
             * @brief Creates a new UnivariateMaternNonStat object.
             * @details This method creates a new UnivariateMaternNonStat object and returns a pointer to it.
             * @return A pointer to the new UnivariateMaternNonStat object.
             *
             */
            static Kernel<T> *Create();

            /**
             * Function for smoothness parameter
             * @param[in] aX x co-ordinate
             * @param[in] aY y co-ordinate
             * @param[in] aG parameter for function
             * @param[in] aH parameter for function
             * @param[in] aTi parameter for function
             * @return The function ge^(h(x+y)) + ti
             *
             */
            static double Neu(double aX, double aY, double aG, double aH, double aTi);

            /**
             * Function for partial sill
             * @param[in] aX x co-ordinate
             * @param[in] aY y co-ordinate
             * @param[in] aD parameter for function
             * @param[in] aE parameter for function
             * @param[in] aF parameter for function
             * @return The function de^(e(x+y)) + f
             *
             */
            static double Sigma(double aX, double aY, double aD, double aE, double aF);

            /**
             * Function for spatial range
              * @param[in] aX x co-ordinate
             * @param[in] aY y co-ordinate
             * @param[in] aA parameter for function
             * @param[in] aB parameter for function
             * @return The function ae^(sin bx + sin by)
             *
             */
            static double Lambda(double aX, double aY, double aA, double aB);

            /**
             * Returns the Mahalanobis distance between two points to account for anisotropy
             * @param[in] aX1 x co-ordinate of first point
             * @param[in] aY1 y co-ordinate of first point
             * @param[in] aX2 x co-ordinate of second point
             * @param[in] aY2 y co-ordinate of second point
             * @param[in] aA11 First element of the positive definite matrix that defines the Mahalanobis Distance
             * @param[in] aA12 Second element of the positive definite matrix that defines the Mahalanobis Distance
             * @param[in] aA21 Third element of the positive definite matrix that defines the Mahalanobis Distance
             * @param[in] aA22 Fourth element of the positive definite matrix that defines the Mahalanobis Distance
             * @return The Mahalanobis Distance
             *
             */
            static double CalculateMahalanobisDistanceSquared(double aX1, double aY1, double aX2,
                                                              double aY2, double aA11, double aA12,
                                                              double aA21, double aA22);

            /**
             * Utility function that evaluates the matern. Similiar to (https://www.rdocumentation.org/packages/PrevMap/versions/1.5.3/topics/matern.kernel) in R
             * @param[in] aRange Spatial Range parameter (Also known as rho)
             * @param[in] aSmoothness Smoothness parameter (Also known as neu)
             * @param[in] aDistance Distance between the two locations
             * @return Matern function evaluation
             *
             */
            static double MaternUtil(double aRange, double aSmoothness, double aDistance);

        private:
            //// Used plugin name for static registration
            static bool plugin_name;
        };

        /**
         * @brief Instantiates the Data Generator class for float and double types.
         * @tparam T Data Type: float or double
         *
         */
        EXAGEOSTAT_INSTANTIATE_CLASS(UnivariateMaternNonStat)
    }//namespace Kernels
}//namespace exageostat

#endif //EXAGEOSTATCPP_UNIVARIATEMATERNNONSTAT_HPP
