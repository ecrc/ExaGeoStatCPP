
// Copyright (c) 2017-2023 King Abdullah University of Science and Technology,
// All rights reserved.
// ExaGeoStat is a software package, provided by King Abdullah University of Science and Technology (KAUST).

/**
 * @file UnivariateMaternDdsigmaSquareNu.hpp
 * @brief Defines the UnivariateMaternDdsigmaSquareNu class, a Univariate Matern Ddsigma Square Nu kernel.
 * @version 1.0.0
 * @author Suhas Shankar
 * @author Mary Lai Salvana
 * @author Sameh Abdulah
 * @author Mahmoud ElKarargy
 * @date 2023-04-14
**/

#ifndef EXAGEOSTATCPP_UNIVARIATEMATERNDDSIGMASQUARENU_HPP
#define EXAGEOSTATCPP_UNIVARIATEMATERNDDSIGMASQUARENU_HPP

#include <kernels/Kernel.hpp>

namespace exageostat {
    namespace kernels {

        /**
         * @class UnivariateMaternDdsigmaSquareNu
         * @brief A class representing a Univariate Matern Ddsigma Square Nu kernel.
         * @details This class represents a Univariate Matern Ddsigma Square Nu, which is a subclass of the Kernel class.
         * It provides a method for generating a covariance matrix using a set of input locations and kernel parameters.
         *
         */
        template<typename T>
        class UnivariateMaternDdsigmaSquareNu : public Kernel<T> {

        public:

            /**
             * @brief Constructs a new UnivariateMaternDdsigmaSquareNu object.
             * @details Initializes a new UnivariateMaternDdsigmaSquareNu object with default values.
             */
            UnivariateMaternDdsigmaSquareNu();

            /**
             * @brief Virtual destructor to allow calls to the correct concrete destructor.
             *
             */
            ~UnivariateMaternDdsigmaSquareNu() override = default;

            /**
             * @brief Generates a covariance matrix using a set of locations and kernel parameters.
             * @copydoc Kernel::GenerateCovarianceMatrix()
             */
            void GenerateCovarianceMatrix(T *apMatrixA, int &aRowsNumber, int &aColumnsNumber, int &aRowOffset,
                                          int &aColumnOffset, dataunits::Locations<T> *apLocation1,
                                          dataunits::Locations<T> *apLocation2, dataunits::Locations<T> *apLocation3,
                                          T *aLocalTheta, int &aDistanceMetric) override;

            /**
             * @brief Creates a new UnivariateMaternDdsigmaSquareNu object.
             * @details This method creates a new UnivariateMaternDdsigmaSquareNu object and returns a pointer to it.
             * @return A pointer to the new UnivariateMaternDdsigmaSquareNu object.
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
        EXAGEOSTAT_INSTANTIATE_CLASS(UnivariateMaternDdsigmaSquareNu)

    }//namespace Kernels
}//namespace exageostat

#endif //EXAGEOSTATCPP_UNIVARIATEMATERNDDSIGMASQUARENU_HPP
