
// Copyright (c) 2017-2023 King Abdullah University of Science and Technology,
// All rights reserved.
// ExaGeoStat is a software package, provided by King Abdullah University of Science and Technology (KAUST).

/**
 * @file UnivariateMaternDdbetaNu.hpp
 * @brief Defines the UnivariateMaternDdbetaNu class, a Univariate Matern Ddbeta Nu kernel.
 * @version 1.0.0
 * @author Sameh Abdulah
 * @date 2023-04-14
**/

#ifndef EXAGEOSTATCPP_UNIVARIATEMATERNDDBETANU_HPP
#define EXAGEOSTATCPP_UNIVARIATEMATERNDDBETANU_HPP

#include <kernels/Kernel.hpp>

namespace exageostat {
    namespace kernels {

        /**
         * @class UnivariateMaternDdbetaNu
         * @brief A class representing a Univariate Matern Ddbeta Nu kernel.
         * @details This class represents a Univariate Matern Ddbeta Nu, which is a subclass of the Kernel class.
         * It provides a method for generating a covariance matrix using a set of input locations and kernel parameters.
         * 
         */
        template<typename T>
        class UnivariateMaternDdbetaNu : public Kernel<T> {

        public:

            /**
             * @brief Constructs a new UnivariateMaternDdbetaNu object.
             * @details Initializes a new UnivariateMaternDdbetaNu object with default values.
             */
            UnivariateMaternDdbetaNu();

            /**
             * @brief Virtual destructor to allow calls to the correct concrete destructor.
             * 
             */
            ~UnivariateMaternDdbetaNu() override = default;

            /**
             * @brief Generates a covariance matrix using a set of locations and kernel parameters.
             * @copydoc Kernel::GenerateCovarianceMatrix()
             */
            void GenerateCovarianceMatrix(T *apMatrixA, int &aRowsNumber, int &aColumnsNumber, int &aRowOffset,
                                          int &aColumnOffset, dataunits::Locations<T> *apLocation1,
                                          dataunits::Locations<T> *apLocation2, dataunits::Locations<T> *apLocation3,
                                          T *aLocalTheta, int &aDistanceMetric) override;

            /**
             * @brief Creates a new UnivariateMaternDdbetaNu object.
             * @details This method creates a new UnivariateMaternDdbetaNu object and returns a pointer to it.
             * @return A pointer to the new UnivariateMaternDdbetaNu object.
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
        EXAGEOSTAT_INSTANTIATE_CLASS(UnivariateMaternDdbetaNu)
    }//namespace Kernels
}//namespace exageostat

#endif //EXAGEOSTATCPP_UNIVARIATEMATERNDDBETANU_HPP
