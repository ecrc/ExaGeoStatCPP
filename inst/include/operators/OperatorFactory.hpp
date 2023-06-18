// Copyright (c) 2017-2023 King Abdullah University of Science and Technology,
// All rights reserved.
// ExaGeoStat is a software package, provided by King Abdullah University of Science and Technology (KAUST).

/**
 * @file OperatorFactory.hpp
 * @brief Header file for the OperatorFactory class, which creates operators based on the input operator type.
 * @version 1.0.0
 * @author Sameh Abdulah
 * @date 2023-04-30
 *
 * This file contains the declaration of the OperatorFactory class, which is responsible for creating operators based on the input operator type.
 * The OperatorFactory class is templated on the data type of the operator, and provides a static CreateOperator() method for creating new operator instances.
 *
**/

#ifndef EXAGEOSTATCPP_OPERATORFACTORY_HPP
#define EXAGEOSTATCPP_OPERATORFACTORY_HPP

#include <common/Definitions.hpp>
#include <operators/OperatorMethods.hpp>

namespace exageostat {
    namespace operators {

        /**
         * @class OperatorFactory
         * @brief A class that creates operators based on the input operator type.
         * @tparam T Data Type: float or double
         */
        template<typename T>
        class OperatorFactory {

        public:

            /**
             * @brief Creates an operator based on the input operator type.
             * @param[in] aOperator The operator type to create.
             * @return A unique pointer to the created operator.
             *
             */
            static std::unique_ptr<OperatorMethods<T>>
            CreateOperator(common::Operators aOperator);
        };

        /**
        * @brief Instantiates the operator methods class for float and double types.
        * @tparam T Data Type: float or double
        *
        */
        EXAGEOSTAT_INSTANTIATE_CLASS(OperatorFactory)

    }//namespace operators
}//namespace exageostat


#endif //EXAGEOSTATCPP_OPERATORFACTORY_HPP
