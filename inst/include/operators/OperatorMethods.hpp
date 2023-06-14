// Copyright (c) 2017-2023 King Abdullah University of Science and Technology,
// Copyright (C) 2023 by Brightskies inc,
// All rights reserved.
// ExaGeoStat is a software package, provided by King Abdullah University of Science and Technology (KAUST).

/**
 * @file OperatorMethods.hpp
 * @brief Header file for the OperatorMethods class, which defines the interface for operators.
 * @version 1.0.0
 * @author Sameh Abdulah
 * @date 2023-04-30
 *
 * This file contains the declaration of the OperatorMethods class, which defines the interface for operators.
 * The OperatorMethods class is templated on the data type of the operator, and does not provide any method implementations.
**/

#ifndef EXAGEOSTATCPP_OPERATORMETHODS_HPP
#define EXAGEOSTATCPP_OPERATORMETHODS_HPP

#include <common/Definitions.hpp>

namespace exageostat {
    namespace operators {
        /**
         * @class OperatorMethods
         * @brief A class that defines the interface for operators.
         * @tparam T The data type of the operator.
         */
        template<typename T>
        class OperatorMethods {

        public:

        };

        EXAGEOSTAT_INSTANTIATE_CLASS(OperatorMethods)
    }//namespace operators
}//namespace exageostat
#endif //EXAGEOSTATCPP_OPERATORMETHODS_HPP
