
// Copyright (c) 2017-2024 King Abdullah University of Science and Technology,
// All rights reserved.
// ExaGeoStat is a software package, provided by King Abdullah University of Science and Technology (KAUST).

/**
 * @file StarPuHelpersFactory.hpp
 * @brief Factory for StarPu helpers.
 * @version 1.1.0
 * @author Mahmoud ElKarargy
 * @date 2024-02-25
**/

#ifndef EXAGEOSTATCPP_STARPUHELPERSFACTORY_HPP
#define EXAGEOSTATCPP_STARPUHELPERSFACTORY_HPP

#include <runtime/starpu/helpers/StarPuHelpers.hpp>

namespace exageostat::runtime {

    /**
     * @class StarPuHelpersFactory
     * @brief A class that creates StarPu helpers based on the input computation type.
     *
     */
    class StarPuHelpersFactory {

    public:

        /**
         * @brief Creates a StarPu helper.
         * @param[in] aComputation The computation type to create the solver for.
         * @return Unique pointer to the created StarPu helper.
         *
         */
        static std::unique_ptr<StarPuHelpers> CreateStarPuHelper(const common::Computation &aComputation);

    };

}//namespace exageostat

#endif //EXAGEOSTATCPP_STARPUHELPERSFACTORY_HPP
