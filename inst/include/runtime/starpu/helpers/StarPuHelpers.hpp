
// Copyright (c) 2017-2024 King Abdullah University of Science and Technology,
// All rights reserved.
// ExaGeoStat is a software package, provided by King Abdullah University of Science and Technology (KAUST).

/**
 * @file StarPuHelpers.hpp
 * @brief An interface for StarPu helpers.
 * @version 1.1.0
 * @author Mahmoud ElKarargy
 * @date 2024-02-25
**/

#ifndef EXAGEOSTATCPP_STARPUHELPERS_HPP
#define EXAGEOSTATCPP_STARPUHELPERS_HPP

#include <complex>

#include <common/Definitions.hpp>
#include <hardware/ExaGeoStatHardware.hpp>

namespace exageostat::runtime {

    /**
     * @class StarPuHelpers
     * @brief A class that defines the interface for StarPu helpers.
     * @tparam T Data Type: float or double.
     *
     */
    class StarPuHelpers {
    public:

        /**
         * @brief Initialize the runtime option structure for either HiCMA or CHAMELEON.
         * @param[in, out] apOptions The options structure that needs to be initialized.
         * @param[in] apSequence The sequence structure to associate in the options.
         * @param[in] apRequest The request structure to associate in the options.
         * @return void
         *
         */
        virtual void
        ExaGeoStatOptionsInit(void *apOptions, void *apSequence, void *apRequest) = 0;

        /**
         * @brief Submit the release of the workspaces associated to the options structure.
         * @param[in,out] apOptions The options structure for which to workspaces will be released
         * @return void
         *
         */
        virtual void ExaGeoStatOptionsFree(void *apOptions) = 0;

        /**
         * @brief Finalize the runtime option structure for either HiCMA or CHAMELEON.
         * @param[in,out] apOptions The options structure that needs to be finalized.
         * @return void
         *
         */
        virtual void ExaGeoStatOptionsFinalize(void *apOptions) = 0;

        /**
         * @brief Get the pointer to the data or the runtime handler associated to the piece of data (m, n) in desc.
         * @param[in] apDescriptor The descriptor to which belongs the piece of data
         * @param[in] aDescRow The row coordinate of the piece of data in the matrix
         * @param[in] aDescCol The column coordinate of the piece of data in the matrix
         * @return void
         *
         */
        virtual void *ExaGeoStatDataGetAddr(void *apDescriptor, const int &aDescRow, const int &aDescCol) = 0;

        /**
        * @brief Get the number of tile rows of the sub-matrix
        * @param[in] apDescriptor
        * @return int
        *
        */
        virtual int GetMT(void *apDescriptor) = 0;

        /**
         * @brief Get the descriptor number of rows
         * @param[in] apDescriptor
         * @return int
         *
         */
        virtual int GetM(void *apDescriptor) = 0;

        /**
        * @brief Get the descriptor number of rows in a tile
        * @param[in] apDescriptor
        * @return int
        *
        */
        virtual int GetMB(void *apDescriptor) = 0;

        /**
         * @brief Get the descriptor options
         * @return void pointer to descriptor_option
         * @return void
         *
         */
        virtual void *GetOptions() = 0;

        /**
         * @brief Delete the options object
         * @param apOptions
         * @return void
         *
         */
        virtual void DeleteOptions(void *apOptions) = 0;

    };

}//namespace exageostat

#endif //EXAGEOSTATCPP_STARPUHELPERS_HPP
