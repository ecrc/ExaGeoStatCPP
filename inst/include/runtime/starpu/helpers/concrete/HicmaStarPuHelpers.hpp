
// Copyright (c) 2017-2023 King Abdullah University of Science and Technology,
// All rights reserved.
// ExaGeoStat is a software package, provided by King Abdullah University of Science and Technology (KAUST).

/**
 * @file HicmaStarPuHelpers.hpp
 * @brief A class for Hicma implementation of StarPu helpers interface StarPuHelpers.hpp.
 * @version 1.1.0
 * @author Mahmoud ElKarargy
 * @date 2024-02-25
**/

#ifndef EXAGEOSTATCPP_HICMASTARPUHELPERS_HPP
#define EXAGEOSTATCPP_HICMASTARPUHELPERS_HPP

#include <runtime/starpu/helpers/StarPuHelpers.hpp>

namespace exageostat::runtime {

    /**
     * @brief HicmaStarPuHelpers is a concrete implementation of StarPuHelpers interface for Hicma library.
     *
     */
    class HicmaStarPuHelpers : public StarPuHelpers {
    public:

        /**
        * @brief Default constructor.
        *
        */
        HicmaStarPuHelpers() = default;

        /**
         * @brief Default destructor.
         *
         */
        ~HicmaStarPuHelpers() = default;

        /**
        * @brief Initialize the runtime option structure for HiCMA
        * @copydoc StarPuHelpers::ExaGeoStatOptionsInit()
        *
        */
        void ExaGeoStatOptionsInit(void *apOptions, void *apSequence, void *apRequest) override;

        /**
         * @brief Submit the release of the workspaces associated to the options structure.
         * @copydoc StarPuHelpers::ExaGeoStatOptionsFree()
         *
         */
        void ExaGeoStatOptionsFree(void *apOptions) override;

        /**
         * @brief Finalize the runtime option structure for HiCMA.
         * @copydoc StarPuHelpers::ExaGeoStatOptionsFinalize()
         *
         */
        void ExaGeoStatOptionsFinalize(void *apOptions) override;

        /**
         * @brief Get the pointer to the data or the runtime handler associated to the piece of data (m, n) in desc.
         * @copydoc StarPuHelpers::ExaGeoStatDataGetAddr()
         *
         */
        void *ExaGeoStatDataGetAddr(void *apDescriptor, const int& aDescRow, const int& aDescCol) override;

        /**
         *  @brief Get the number of tile rows of the sub-matrix
         *  @copydoc StarPuHelpers::GetMT()
         *
         **/
        int GetMT(void *apDescriptor) override;

        /**
         * @brief Get the descriptor number of rows
         * @copydoc StarPuHelpers::GetM()
         *
         */
        int GetM(void *apDescriptor) override;

        /**
        * @brief Get the descriptor number of rows in a tile
        * @copydoc StarPuHelpers::GetMB()
        *
        */
        int GetMB(void *apDescriptor) override;

        /**
         * @brief Get the descriptor options
         * @copydoc StarPuHelpers::GetOptions)
         * 
         */
        void *GetOptions() override;

        /**
        * @brief Delete the options object
        * @copydoc StarPuHelpers::DeleteOptions()
        *
        */
        void DeleteOptions(void *apOptions) override;

    };

}//namespace exageostat

#endif //EXAGEOSTATCPP_HICMASTARPUHELPERS_HPP
