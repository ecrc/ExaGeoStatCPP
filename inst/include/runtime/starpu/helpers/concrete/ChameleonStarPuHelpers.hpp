
// Copyright (c) 2017-2024 King Abdullah University of Science and Technology,
// All rights reserved.
// ExaGeoStat is a software package, provided by King Abdullah University of Science and Technology (KAUST).

/**
 * @file ChameleonStarPuHelpers.hpp
 * @brief A class for Chameleon implementation of StarPu helpers interface StarPuHelpers.hpp.
 * @version 1.1.0
 * @author Mahmoud ElKarargy
 * @date 2024-02-25
**/

#ifndef EXAGEOSTATCPP_CHAMELEONSTARPUHELPERS_HPP
#define EXAGEOSTATCPP_CHAMELEONSTARPUHELPERS_HPP

#include <runtime/starpu/helpers/StarPuHelpers.hpp>

namespace exageostat::runtime {

    /**
     * @brief ChameleonStarPuHelpers is a concrete implementation of StarPuHelpers interface for Chameleon library.
     *
     */
    class ChameleonStarPuHelpers : public StarPuHelpers {
    public:

        /**
        * @brief Default constructor.
        */
        ChameleonStarPuHelpers() = default;

        /**
         * @brief Default destructor.
         */
        ~ChameleonStarPuHelpers() = default;

        /**
        * @brief Initialize the runtime option structure for CHAMELEON
        * @copydoc StarPuHelpers::ExaGeoStatOptionsInit()
        *
        */
        void
        ExaGeoStatOptionsInit(void *apOptions, void *apSequence, void *apRequest) override;

        /**
         * @brief Submit the release of the workspaces associated to the options structure.
         * @copydoc StarPuHelpers::ExaGeoStatOptionsFree()
         *
         */
        void
        ExaGeoStatOptionsFree(void *apOptions) override;

        /**
         * @brief Finalize the runtime option structure for CHAMELEON.
         * @copydoc StarPuHelpers::ExaGeoStatOptionsFinalize()
         *
         */
        void
        ExaGeoStatOptionsFinalize(void *apOptions) override;

        /**
         * @brief Get the pointer to the data or the runtime handler associated to the piece of data (m, n) in desc.
         * @copydoc StarPuHelpers::ExaGeoStatDataGetAddr()
         *
         */
        void *ExaGeoStatDataGetAddr(void *apDescriptor, const int &aDescRow, const int &aDescCol) override;

        /**
        * @brief Get the number of tile rows of the sub-matrix
        * @copydoc StarPuHelpers::GetMT()
        *
        */
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
         * @copydoc StarPuHelpers::GetOptions()
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

#endif //EXAGEOSTATCPP_CHAMELEONSTARPUHELPERS_HPP
