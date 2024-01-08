
// Copyright (c) 2017-2023 King Abdullah University of Science and Technology,
// All rights reserved.
// ExaGeoStat is a software package, provided by King Abdullah University of Science and Technology (KAUST).

/**
 * @file ExaGeoStatHardware.hpp
 * @brief Contains the definition of the ExaGeoStatHardware class.
 * @version 1.0.0
 * @author Mahmoud ElKarargy
 * @author Sameh Abdulah
 * @date 2023-08-07
**/

#ifndef EXAGEOSTATCPP_EXAGEOSTATHARDWARE_HPP
#define EXAGEOSTATCPP_EXAGEOSTATHARDWARE_HPP

#include <common/Definitions.hpp>

namespace exageostat::hardware {

    /**
     * @brief Class representing the hardware configuration for the ExaGeoStat solver.
     */
    class ExaGeoStatHardware {

    public:
        /**
         * @brief Constructor for ExaGeoStatHardware.
         * @param[in] aComputation The computation mode for the solver.
         * @param[in] aCoreNumber The number of CPU cores to use for the solver.
         * @param[in] aGpuNumber The number of GPUs to use for the solver.
         *
         */
        ExaGeoStatHardware(const common::Computation &aComputation, const int &aCoreNumber, const int &aGpuNumber);

        /**
         * @brief Destructor for ExaGeoStatHardware.
         */
        virtual ~ExaGeoStatHardware();

        /**
         * @brief Get the Chameleon hardware context.
         * @return Pointer to the hardware context.
         *
         */
        [[nodiscard]] void *GetChameleonContext() const;

#ifdef USE_HICMA

/**
             * @brief Get the Hicma hardware context.
             * @return Pointer to the hardware context.
             *
             */
        [[nodiscard]] void *GetHicmaContext() const;

#endif

        /**
         * @brief Get the hardware context.
         * @param[in] aComputation Used computation to decide whether to use Hicma or Chameleon context.
         * @return Pointer to the hardware context.
         *
         */
        [[nodiscard]] void *GetContext(common::Computation aComputation) const;

    private:
        //// Used Pointer to the Chameleon hardware context.
        void *mpChameleonContext = nullptr;
#ifdef USE_HICMA
        //// Used Pointer to the Hicma hardware context.
        void *mpHicmaContext = nullptr;
#endif
        //// Used Computation mode for the solver.
        common::Computation mComputation;
    };
} // namespace exageostat

#endif // EXAGEOSTATCPP_EXAGEOSTATHARDWARE_HPP