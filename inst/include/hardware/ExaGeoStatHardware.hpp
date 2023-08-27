
// Copyright (c) 2017-2023 King Abdullah University of Science and Technology,
// All rights reserved.
// ExaGeoStat is a software package, provided by King Abdullah University of Science and Technology (KAUST).

/**
 * @file ExaGeoStatHardware.hpp
 * @brief Contains the definition of the ExaGeoStatHardware class.
 * @version 1.0.0
 * @author Sameh Abdulah
 * @author Mahmoud ElKarargy
 * @date 2023-08-07
**/

#ifndef EXAGEOSTATCPP_EXAGEOSTATHARDWARE_HPP
#define EXAGEOSTATCPP_EXAGEOSTATHARDWARE_HPP

#include <common/Definitions.hpp>

namespace exageostat {
    namespace hardware {

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
             * @brief Get the hardware context.
             * @return Pointer to the hardware context.
             *
             */
            void *GetContext() const;

        private:
            //// Used Pointer to the hardware context.
            void *mpContext = nullptr;
            //// Used Computation mode for the solver.
            common::Computation mComputation;
        };
    } // namespace hardware
} // namespace exageostat

#endif // EXAGEOSTATCPP_EXAGEOSTATHARDWARE_HPP