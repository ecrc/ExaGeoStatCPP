
// Copyright (c) 2017-2024 King Abdullah University of Science and Technology,
// All rights reserved.
// ExaGeoStat is a software package, provided by King Abdullah University of Science and Technology (KAUST).

/**
 * @file CommunicatorMPI.hpp
 * @brief Defines the CommunicatorMPI class for MPI rank communication.
 * @version 1.1.0
 * @author Sameh Abdulah
 * @date 2023-11-10
**/

#ifndef EXAGEOSTATCPP_COMMUNICATORMPI_HPP
#define EXAGEOSTATCPP_COMMUNICATORMPI_HPP

namespace exageostat::helpers {

    /**
     * @class CommunicatorMPI
     * @brief A class for Communicating MPI rank.
     * @details The CommunicatorMPI class provides functionality to communicate MPI rank information.
     */
    class CommunicatorMPI {
    public:

        /**
         * @brief Get a pointer to the singleton instance of the CommunicatorMPI class.
         * @return A pointer to the instance of the CommunicatorMPI class.
         *
         */
        static CommunicatorMPI *GetInstance();

        /**
         * @brief Get the rank of the MPI process.
         * @return The rank of the MPI process.
         *
         */
        [[nodiscard]] int GetRank() const;

        /**
         * @brief Set the hardware initialization flag.
         * @details This function sets the flag to indicate that hardware has been initialized.
         * @return void
         *
         */
        void SetHardwareInitialization();

        /**
         * @brief Unset the hardware initialization flag.
         * @details This function remove the flag to indicate that hardware has been initialized.
         * @return void
         *
         */
        void RemoveHardwareInitialization();

    private:

        /**
         * @brief Prevent Class Instantiation for Communicator MPI Class.
         *
         */
        CommunicatorMPI() = default;

        /**
         * @brief Pointer to the singleton instance of the CommunicatorMPI class.
         *
         */
        static CommunicatorMPI *mpInstance;
        /// Used boolean to check if hardware is initialized.
        bool mIsHardwareInitialized;
    };
}
#endif //EXAGEOSTATCPP_COMMUNICATORMPI_HPP
