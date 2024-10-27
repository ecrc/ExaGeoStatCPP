
// Copyright (c) 2017-2024 King Abdullah University of Science and Technology,
// All rights reserved.
// ExaGeoStat is a software package, provided by King Abdullah University of Science and Technology (KAUST).

/**
 * @file ExaGeoStatHardware.hpp
 * @brief Contains the definition of the ExaGeoStatHardware class.
 * @version 1.1.0
 * @author Mahmoud ElKarargy
 * @author Sameh Abdulah
 * @date 2024-01-24
**/

#ifndef EXAGEOSTATCPP_EXAGEOSTATHARDWARE_HPP
#define EXAGEOSTATCPP_EXAGEOSTATHARDWARE_HPP

#include <common/Definitions.hpp>
#include <configurations/Configurations.hpp>

/**
 * @brief Class represents the hardware configuration for the ExaGeoStat solver.
 *
 */
class ExaGeoStatHardware {

public:

    /**
     * @brief Constructor for ExaGeoStatHardware.
     * @param[in] aConfigurations The set of arguments from the configurations.
     *
     */
    explicit ExaGeoStatHardware(exageostat::configurations::Configurations &aConfigurations);

    /**
     * @brief Constructor for ExaGeoStatHardware.
     * @param[in] aComputation The computation mode for the solver.
     * @param[in] aCoreNumber The number of CPU cores to use for the solver.
     * @param[in] aGpuNumber The number of GPUs to use for the solver.
     * @param[in] aP The P grid dimension setting, default is 1.
     * @param[in] aQ The Q grid dimension setting, default is 1.
     *
     */
    explicit ExaGeoStatHardware(const exageostat::common::Computation &aComputation, const int &aCoreNumber,
                                const int &aGpuNumber, const int &aP = 1, const int &aQ = 1);

    /**
     * @brief Constructor for ExaGeoStatHardware.
     * @param[in] aComputation The computation mode for the solver as a string.
     * @param[in] aCoreNumber The number of CPU cores to use for the solver.
     * @param[in] aGpuNumber The number of GPUs to use for the solver.
     * @param[in] aP The P grid dimension setting.
     * @param[in] aQ The Q grid dimension setting.
     *
     */
    explicit ExaGeoStatHardware(const std::string &aComputation, const int &aCoreNumber, const int &aGpuNumber,
                                const int &aP = 1, const int &aQ = 1);

    /**
     * @brief A Finalize caller for Hardware.
     * @return void.
     *
     */
    void FinalizeHardware();

    /**
     * @brief Destructor for ExaGeoStatHardware.
     *
     */
    ~ExaGeoStatHardware();

    /**
     * @brief Initializes hardware configuration.
     * @param[in] aComputation The computation mode for the solver.
     * @param[in] aCoreNumber The number of CPU cores to use for the solver.
     * @param[in] aGpuNumber The number of GPUs to use for the solver.
     * @param[in] aP The P grid dimension setting.
     * @param[in] aQ The Q grid dimension setting.
     * @return void
     *
     */
    static void
    InitHardware(const exageostat::common::Computation &aComputation, const int &aCoreNumber, const int &aGpuNumber,
                 const int &aP, const int &aQ);

    /**
     * @brief Get the Chameleon hardware context.
     * @return Pointer to the hardware context.
     *
     */
    [[nodiscard]] static void *GetChameleonContext();

    /**
     * @brief Get the HiCMA hardware context.
     * @return Pointer to the hardware context.
     *
     */
    [[nodiscard]] static void *GetHicmaContext();

    /**
     * @brief Get the PaRSEC hardware context.
     * @return Pointer to the hardware context.
     *
     */
    [[nodiscard]] static void *GetParsecContext();

    /**
     * @brief Get the hardware context.
     * @param[in] aComputation Used computation to decide whether to use Hicma or Chameleon context.
     * @return Pointer to the hardware context.
     *
     */
    [[nodiscard]] static void *GetContext(exageostat::common::Computation aComputation);

    /**
     * @brief Sets the rank of MPI for PaRSEC.
     * @param[in] aRank The new value for the rank.
     * @return void
     *
    **/
    static void SetParsecMPIRank(int aRank);

    /**
     * @brief Retrieves the rank of MPI for PaRSEC.
     * @return The current rank of MPI PaRSEC.
     *
    **/
    static int GetParsecMPIRank();

    /**
     * @brief Retrieves the P dimension of the grid.
     * @details This function returns the current setting of the P dimension of the grid, which is part of the grid configuration used in various computational processes.
     * @return The current P dimension setting.
     *
    **/
    static int GetPGrid();

    /**
     * @brief Retrieves the Q dimension of the grid.
     * @details This function returns the current setting of the Q dimension of the grid, which is part of the grid configuration used in various computational processes.
     * @return int The current Q dimension setting.
     *
    **/
    static int GetQGrid();

    /**
     * @brief Sets the P dimension of the grid.
     * @details This function updates the P dimension setting of the grid. This dimension is critical in configuring the grid's layout for simulations or calculations.
     * @param[in] aP The new value for the P dimension.
     * @return void
     *
    **/
    static void SetPGrid(int aP);

    /**
     * @brief Sets the Q dimension of the grid.
     * @details This function updates the Q dimension setting of the grid. This dimension is crucial in configuring the grid's layout for simulations or calculations.
     * @param[in] aQ The new value for the Q dimension.
     * @return void
     *
    **/
    static void SetQGrid(int aQ);

private:
    //// Used Pointer to the Chameleon hardware context.
    static void *mpChameleonContext;
    //// Used Pointer to the Hicma hardware context.
    static void *mpHicmaContext;
    //// Used Pointer to the PaRSEC hardware context.
    static void *mpParsecContext;
    //// Used P-Grid
    static int mParsecMPIRank;
    //// Used P-Grid
    static int mPGrid;
    //// Used Q-Grid
    static int mQGrid;
    //// Used boolean to avoid re-init mpi
    static bool mIsMPIInit;
};

#endif // EXAGEOSTATCPP_EXAGEOSTATHARDWARE_HPP