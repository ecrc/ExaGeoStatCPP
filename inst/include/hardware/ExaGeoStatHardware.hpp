
// Copyright (c) 2017-2023 King Abdullah University of Science and Technology,
// All rights reserved.
// ExaGeoStat is a software package, provided by King Abdullah University of Science and Technology (KAUST).

/**
 * @file ExaGeoStatHardware.hpp
 * @brief Contains the definition of the ExaGeoStatHardware class.
 * @version 1.1.0
 * @author Mahmoud ElKarargy
 * @author Sameh Abdulah
 * @date 2023-01-20
**/

#ifndef EXAGEOSTATCPP_EXAGEOSTATHARDWARE_HPP
#define EXAGEOSTATCPP_EXAGEOSTATHARDWARE_HPP

#include <common/Definitions.hpp>

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
    ExaGeoStatHardware(const exageostat::common::Computation &aComputation, const int &aCoreNumber, const int &aGpuNumber);

    /**
     * @brief Constructor for ExaGeoStatHardware.
     * @param[in] aComputation The computation mode for the solver as a string.
     * @param[in] aCoreNumber The number of CPU cores to use for the solver.
     * @param[in] aGpuNumber The number of GPUs to use for the solver.
     */
    ExaGeoStatHardware(const std::string &aComputation, const int &aCoreNumber, const int &aGpuNumber);

    /**
     * @brief Destructor for ExaGeoStatHardware.
     */
    virtual ~ExaGeoStatHardware();

    /**
     * @brief Initializes hardware configuration.
     * @param[in] aComputation The computation mode for the solver.
     * @param[in] aCoreNumber The number of CPU cores to use for the solver.
     * @param[in] aGpuNumber The number of GPUs to use for the solver.
     */
    void InitHardware(const exageostat::common::Computation &aComputation, const int &aCoreNumber, const int &aGpuNumber);

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
    [[nodiscard]] void *GetContext(exageostat::common::Computation aComputation) const;

private:
    //// Used Pointer to the Chameleon hardware context.
    void *mpChameleonContext = nullptr;
    //// Used Pointer to the Hicma hardware context.
    void *mpHicmaContext = nullptr;
};

#endif // EXAGEOSTATCPP_EXAGEOSTATHARDWARE_HPP