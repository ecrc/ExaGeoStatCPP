
// Copyright (c) 2017-2024 King Abdullah University of Science and Technology,
// All rights reserved.
// ExaGeoStat is a software package, provided by King Abdullah University of Science and Technology (KAUST).

/**
 * @file ExaGeoStatData.hpp
 * @brief Contains the definition of the ExaGeoStatData class.
 * @version 1.1.0
 * @author Mahmoud ElKarargy
 * @author Sameh Abdulah
 * @date 2024-02-04
**/

#ifndef EXAGEOSTATCPP_EXAGEOSTATDATA_HPP
#define EXAGEOSTATCPP_EXAGEOSTATDATA_HPP

#include <data-units/DescriptorData.hpp>
#include <data-units/Locations.hpp>

/**
 * @Class ExaGeoStatData
 * @brief  Manages geo-statistical data with functions for location and descriptor manipulation
 * @tparam T Data Type: float or double
 */
template<typename T>
class ExaGeoStatData {

public:
    /**
     * @brief Constructor for ExaGeoStatData.
     * @param[in] aSize The size of the data.
     * @param[in] aDimension The dimension of the data.
     *
     */
    ExaGeoStatData(const int &aSize, const exageostat::common::Dimension &aDimension);

    /**
     * @brief Constructor for ExaGeoStatData.
     * @param[in] aSize The size of the data.
     * @param[in] aDimension The dimension of the data.
     *
     */
    ExaGeoStatData(const int &aSize, const std::string &aDimension);

    /**
     * @brief Default constructor for ExaGeoStatData.
     *
     */
    ExaGeoStatData() = default;

    /**
     * @brief Destructor for ExaGeoStatData.
     *
     */
    ~ExaGeoStatData();

    /**
     * @brief Get the locations.
     * @return Pointer to the Locations object.
     *
     */
    exageostat::dataunits::Locations<T> *GetLocations();

    /**
     * @brief Set the locations.
     * @param[in] aLocation Pointer to the Locations object.
     * @return void
     *
     */
    void SetLocations(exageostat::dataunits::Locations<T> &aLocation);

    /**
     * @brief Get the descriptor data.
     * @return Pointer to the DescriptorData object.
     *
     */
    exageostat::dataunits::DescriptorData<T> *GetDescriptorData();

    /**
     * @brief Setter for the number of performed MLE iterations.
     * @param[in] aMleIterations number of performed MLE iterations.
     * @return void
     *
     */
    void SetMleIterations(const int &aMleIterations);

    /**
     * @brief Get the number of performed MLE iterations.
     * @return Pointer to the DescriptorData object.
     *
     */
    int GetMleIterations();

    /**
     * @brief Calculates Median Locations.
     * @param[in] aKernelName Name of the Kernel used.
     * @param[out] aLocations Location object to save medianLocations in.
     * @return void
     *
     */
    void CalculateMedianLocations(const std::string &aKernelName, exageostat::dataunits::Locations<T> &aLocations);

private:
    //// Used descriptor data.
    exageostat::dataunits::DescriptorData<T> *mpDescriptorData = nullptr;
    //// Used locations data.
    exageostat::dataunits::Locations<T> *mpLocations = nullptr;
    //// Current number of performed MLE iterations.
    int mMleIterations = 0;
};

/**
 * @brief Instantiates the ExaGeoStatData class for float and double types.
 * @tparam T Data Type: float or double
 */
EXAGEOSTAT_INSTANTIATE_CLASS(ExaGeoStatData)

#endif //EXAGEOSTATCPP_EXAGEOSTATDATA_HPP