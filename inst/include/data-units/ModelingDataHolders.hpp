/**
 * @file ModelingDataHolders.hpp
 * @brief This file contains the definition of the mModelingData struct, which contains all the data needed for modeling.
 * @version 1.1.0
 * @author Mahmoud ElKarargy
 * @date 2023-08-24
**/

#ifndef EXAGEOSTATCPP_MODELINGDATAHOLDERS_HPP
#define EXAGEOSTATCPP_MODELINGDATAHOLDERS_HPP

namespace exageostat::dataunits {

    /**
     * @brief Struct containing all the data needed for modeling.
     * @tparam T The data type of the data.
     */
    template<typename T>
    struct mModelingData {
        /// ExaGeoStatData<T> object containing needed descriptors, and locations.
        std::unique_ptr<ExaGeoStatData<T>> *mpData;
        /// Configurations object containing user input data.
        Configurations *mpConfiguration;
        /// Hardware configuration for the ExaGeoStat solver.
        const ExaGeoStatHardware *mpHardware;
        /// Used Kernel for ExaGeoStat Modeling Data.
        const kernels::Kernel<T> *mpKernel;

        /// User Input Measurements Matrix
        T *mpMeasurementsMatrix;

        /**
         * @brief Constructor.
         * @param aData The ExaGeoStatData object.
         * @param aConfiguration The Configurations object.
         * @param aHardware The hardware configuration object.
         * @param aKernel The Kernel object.
         */
        mModelingData(std::unique_ptr<ExaGeoStatData<T>> &aData, Configurations &aConfiguration,
                      const ExaGeoStatHardware &aHardware, T &aMatrix, const kernels::Kernel<T> &aKernel) :
                mpData(std::move(&aData)), mpConfiguration(&aConfiguration), mpHardware(&aHardware),
                mpMeasurementsMatrix(&aMatrix), mpKernel(&aKernel) {}
    };

}//namespace exageostat
#endif //EXAGEOSTATCPP_MODELINGDATAHOLDERS_HPP