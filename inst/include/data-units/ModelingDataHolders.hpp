/**
 * @file ModelingDataHolders.hpp
 * @brief This file contains the definition of the mModelingData struct, which contains all the data needed for modeling.
 * @version 1.0.0
 * @author Sameh Abdulah
 * @date 2023-08-24
**/

#ifndef EXAGEOSTATCPP_MODELINGDATAHOLDERS_HPP
#define EXAGEOSTATCPP_MODELINGDATAHOLDERS_HPP


namespace exageostat {
    namespace dataunits {

        /**
         * @brief Struct containing all the data needed for modeling.
         * @tparam T The data type of the data.
         */
        template<typename T>
        struct mModelingData {
            /// ExaGeoStatData<T> object containing needed descriptors, and locations.
            dataunits::ExaGeoStatData<T> *mpData;
            /// Configurations object containing user input data.
            configurations::Configurations *mpConfiguration;
            /// Hardware configuration for the ExaGeoStat solver.
            const hardware::ExaGeoStatHardware *mpHardware;
            /// User Input Measurements Matrix
            T *mpMeasurementsMatrix;

            /**
             * @brief Constructor.
             * @param data The ExaGeoStatData object.
             * @param configuration The Configurations object.
             * @param hardware The hardware configuration object.
             */
            mModelingData(dataunits::ExaGeoStatData<T> *data, configurations::Configurations *configuration,
                          const hardware::ExaGeoStatHardware *hardware, T *matrix) : mpData(std::move(data)),
                                                                                     mpConfiguration(configuration),
                                                                                     mpHardware(hardware),
                                                                                     mpMeasurementsMatrix(matrix) {}
        };

    }//namespace dataunits
}//namespace exageostat
#endif //EXAGEOSTATCPP_MODELINGDATAHOLDERS_HPP
