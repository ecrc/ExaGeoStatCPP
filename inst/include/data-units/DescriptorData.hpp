
// Copyright (c) 2017-2023 King Abdullah University of Science and Technology,
// All rights reserved.
// ExaGeoStat is a software package, provided by King Abdullah University of Science and Technology (KAUST).

/**
 * @file DescriptorData.hpp
 * @brief Contains the definition of the DescriptorData class.
 * @version 1.0.0
 * @author Sameh Abdulah
 * @author Mahmoud ElKarargy
 * @date 2023-07-18
**/

#ifndef EXAGEOSTATCPP_DESCRIPTORDATA_HPP
#define EXAGEOSTATCPP_DESCRIPTORDATA_HPP

#include <iostream>
#include <any>
#include <vector>
#include <unordered_map>

#ifdef EXAGEOSTAT_USE_CHAMELEON
extern "C" {
#include <chameleon/struct.h>
#include <chameleon.h>
}
#endif
#ifdef EXAGEOSTAT_USE_HiCMA
extern "C"{
#include <hicma_struct.h>
#include <hicma.h>
}
#endif

#include <common/Definitions.hpp>
#include <data-units/descriptor/ExaGeoStatDescriptor.hpp>
#include <hardware/ExaGeoStatHardware.hpp>

namespace exageostat {
    namespace dataunits {

        /**
         * @brief Union representing the base descriptor.
         * @details This union is used to store different types of descriptors based on the configuration.
         */
        union BaseDescriptor {
#ifdef EXAGEOSTAT_USE_CHAMELEON
            CHAM_desc_t *chameleon_desc;
#endif
#ifdef EXAGEOSTAT_USE_HiCMA
            HICMA_desc_t* hicma_desc;
#endif
        };

        /**
         * @Class DescriptorData
         * @brief Manages geo-statistical descriptor data with functions for retrieving and manipulating descriptors
         * @tparam T Data Type: float or double
         */
        template<typename T>
        class DescriptorData {

        public:
            /**
             * @brief Constructor for DescriptorData.
             * @param[in] apHardware Reference to the ExaGeoStatHardware object.
             * @throws std::runtime_error if hardware is not initialized.
             */

            explicit DescriptorData(hardware::ExaGeoStatHardware &apHardware);

            /**
             * @brief Destructor for DescriptorData.
             */
            virtual ~DescriptorData();

            /**
             * @brief Get the base descriptor.
             * @param[in] aDescriptorType The type of the descriptor.
             * @param[in] aDescriptorName The name of the descriptor.
             * @return The base descriptor.
             * @throws std::runtime_error if the corresponding library is not enabled (EXAGEOSTAT_USE_CHAMELEON or EXAGEOSTAT_USE_HiCMA).
             */
            BaseDescriptor
            GetDescriptor(common::DescriptorType aDescriptorType, common::DescriptorName aDescriptorName);

            /**
             * @brief Get the sequence.
             * @return Pointer to the sequence.
             *
             */
            void *GetSequence();

            /**
             * @brief Set the sequence.
             * @param[in] apSequence Pointer to the sequence.
             *
             */
            void SetSequence(void *apSequence);

            /**
             * @brief Get the request.
             * @return Pointer to the request.
             *
             */
            void *GetRequest();

            /**
             * @brief Set the request.
             * @param[in] apRequest Pointer to the request.
             *
             */
            void SetRequest(void *apRequest);

            /**
             * @brief Set the descriptor.
             * @param[in] aDescriptorType The type of the descriptor.
             * @param[in] aDescriptorName The name of the descriptor.
             * @param[in] aIsOOC Boolean indicating if the descriptor is out-of-core.
             * @param[in] apMatrix Pointer to the matrix.
             * @param[in] aFloatPoint The floating-point precision.
             * @param[in] aMB The number of rows in a block.
             * @param[in] aNB The number of columns in a block.
             * @param[in] aSize The size of the matrix.
             * @param[in] aLM The leading dimension of the matrix.
             * @param[in] aLN The trailing dimension of the matrix.
             * @param[in] aI The row index of the submatrix.
             * @param[in] aJ The column index of the submatrix.
             * @param[in] aM The number of rows in the submatrix.
             * @param[in] aN The number of columns in the submatrix.
             * @param[in] aP The number of rows in the complete matrix.
             * @param[in] aQ The number of columns in the complete matrix.
             * @return void
             * @throws std::runtime_error if the corresponding library is not enabled (EXAGEOSTAT_USE_CHAMELEON or EXAGEOSTAT_USE_HiCMA).
             */
            void
            SetDescriptor(common::DescriptorType aDescriptorType, common::DescriptorName aDescriptorName, bool aIsOOC,
                          void *apMatrix, common::FloatPoint aFloatPoint, int aMB, int aNB, int aSize, int aLM, int aLN,
                          int aI, int aJ, int aM, int aN, int aP, int aQ);

            /**
             * @brief Getter for the Descriptor matrix.
             * @param[in] aDescriptorType Type of the descriptor, whether it's CHAMELEON or HiCMA.
             * @param[in] apDescriptor Pointer to the descriptor.
             * @return pointer to the Descriptor matrix.
             * @throws std::runtime_error if the corresponding library is not enabled (EXAGEOSTAT_USE_CHAMELEON or EXAGEOSTAT_USE_HiCMA).
             */
            T *GetDescriptorMatrix(common::DescriptorType aDescriptorType, void *apDescriptor);

        private:
            /**
             * @brief Get the descriptor name.
             * @param[in] aDescriptorName The descriptor name.
             * @return The descriptor name as a string.
             * @throws std::invalid_argument if the provided descriptor name is not available.
             */
            std::string GetDescriptorName(common::DescriptorName aDescriptorName);

        private:
            //// Used Dictionary including the used descriptors.
            std::unordered_map<std::string, void *> mDictionary;
            //// Used sequence.
            void *mpSequence = nullptr;
            //// Used request
            void *mpRequest = nullptr;
            //// Used context
            void *mpContext = nullptr;
        };

        /**
         * @brief Instantiates the DescriptorData class for float and double types.
         * @tparam T Data Type: float or double
         */
        EXAGEOSTAT_INSTANTIATE_CLASS(DescriptorData)
    } // namespace dataunits
} // namespace exageostat


#endif //EXAGEOSTATCPP_DESCRIPTORDATA_HPP
