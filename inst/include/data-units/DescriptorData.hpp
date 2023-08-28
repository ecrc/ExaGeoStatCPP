
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

#include <linear-algebra-solvers/concrete/ChameleonHeaders.hpp>
#include <linear-algebra-solvers/concrete/HicmaHeaders.hpp>
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
#ifdef EXAGEOSTAT_USE_HICMA
            HICMA_desc_t *hicma_desc;
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
             * @param[in] aHardware Reference to the ExaGeoStatHardware object.
             * @throws std::runtime_error if hardware is not initialized.
             */

            explicit DescriptorData(const hardware::ExaGeoStatHardware &aHardware);

            /**
             * @brief Destructor for DescriptorData.
             */
            virtual ~DescriptorData();

            /**
             * @brief Get the base descriptor.
             * @param[in] aDescriptorType The type of the descriptor.
             * @param[in] aDescriptorName The name of the descriptor.
             * @return The base descriptor.
             * @throws std::runtime_error if the corresponding library is not enabled (EXAGEOSTAT_USE_CHAMELEON or EXAGEOSTAT_USE_HICMA).
             */
            BaseDescriptor
            GetDescriptor(const common::DescriptorType &aDescriptorType, const common::DescriptorName &aDescriptorName);

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
             * @throws std::runtime_error if the corresponding library is not enabled (EXAGEOSTAT_USE_CHAMELEON or EXAGEOSTAT_USE_HICMA).
             */
            void
            SetDescriptor(const common::DescriptorType &aDescriptorType, const common::DescriptorName &aDescriptorName,
                          const bool &aIsOOC, void *apMatrix, const common::FloatPoint &aFloatPoint, const int &aMB,
                          const int &aNB, const int &aSize, const int &aLM, const int &aLN, const int &aI,
                          const int &aJ, const int &aM, const int &aN, const int &aP, const int &aQ);

            /**
             * @brief Getter for the Descriptor matrix.
             * @param[in] aDescriptorType Type of the descriptor, whether it's CHAMELEON or HiCMA.
             * @param[in] apDescriptor Pointer to the descriptor.
             * @return pointer to the Descriptor matrix.
             * @throws std::runtime_error if the corresponding library is not enabled (EXAGEOSTAT_USE_CHAMELEON or EXAGEOSTAT_USE_HICMA).
             */
            T *GetDescriptorMatrix(const common::DescriptorType &aDescriptorType, void *apDescriptor);

            /**
             * @brief Getter for the mIsDescriptorInitiated field.
             * @return mIsDescriptorInitiated
             */
            bool GetIsDescriptorInitiated();

            /**
             * @brief Setter for mIsDescriptorInitiated field.
             * @param aIsInitiated Boolean for setting the mIsDescriptorInitiated field.
             */
            void SetIsDescriptorInitiated(bool aIsInitiated);

        private:
            /**
             * @brief Get the descriptor name.
             * @param[in] aDescriptorName The descriptor name.
             * @return The descriptor name as a string.
             * @throws std::invalid_argument if the provided descriptor name is not available.
             */
            std::string GetDescriptorName(const common::DescriptorName &aDescriptorName);

        private:
            //// Used Dictionary including the used descriptors.
            std::unordered_map<std::string, void *> mDictionary;
            //// Used sequence.
            void *mpSequence = nullptr;
            //// Used request
            void *mpRequest = nullptr;
            //// Used context
            void *mpContext = nullptr;
            //// Specifies whether the descriptors have been initiated.
            bool mIsDescriptorInitiated = false;
        };

        /**
         * @brief Instantiates the DescriptorData class for float and double types.
         * @tparam T Data Type: float or double
         */
        EXAGEOSTAT_INSTANTIATE_CLASS(DescriptorData)
    } // namespace dataunits
} // namespace exageostat


#endif //EXAGEOSTATCPP_DESCRIPTORDATA_HPP
