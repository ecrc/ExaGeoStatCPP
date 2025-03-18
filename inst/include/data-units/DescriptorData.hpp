
// Copyright (c) 2017-2024 King Abdullah University of Science and Technology,
// All rights reserved.
// ExaGeoStat is a software package, provided by King Abdullah University of Science and Technology (KAUST).

/**
 * @file DescriptorData.hpp
 * @brief Contains the definition of the DescriptorData class.
 * @version 1.1.0
 * @author Mahmoud ElKarargy
 * @author Sameh Abdulah
 * @date 2023-07-18
**/

#ifndef EXAGEOSTATCPP_DESCRIPTORDATA_HPP
#define EXAGEOSTATCPP_DESCRIPTORDATA_HPP

#include <unordered_map>

#include <data-units/descriptor/ExaGeoStatDescriptor.hpp>
#include <hardware/ExaGeoStatHardware.hpp>

namespace exageostat::dataunits {

    /**
     * @brief Union represents the base descriptor.
     * @details This union is used to store different types of descriptors based on the configuration.
     *
     */
    union BaseDescriptor {
#if DEFAULT_RUNTIME
        CHAM_desc_t *chameleon_desc;
    #ifdef USE_HICMA
        HICMA_desc_t *hicma_desc;
    #endif
#else
        parsec_matrix_block_cyclic_t *parsec_desc;
#endif
    };

    /**
     * @Class DescriptorData
     * @brief Manages geo-statistical descriptor data with functions for retrieving and manipulating descriptors
     * @tparam T Data Type: float or double
     *
     */
    template<typename T>
    class DescriptorData {

    public:
        /**
         * @brief Default Constructor for DescriptorData.
         *
         */

        explicit DescriptorData() = default;

        /**
         * @brief Destructor for DescriptorData.
         *
         */
        ~DescriptorData();

        /**
         * @brief Get the base descriptor.
         * @param[in] aDescriptorType The type of the descriptor.
         * @param[in] aDescriptorName The name of the descriptor.
         * @return The base descriptor.
         * @throws std::runtime_error if the corresponding library is not enabled (USE_HICMA).
         *
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
         * @return void
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
         * @return void
         *
         */
        void SetRequest(void *apRequest);

#ifdef USE_HICMA

        /**
         * @brief Converts a CHAMELEON descriptor to a HICMA descriptor.
         * @param[in] apChameleonDesc Pointer to the CHAMELEON descriptor to be converted.
         * @param[in] aDescriptorName The name of the descriptor.
         * @return Pointer to the converted HICMA descriptor.
         *
         */
        HICMA_desc_t *
        ConvertChameleonToHicma(CHAM_desc_t *apChameleonDesc, const common::DescriptorName &aDescriptorName);

#endif

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
         * @param[in] aI The row index of the sub-matrix.
         * @param[in] aJ The column index of the sub-matrix.
         * @param[in] aM The number of rows in the sub-matrix.
         * @param[in] aN The number of columns in the sub-matrix.
         * @param[in] aP The number of rows in the complete matrix.
         * @param[in] aQ The number of columns in the complete matrix.
         * @param[in] aValidOOC Boolean refer to whether this descriptor can be created with OOC technology or not, default is true
         * @param[in] aConverted Boolean indicates whether it's a converted descriptor from chameleon to hicma or not, default is false
         * @return void
         * @throws std::runtime_error if the corresponding library is not enabled (USE_HICMA).
         *
         */
        void SetDescriptor(const common::DescriptorType &aDescriptorType, const common::DescriptorName &aDescriptorName,
                           const bool &aIsOOC = false, void *apMatrix = nullptr,
                           const common::FloatPoint &aFloatPoint = common::EXAGEOSTAT_REAL_DOUBLE, const int &aMB = 0,
                           const int &aNB = 0, const int &aSize = 0, const int &aLM = 0, const int &aLN = 0,
                           const int &aI = 0, const int &aJ = 0, const int &aM = 0, const int &aN = 0,
                           const int &aP = 0, const int &aQ = 0, const bool &aValidOOC = true, const bool &aConverted = false);

        /**
         * @brief Getter for the Descriptor matrix.
         * @param[in] aDescriptorType Type of the descriptor, whether it's CHAMELEON or HiCMA.
         * @param[in] aDescriptorName The name of the descriptor.
         * @return pointer to the Descriptor matrix.
         * @throws std::runtime_error if the corresponding library is not enabled (USE_HICMA).
         *
         */
        T *GetDescriptorMatrix(const common::DescriptorType &aDescriptorType,
                               const common::DescriptorName &aDescriptorName);

        /**
         * @brief Getter for the mIsDescriptorInitiated field.
         * @return mIsDescriptorInitiated
         *
         */
        bool GetIsDescriptorInitiated();

        /**
         * @brief Setter for mIsDescriptorInitiated field.
         * @param aIsInitiated Boolean for setting the mIsDescriptorInitiated field.
         *
         */
        void SetIsDescriptorInitiated(bool aIsInitiated);

    private:
        /**
         * @brief Get the descriptor name.
         * @param[in] aDescriptorName The descriptor name.
         * @return The descriptor name as a string.
         * @throws std::invalid_argument if the provided descriptor name is not available.
         *
         */
        std::string GetDescriptorName(const common::DescriptorName &aDescriptorName);

    private:
        //// Used Dictionary including the used descriptors.
        std::unordered_map<std::string, void *> mDictionary;
        //// Used sequence.
        void *mpSequence = nullptr;
        //// Used request
        void *mpRequest = nullptr;
        //// Specifies whether the descriptors have been initiated.
        bool mIsDescriptorInitiated = false;
    };

    /**
     * @brief Instantiates the DescriptorData class for float and double types.
     * @tparam T Data Type: float or double
     */
    EXAGEOSTAT_INSTANTIATE_CLASS(DescriptorData)
} // namespace exageostat


#endif //EXAGEOSTATCPP_DESCRIPTORDATA_HPP
