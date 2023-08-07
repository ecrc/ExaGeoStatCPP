/**
 * @file DescriptorData.hpp
 * @brief 
 * @version 1.0.0
 * @author Sameh Abdulah
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
}
#endif

#include <common/Definitions.hpp>
#include <data-units/ExaGeoStatDescriptor.hpp>


namespace exageostat {
    namespace dataunits {

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
             * @brief Default constructor for DescriptorData.
             */
            DescriptorData() = default;

            /**
             * @brief Default destructor for DescriptorData.
             */
            virtual ~DescriptorData() = default;

            /**
             * @brief Getter for descriptors.
             * @param aDescriptorType Either Hicma or Chameleon descriptor.
             * @param aDescriptorName Name of the descriptor needed.
             * @return
             */
            BaseDescriptor
            GetDescriptor(common::DescriptorType aDescriptorType, common::DescriptorName aDescriptorName);

            /**
             * @brief Getter for mpSequence.
             * @return Sequence
             */
            void *GetSequence();

            /**
             * @brief Setter for mpSequence.
             * @param apSequence
             * @return void
             */
            void SetSequence(void *apSequence);

            /**
             * @brief Getter for mpRequest.
             * @return Request
             */
            void *GetRequest();

            /**
             * @brief Setter for mpRequest.
             * @param apRequest
             */
            void SetRequest(void *apRequest);

            /**
             * @brief setter for descriptors
             * @param aDescriptorType Either Hicma or Chameleon descriptor.
             * @param aDescriptorName Name of the descriptor to be set.
             * @param aIsOOC A boolean value indicating whether the matrix is out-of-core or not.
             * @param apMatrix A pointer to the beginning of the matrix.
             * @param aFloatPoint The precision of the matrix.
             * @param aMB The number of rows in a tile.
             * @param aNB The number of columns in a tile.
             * @param aSize The size of the matrix in elements including padding.
             * @param aLM The number of rows of the entire matrix.
             * @param aLN The number of columns of the entire matrix.
             * @param aI The row index to the beginning of the submatrix.
             * @param aJ The column index to the beginning of the submatrix.
             * @param aM The number of rows of the submatrix.
             * @param aN The number of columns of the submatrix.
             * @param aP The number of rows of the 2D distribution grid.
             * @param aQ The number of columns of the 2D distribution grid.
             * @return A pointer to the newly created CHAM_desc_t descriptor.
             */
            void
            SetDescriptor(common::DescriptorType aDescriptorType, common::DescriptorName aDescriptorName, bool aIsOOC,
                          void *apMatrix, common::FloatPoint aFloatPoint, int aMB, int aNB, int aSize, int aLM, int aLN,
                          int aI, int aJ, int aM, int aN, int aP, int aQ);

            /**
             * @brief Descriptors' Submatrix Creator
             * @param aDescriptorType Either Hicma or Chameleon descriptor.
             * @param aDescriptorName Name of the descriptor to be set.
             * @param apDescA The descriptor for which the submatrix will be created.
             * @param aI The row index to the beginning of the submatrix.
             * @param aJ The column index to the beginning of the submatrix.
             * @param aM The number of rows of the submatrix.
             * @param aN The number of columns of the submatrix.
             * @return
             */
            BaseDescriptor
            CreateSubMatrixDescriptor(common::DescriptorType aDescriptorType, common::DescriptorName aDescriptorName,
                                      void *apDescA, int aI, int aJ, int aM, int aN);

            /**
             * @brief Getter for the Descriptor matrix.
             * @param apDescriptor
             * @return pointer to the Descriptor matrix.
             */
            T *GetDescriptorMatrix(void *apDescriptor);

        private:
            /**
             * @brief Getter for Descriptor Name.
             * @param aDescriptorName
             * @return string of the descriptor's name.
             */
            std::string GetDescriptorName(common::DescriptorName aDescriptorName);

        private:
            /// Used Dictionary
            std::unordered_map<std::string, void *> mDictionary;
            /// Used Sequence
            void *mpSequence = nullptr;
            /// Used Request
            void *mpRequest = nullptr;
        };

        /**
        * @brief Instantiates the Linear Algebra methods class for float and double types.
        * @tparam T Data Type: float or double
        *
        */
        EXAGEOSTAT_INSTANTIATE_CLASS(DescriptorData)
    }//namespace configurations
}//namespace exageostat



#endif //EXAGEOSTATCPP_DESCRIPTORDATA_HPP
