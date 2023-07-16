
/*
 * Copyright (c) 2017-2023 King Abdullah University of Science and Technology,
 * All rights reserved.
 * ExaGeoStat is a software package, provided by King Abdullah University of Science and Technology (KAUST).
 */

/**
 * @file DataGenerator.hpp
 * @brief Contains the declaration of the SyntheticDataConfigurations class.
 * @version 1.0.0
 * @author Sameh Abdulah
 * @date 2023-02-14
**/

#ifndef EXAGEOSTAT_CPP_DATAGENERATOR_HPP
#define EXAGEOSTAT_CPP_DATAGENERATOR_HPP

#include <memory>

#include <data-units/Locations.hpp>
#include <configurations/data-generation/concrete/SyntheticDataConfigurations.hpp>
#include <linear-algebra-solvers/LinearAlgebraMethods.hpp>
#include <linear-algebra-solvers/LinearAlgebraFactory.hpp>
#include <kernels/Kernel.hpp>

namespace exageostat {
    namespace generators {

        /**
         * @class DataGenerator
         * @brief Abstract base class for generating synthetic or real data.
         * @tparam T Data Type: float or double
         *
         */
        template<typename T>
        class DataGenerator {

        public:

            /**
             * @brief Generates the data locations.
             * @details This method generates the X, Y, and Z variables used to define the locations of the data points.
             * @return void
             *
             */
            virtual void
            GenerateLocations() = 0;

            /**
             * @brief Generates the data descriptors.
             * @details This method generates the descriptors used to define the properties of the data points.
             * @return void
             *
             */
            virtual void
            GenerateDescriptors() = 0;

            /**
             * @brief Destroys the data descriptors.
             * @details This method frees the memory used by the data descriptors.
             * @return void
             *
             */
            virtual void DestoryDescriptors() = 0;

            /**
             * @brief Generates the data observations.
             * @details This method generates the observations of the data points, which are used to train and test the model.
             * @return void
             *
             */
            virtual void
            GenerateObservations() = 0;

            /**
             * @brief Factory method for creating a data generator object.
             * @details This method creates a data generator object based on the specified configurations.
             * @param[in] apConfigurations Pointer to the synthetic data configurations.
             * @return A unique pointer to the created data generator object.
             *
             */
            static std::unique_ptr<DataGenerator>
            CreateGenerator(configurations::data_configurations::SyntheticDataConfigurations *apConfigurations);

            /**
              * @brief Gets the data locations object.
              * @return A pointer to the locations object.
              *
              */
            dataunits::Locations *
            GetLocations();

            /**
             * @brief Gets the kernel object used to compute the covariance matrix.
             * @return A pointer to the kernel object.
             *
             */
            exageostat::kernels::Kernel *
            GetKernel();

            /**
             * @brief Destructor for the data generator object.
             * @details This method frees the memory used by the data generator object.
             *
             */
            virtual ~DataGenerator();

            /**
             * @brief Gets the linear algebra solver object.
             * @return A pointer to the linear algebra solver object.
             *
             */
            static linearAlgebra::LinearAlgebraMethods<T> *GetLinearAlgberaSolver();

        protected:
            /// Pointer to SyntheticDataConfigurations object
            configurations::data_configurations::SyntheticDataConfigurations * mpConfigurations{}; // [in] Used Synthetic Configuration.
            /// Pointer to Locations object
            dataunits::Locations * mpLocations = nullptr; // [out] Used Locations
            /// Pointer to Kernel object
            exageostat::kernels::Kernel * mpKernel = nullptr; // [out] Used Kernel
            /// Pointer to LinearAlgebraMethods object
            static linearAlgebra::LinearAlgebraMethods<T> * mpLinearAlgebraSolver; // [out] Used linear Algebra solver
        };

        /**
         * @brief Instantiates the Data Generator class for float and double types.
         * @tparam T Data Type: float or double
         *
         */
        EXAGEOSTAT_INSTANTIATE_CLASS(DataGenerator)
    }//namespace generators
}//namespace exageostat

#endif //EXAGEOSTAT_CPP_DATAGENERATOR_HPP