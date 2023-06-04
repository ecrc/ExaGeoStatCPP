
/**
 * @file DataGenerator.hpp
 * @brief Contains the declaration of the SyntheticDataConfigurations class.
 * @version 1.0.0
 * @author Sameh Abdulah
 * @date 2023-02-14
**/

#ifndef EXAGEOSTAT_CPP_DATAGENERATOR_HPP
#define EXAGEOSTAT_CPP_DATAGENERATOR_HPP

#include <data-units/Locations.hpp>
#include <configurations/data-generation/concrete/SyntheticDataConfigurations.hpp>
#include <linear-algebra-solvers/LinearAlgebraMethods.hpp>
#include <linear-algebra-solvers/LinearAlgebraFactory.hpp>
#include <memory>

namespace exageostat {
    namespace generators {

        /**
         * @class DataGenerator
         * @brief Abstract base class for generating synthetic or real data.
         * @tparam T The data type of the data generator.
         */
        template<typename T>
        class DataGenerator {

        public:

            /**
             * @brief
             * Generates the data locations.
             * This method generates the X, Y, and Z variables used to define the locations of the data points.
             *
             */
            virtual void
            GenerateLocations() = 0;

            /**
             * @brief
             * Generates the data descriptors.
             * This method generates the descriptors used to define the properties of the data points.
             *
             */
            virtual void
            GenerateDescriptors() = 0;

            virtual void DestoryDescriptors() = 0;
            /**
             * @brief
             * Generates the data observations.
             *
             * This method generates the observations of the data points, which are used to train and test the model.
             *
             * @return void
             */
            virtual void
            GenerateObservations() = 0;

            /**
             * @brief
             * Factory method for creating a data generator object.
             * This method creates a data generator object based on the specified configurations.
             *
             * @param[in] apConfigurations
             * Pointer to the synthetic data configurations.
             *
             * @return std::unique_ptr<DataGenerator>
             * A unique pointer to the created data generator object.
             *
             */
            static std::unique_ptr<DataGenerator>
            CreateGenerator(configurations::data_configurations::SyntheticDataConfigurations *apConfigurations);

            /**
              * @brief
              * Gets the data locations object.
              *
              * @return Locations *
              * A pointer to the locations object.
              */
            dataunits::Locations *
            GetLocations();

            /**
             * @brief
             * Gets the kernel object used to compute the covariance matrix.
             *
             * @return Kernel *
             * A pointer to the kernel object.
             */
            exageostat::kernels::Kernel *
            GetKernel();

            /**
             * @brief
             * Destructor for the data generator object.
             * This method frees the memory used by the data generator object.
             *
             */
            virtual ~DataGenerator() = default;

            static linearAlgebra::LinearAlgebraMethods<T> * GetLinearAlgberaSolver();

        protected:
            /// Used Synthetic Configuration.
            configurations::data_configurations::SyntheticDataConfigurations *mpConfigurations{}; // Pointer to SyntheticDataConfigurations object
            /// Used Locations
            dataunits::Locations * mpLocations = nullptr; // Pointer to Locations object
            /// Used Kernel
            exageostat::kernels::Kernel * mpKernel = nullptr;
            /// Used linear Algebra solver
            static linearAlgebra::LinearAlgebraMethods<T> *mpLinearAlgebraSolver;
        };

        /**
         * @brief Instantiates the LinearAlgebraMethods class for float and double types.
         */
        EXAGEOSTAT_INSTANTIATE_CLASS(DataGenerator)
    }//namespace generators
}//namespace exageostat

#endif //EXAGEOSTAT_CPP_DATAGENERATOR_HPP