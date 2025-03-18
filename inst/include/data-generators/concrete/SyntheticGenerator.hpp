
// Copyright (c) 2017-2024 King Abdullah University of Science and Technology,
// All rights reserved.
// ExaGeoStat is a software package, provided by King Abdullah University of Science and Technology (KAUST).

/**
 * @file SyntheticGenerator.hpp
 * @brief A class for generating synthetic data.
 * @version 1.1.0
 * @author Mahmoud ElKarargy
 * @author Sameh Abdulah
 * @date 2023-02-14
**/

#ifndef EXAGEOSTAT_CPP_SYNTHETICGENERATOR_HPP
#define EXAGEOSTAT_CPP_SYNTHETICGENERATOR_HPP

#include <data-generators/DataGenerator.hpp>

namespace exageostat::generators::synthetic {

    /**
     * @class SyntheticGenerator
     * @brief A class for generating synthetic data.
     * @tparam T Data Type: float or double
     * @details This class generates synthetic data for use in testing machine learning models.
     *
     */
    template<typename T>
    class SyntheticGenerator : public DataGenerator<T> {

    public:

        /**
         * @brief Creates the data by synthetically generating it.
         * @copydoc DataGenerator::CreateData()
         *
         */
        std::unique_ptr<ExaGeoStatData<T>>
        CreateData(configurations::Configurations &aConfigurations, exageostat::kernels::Kernel<T> &aKernel) override;

         /**
         * @brief Abstract method for synthetic data generation based on provided configurations and kernel.
         * @param[in] aConfigurations Reference to the configurations object that contains parameters for generating data.
         * @param[in] aKernel Reference to the kernel object that defines the operations to be applied while generating the data.
         * @return A unique pointer to the generated ExaGeoStatData object.
         *
         */
        virtual std::unique_ptr<ExaGeoStatData<T>>
        CreateSyntheticData(configurations::Configurations &aConfigurations, exageostat::kernels::Kernel<T> &aKernel)  = 0;

        /**
         * @brief Factory method for creating a synthetic data generator instance.
         * This method dynamically determines the type of synthetic generator to instantiate based on compile-time conditions.
         * @return A unique pointer to a synthetic data generator instance configured as per the specified runtime conditions.
         *
         */
        static std::unique_ptr<SyntheticGenerator<T>> CreateSyntheticGenerator();

        /**
         * @brief Releases the singleton instance of the currently active synthetic data generator.
         * This method ensures proper deallocation of the singleton instance of the synthetic data generator,
         * depending on the selected runtime.
         *
         */
        static void ReleaseSyntheticGenerator();

    };

    /**
     * @brief Instantiates the Synthetic Data Generator class for float and double types.
     * @tparam T Data Type: float or double
     *
     */
    EXAGEOSTAT_INSTANTIATE_CLASS(SyntheticGenerator)
} // namespace exageostat

#endif //EXAGEOSTAT_CPP_SYNTHETICGENERATOR_HPP