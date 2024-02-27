
// Copyright (c) 2017-2023 King Abdullah University of Science and Technology,
// All rights reserved.
// ExaGeoStat is a software package, provided by King Abdullah University of Science and Technology (KAUST).

/**
 * @file ChameleonDescriptor.cpp
 * @brief Example file for the Chameleon to Hicma Converter.
 * @version 1.0.1
 * @author Mahmoud ElKarargy
 * @date 2024-02-14
**/

#include <kernels/Kernel.hpp>
#include <data-units/DescriptorData.hpp>
#include <hardware/ExaGeoStatHardware.hpp>
#include <configurations/Configurations.hpp>
#include <utilities/Logger.hpp>

using namespace std;

using namespace exageostat::common;
using namespace exageostat::kernels;
using namespace exageostat::plugins;
using namespace exageostat::dataunits;

int main(int argc, char **argv) {

    LOGGER("** Example of Chameleon To Hicma Converter **")

    // Initialize Configuration
    Configurations configuration;
    configuration.InitializeArguments(argc, argv);

    // Initialize Hardware and Data
    auto hardware = ExaGeoStatHardware(configuration.GetComputation(), configuration.GetCoresNumber(),
                                       configuration.GetGPUsNumbers());
    unique_ptr<DescriptorData<double>> data = make_unique<DescriptorData<double>>();
    Kernel<double> *kernel = PluginRegistry<Kernel<double>>::Create(
            configuration.GetKernelName(), configuration.GetTimeSlot());

    // Get arguments for Descriptors Initialization
    int kernel_variables_number = kernel->GetVariablesNumber();
    int config_problem_size = configuration.GetProblemSize();
    int config_full_problem_size = config_problem_size * kernel_variables_number;
    int config_dts = configuration.GetDenseTileSize();
    int config_p_grid = configuration.GetPGrid();
    int config_q_grid = configuration.GetQGrid();
    bool config_is_OOC = configuration.GetIsOOC();

    // Randomly Initialized Matrix of Data
    std::vector<double> matrix(config_problem_size);
    for (int i = 0; i < config_problem_size; ++i) {
        matrix[i] = i;
    }

    // Set Data Descriptor
    data->SetDescriptor(CHAMELEON_DESCRIPTOR, DESCRIPTOR_C, config_is_OOC, matrix.data(), EXAGEOSTAT_REAL_DOUBLE,
                        config_dts,
                        config_dts,
                        config_dts * config_dts, config_full_problem_size, config_full_problem_size, 0, 0,
                        config_full_problem_size,
                        config_full_problem_size, config_p_grid, config_q_grid);

    auto *CHAM_descriptorC = data->GetDescriptor(CHAMELEON_DESCRIPTOR, DESCRIPTOR_C).chameleon_desc;

    //Print Descriptor Attributes Before and After conversion
    LOGGER("** Chameleon Descriptor, Before Conversion:")
    LOGGER(" Problem Size: " << CHAM_descriptorC->m)
    LOGGER(" Dense Tile Size: " << CHAM_descriptorC->mb)

    LOGGER(" Entire Number of Rows :" << CHAM_descriptorC->lm)
    LOGGER(" Entire Number of Columns :" << CHAM_descriptorC->ln)

    LOGGER(" Number of Sub-matrix Tile Rows: " << CHAM_descriptorC->mt)
    LOGGER(" Number of Sub-matrix Tile Columns: " << CHAM_descriptorC->nt)

    LOGGER(" Number of Rows of 2D distribution grid: " << CHAM_descriptorC->p)
    LOGGER(" Number of Rows of 2D distribution grid: " << CHAM_descriptorC->q)

    LOGGER(" Is Matrix Not Fit in Memory: " << CHAM_descriptorC->ooc)
    LOGGER(" Size including Padding: " << CHAM_descriptorC->bsiz)

    // Print Data Matrix of Descriptor
    LOGGER("** Data in Matrix of Chameleon descriptor:")
    LOGGER_2("",0)
    auto *cham_mat = (double *) CHAM_descriptorC->mat;
    for (int i = 0; i < config_problem_size; ++i) {
        LOGGER_PRECISION_1(" " << cham_mat[i],0)
    }

    auto *HICMA_descriptorC = data->ConvertChameleonToHicma(CHAM_descriptorC);

    LOGGER("\n\t\t ** Hicma Descriptor, After Conversion:")
    LOGGER(" Problem Size: " << HICMA_descriptorC->m)
    LOGGER(" Dense Tile Size: " << HICMA_descriptorC->mb)

    LOGGER(" Entire Number of Rows :" << HICMA_descriptorC->lm)
    LOGGER(" Entire Number of Columns :" << HICMA_descriptorC->ln)

    LOGGER(" Number of Sub-matrix Tile Rows: " << HICMA_descriptorC->mt)
    LOGGER(" Number of Sub-matrix Tile Columns: " << HICMA_descriptorC->nt)

    LOGGER(" Number of Rows of 2D distribution grid: " << HICMA_descriptorC->p)
    LOGGER(" Number of Rows of 2D distribution grid: " << HICMA_descriptorC->q)

    LOGGER(" Is Matrix Not Fit in Memory: " << HICMA_descriptorC->ooc)
    LOGGER(" Size including Padding: " << HICMA_descriptorC->bsiz)

    // Print Data Matrix of Descriptor
    LOGGER("** Data in Matrix of Hicma descriptor:")
    auto *hicma_mat = (double *) HICMA_descriptorC->mat;
    LOGGER_2("",0)
    for (int i = 0; i < config_problem_size; ++i) {
        LOGGER_PRECISION_1(" " << hicma_mat[i],0)
    }

    delete kernel;
    delete HICMA_descriptorC;
}