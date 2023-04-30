
// Copyright (c) 2017-2023 King Abdullah University of Science and Technology,
// Copyright (C) 2023 by Brightskies inc,
// All rights reserved.
// ExaGeoStat is a software package, provided by King Abdullah University of Science and Technology (KAUST).

/**
 * @file ChameleonAllocateDescriptors.cpp
 * @brief Sets up the Chameleon descriptors needed for the dense matrix computations in ExaGeoStat.
 * @version 1.0.0
 * @author Sameh Abdulah
 * @date 2023-03-20
**/

#include <linear-algebra-solvers/concrete/dense/ChameleonImplementationDense.hpp>
// Include Chameleon libraries
extern "C" {
#include <chameleon/struct.h>
#include <chameleon.h>
#include <control/descriptor.h>
#include <control/context.h>
}

// Use the following namespaces for convenience
using namespace exageostat::linearAlgebra::dense;
using namespace exageostat::common;
using namespace std;

// Define a method to set up the Chameleon descriptors
template<typename T>
void ChameleonImplementationDense<T>::InitiateDescriptors() {

    // Initialize the Chameleon context
    this->ExaGeoStatInitContext(this->mpConfigurations->GetCoresNumber(), this->mpConfigurations->GetGPUsNumber());

    // Declare variables for Chameleon descriptors
    vector<void *> &pDescriptorC = this->mpConfigurations->GetDescriptorC();
    vector<void *> &pDescriptorZ = this->mpConfigurations->GetDescriptorZ();
    auto pDescriptorZcpy = (CHAM_desc_t **) &this->mpConfigurations->GetDescriptorZcpy();
    vector<void *> &pDescriptorProduct = this->mpConfigurations->GetDescriptorProduct();
    auto pDescriptorDeterminant = (CHAM_desc_t **) &this->mpConfigurations->GetDescriptorDeterminant();

    // Get the problem size and other configuration parameters
    int vectorSize;
    int N = this->mpConfigurations->GetProblemSize() * this->mpConfigurations->GetP();
    int dts = this->mpConfigurations->GetDenseTileSize();
    int pGrid = this->mpConfigurations->GetPGrid();
    int qGrid = this->mpConfigurations->GetQGrid();
    bool isOOC = this->mpConfigurations->GetIsOOC();

    // For distributed system and should be removed
    T *Zcpy = (T *) malloc(N * sizeof(T));
    T dotProductValue;

    // Create a Chameleon sequence
    RUNTIME_sequence_t *pSequence;
    CHAMELEON_Sequence_Create(&pSequence);

    // Set the floating point precision based on the template type
    FloatPoint floatPoint;
    if (sizeof(T) == SIZE_OF_FLOAT) {
        floatPoint = EXAGEOSTAT_REAL_FLOAT;
        vectorSize = 1;
    } else {
        floatPoint = EXAGEOSTAT_REAL_DOUBLE;
        vectorSize = 3;
    }

    // Create the Chameleon descriptors based on the configuration parameters
    // Depending on the passed Precession, the descriptor will resize for value 1 or 3.
    pDescriptorC.resize(vectorSize + 1, nullptr);
    pDescriptorZ.resize(vectorSize, nullptr);
    pDescriptorProduct.resize(vectorSize, nullptr);

    auto **CHAM_descriptorC = (CHAM_desc_t **) &pDescriptorC[0];
    EXAGEOSTAT_ALLOCATE_DENSE_MATRIX_TILE(CHAM_descriptorC, isOOC, nullptr, (cham_flttype_t) floatPoint, dts, dts,
                                          dts * dts, N, N, 0, 0, N, N, pGrid, qGrid)

    if (vectorSize > 1) {

        auto **CHAM_descsubC11 = (CHAM_desc_t **) &pDescriptorC[1];
        auto **CHAM_descsubC12 = (CHAM_desc_t **) &pDescriptorC[2];
        auto **CHAM_descsubC22 = (CHAM_desc_t **) &pDescriptorC[3];

        *CHAM_descsubC11 = chameleon_desc_submatrix(*CHAM_descriptorC, 0, 0, (*CHAM_descriptorC)->m / 2,
                                                    (*CHAM_descriptorC)->n / 2);
        *CHAM_descsubC12 = chameleon_desc_submatrix(*CHAM_descriptorC, (*CHAM_descriptorC)->m / 2, 0,
                                                    (*CHAM_descriptorC)->m / 2, (*CHAM_descriptorC)->n / 2);
        *CHAM_descsubC22 = chameleon_desc_submatrix(*CHAM_descriptorC, (*CHAM_descriptorC)->m / 2,
                                                    (*CHAM_descriptorC)->n / 2,
                                                    (*CHAM_descriptorC)->m / 2, (*CHAM_descriptorC)->n / 2);
    }
    EXAGEOSTAT_ALLOCATE_DENSE_MATRIX_TILE(((CHAM_desc_t **) &pDescriptorZ[0]), isOOC, nullptr,
                                          (cham_flttype_t) floatPoint, dts, dts, dts * dts, N, 1, 0, 0, N, 1, pGrid,
                                          qGrid)
    EXAGEOSTAT_ALLOCATE_DENSE_MATRIX_TILE(pDescriptorZcpy, isOOC, Zcpy, (cham_flttype_t) floatPoint, dts, dts,
                                          dts * dts, N, 1, 0, 0, N, 1, pGrid, qGrid)
    EXAGEOSTAT_ALLOCATE_DENSE_MATRIX_TILE(pDescriptorDeterminant, isOOC, &dotProductValue, (cham_flttype_t) floatPoint,
                                          dts, dts, dts * dts, 1, 1, 0, 0, 1, 1, pGrid, qGrid)

    for (int idx = 1; idx < pDescriptorZ.size(); idx++) {
        auto **CHAM_descriptorZ_ = (CHAM_desc_t **) &pDescriptorZ[idx];
        EXAGEOSTAT_ALLOCATE_DENSE_MATRIX_TILE(CHAM_descriptorZ_, isOOC, nullptr, (cham_flttype_t) floatPoint, dts, dts,
                                              dts * dts, N / 2, 1, 0, 0, N / 2, 1, pGrid, qGrid)
    }

    for (auto & idx : pDescriptorProduct) {
        auto **CHAM_descriptorProduct = (CHAM_desc_t **) &idx;
        EXAGEOSTAT_ALLOCATE_DENSE_MATRIX_TILE(CHAM_descriptorProduct, isOOC, &dotProductValue,
                                              (cham_flttype_t) floatPoint, dts, dts, dts * dts, 1, 1, 0, 0, 1, 1, pGrid,
                                              qGrid)
    }
    this->ExaGeoStatFinalizeContext();
    //stop gsl error handler
    gsl_set_error_handler_off();
}

template<typename T>
void ChameleonImplementationDense<T>::ExaGeoStatInitContext(const int &apCoresNumber, const int &apGPUs) {

    CHAM_context_t *chameleonContext;
    chameleonContext = chameleon_context_self();
    if (chameleonContext != nullptr) {
        printf("Another instance of Chameleon is already running...!");
    } else {
        CHAMELEON_user_tag_size(31, 26);
        CHAMELEON_Init(apCoresNumber, apGPUs)
    }
}

template<typename T>
void ChameleonImplementationDense<T>::ExaGeoStatFinalizeContext() {

    CHAM_context_t *chameleonContext;
    chameleonContext = chameleon_context_self();
    if (chameleonContext == nullptr) {
        printf("No active instance oh Chameleon...please use ExaGeoStatInitContext() function to initiate a new instance!\n");
    } else
        CHAMELEON_Finalize();
}


//#include <kernels/Kernel.hpp>
#include <data-generators/DataGenerator.hpp>

void *EXAGEOSTAT_data_getaddr(const CHAM_desc_t *A, int type, int m, int n) {
    int64_t mm = m + (A->i / A->mb);
    int64_t nn = n + (A->j / A->nb);

    starpu_data_handle_t *ptrtile = static_cast<starpu_data_handle_t *>(A->schedopt);
    ptrtile += ((int64_t) A->lmt) * nn + mm;

    if (*ptrtile == NULL) {
        int home_node = -1;
        void *user_ptr = NULL;
        int myrank = A->myrank;
        int owner = A->get_rankof(A, m, n);
        int64_t eltsze = CHAMELEON_Element_Size(type);
        int tempmm = (mm == A->lmt - 1) ? (A->lm - mm * A->mb) : A->mb;
        int tempnn = (nn == A->lnt - 1) ? (A->ln - nn * A->nb) : A->nb;

        if (myrank == owner) {
            user_ptr = A->get_blkaddr(A, m, n);
            if (user_ptr != NULL) {
                home_node = STARPU_MAIN_RAM;
            }
        }

        starpu_matrix_data_register(ptrtile, home_node, (uintptr_t) user_ptr,
                                    BLKLDD(A, m),
                                    tempmm, tempnn, eltsze);

#ifdef HAVE_STARPU_DATA_SET_COORDINATES
        starpu_data_set_coordinates(*ptrtile, 2, m, n);
#endif

#if defined(CHAMELEON_USE_MPI)
        {
            int64_t block_ind = A->lmt * nn + mm;
            starpu_mpi_data_register(*ptrtile, (A->id << tag_sep) | (block_ind), owner);
        }
#endif /* defined(CHAMELEON_USE_MPI) */
    }
    return *ptrtile;
}

#define EXAGEOSTAT_RTBLKADDR(desc, type, m, n) ( (starpu_data_handle_t)EXAGEOSTAT_data_getaddr( desc, type, m, n ) )


template<typename T>
void ChameleonImplementationDense<T>::testKernelfornow() {
//    auto kernel = new exageostat::kernels::UnivariateMaternStationary();
//
//    auto **CHAM_descriptorC = (CHAM_desc_t **) &this->mpConfigurations->GetDescriptorC()[0];
//
//    // Create a unique pointer to a DataGenerator object
//    unique_ptr<exageostat::generators::DataGenerator> d1;
//
//    // Create the DataGenerator object
//    d1 = d1->CreateGenerator(
//            dynamic_cast<configurations::data_configurations::SyntheticDataConfigurations *>(this->mpConfigurations));
//
//    // Initialize the locations of the generated data
//    d1->GenerateLocations();
//    dataunits::Locations * l1 = d1->GetLocations();
//
//    unique_ptr<exageostat::generators::DataGenerator> d2;
//    // Create the DataGenerator object
//    d2 = d2->CreateGenerator(
//            dynamic_cast<configurations::data_configurations::SyntheticDataConfigurations *>(this->mpConfigurations));
//
//    // Initialize the locations of the generated data
//    d2->GenerateLocations();
//    dataunits::Locations * l2 = d2->GetLocations();
//    double * initial_theta = (double* ) malloc(3 * sizeof(double));
//
//    initial_theta[0] = 1.0;
//    initial_theta[1] = 0.1;
//    initial_theta[2] = 0.5;
//
//    auto * A = (double* ) EXAGEOSTAT_RTBLKADDR((*CHAM_descriptorC), ChamRealDouble, 0, 0);
//    kernel->GenerateCovarianceMatrix(A, 5, 5, 0, 0, l1, l1, nullptr, initial_theta, 0);
    /*
     * 1.000000 0.108306 0.004480 0.012114 0.015264
     * 0.108306 1.000000 0.030836 0.017798 0.062895
     * 0.004480 0.030836 1.000000 0.001011 0.007046
     * 0.012114 0.017798 0.001011 1.000000 0.123679
     * 0.015264 0.062895 0.007046 0.123679 1.000000
     */
//    kernel->GenerateCovarianceMatrix(A, 4, 5, 5, 0, l1, l1, nullptr, initial_theta, 0);
    /*
     * 0.002145 0.004299 0.000404 0.171258 0.055331
     * 0.000468 0.001824 0.000498 0.021786 0.028073
     * 0.000221 0.001990 0.009302 0.000805 0.005396
     * 0.000061 0.000440 0.000713 0.000907 0.003317
     */
//    kernel->GenerateCovarianceMatrix(A, 4, 4, 5, 5, l1, l1, nullptr, initial_theta, 0);
    /*
     * 1.000000 0.085375 0.000986 0.002264
     * 0.085375 1.000000 0.005156 0.023215
     * 0.000986 0.005156 1.000000 0.053542
     * 0.002264 0.023215 0.053542 1.000000
     */
}
