
// Copyright (c) 2017-2023 King Abdullah University of Science and Technology,
// All rights reserved.
// ExaGeoStat is a software package, provided by King Abdullah University of Science and Technology (KAUST).

/**
 * @file ChameleonImplementationDST.cpp
 * @brief Diagonal Super Tile implementation of linear algebra methods.
 * @version 1.0.0
 * @author Mahmoud ElKarargy
 * @author Sameh Abdulah
 * @date 2023-03-20
**/

#include <linear-algebra-solvers/concrete/chameleon/diagonal-super-tile/ChameleonImplementationDST.hpp>

using namespace std;

using namespace exageostat::linearAlgebra::diagonalSuperTile;
using namespace exageostat::common;

template<typename T>
void ChameleonImplementationDST<T>::ExaGeoStatPotrfTile(const UpperLower &aUpperLower, void *apA, int aBand, void *apCD,
                                                        void *apCrk, const int &aMaxRank, const int &aAcc) {

    CHAM_context_t *chameleon_context;
    RUNTIME_sequence_t *sequence = nullptr;
    RUNTIME_request_t request = RUNTIME_REQUEST_INITIALIZER;

    chameleon_context = chameleon_context_self();
    if (chameleon_context == nullptr) {
        throw std::runtime_error("CHAMELEON_dpotrf_diag_Tile() Failed, Hardware not Initialized.");
    }
    chameleon_sequence_create(chameleon_context, &sequence);
    ExaGeoStatPotrfDiagonalTileAsync(aUpperLower, apA, aBand, sequence, &request);
    CHAMELEON_Desc_Flush((CHAM_desc_t *) apA, sequence);
    chameleon_sequence_wait(chameleon_context, sequence);
    chameleon_sequence_destroy(chameleon_context, sequence);
}

template<typename T>
int ChameleonImplementationDST<T>::ExaGeoStatPotrfDiagonalTileAsync(const common::UpperLower &aUpperLower, void *apA,
                                                                    int aBand, void *apSequence, void *apRequest) {

    CHAM_context_t *chameleon_context;
    chameleon_context = chameleon_context_self();
    if (chameleon_context == nullptr) {
        chameleon_fatal_error("CHAMELEON_dpotrf_Tile_Async", "CHAM not initialized");
        return CHAMELEON_ERR_NOT_INITIALIZED;
    }
    if (apSequence == nullptr) {
        chameleon_fatal_error("CHAMELEON_dpotrf_Tile_Async", "NULL sequence");
        return CHAMELEON_ERR_UNALLOCATED;
    }
    if (apRequest == nullptr) {
        chameleon_fatal_error("CHAMELEON_dpotrf_Tile_Async", "NULL request");
        return CHAMELEON_ERR_UNALLOCATED;
    }
    /* Check sequence status */
    if (((RUNTIME_sequence_t *) apSequence)->status == CHAMELEON_SUCCESS) {
        ((RUNTIME_request_t *) apRequest)->status = CHAMELEON_SUCCESS;
    } else {
        return chameleon_request_fail((RUNTIME_sequence_t *) apSequence, (RUNTIME_request_t *) apRequest,
                                      CHAMELEON_ERR_ILLEGAL_VALUE);
    }

    /* Check descriptors for correctness */
    if (chameleon_desc_check((CHAM_desc_t *) apA) != CHAMELEON_SUCCESS) {
        chameleon_error("CHAMELEON_dpotrf_Tile_Async", "invalid descriptor");
        return chameleon_request_fail((RUNTIME_sequence_t *) apSequence, (RUNTIME_request_t *) apRequest,
                                      CHAMELEON_ERR_ILLEGAL_VALUE);
    }
    /* Check input arguments */
    if (((CHAM_desc_t *) apA)->nb != ((CHAM_desc_t *) apA)->mb) {
        chameleon_error("CHAMELEON_dpotrf_Tile_Async", "only square tiles supported");
        return chameleon_request_fail((RUNTIME_sequence_t *) apSequence, (RUNTIME_request_t *) apRequest,
                                      CHAMELEON_ERR_ILLEGAL_VALUE);
    }
    if (aUpperLower != EXAGEOSTAT_UPPER && aUpperLower != EXAGEOSTAT_LOWER) {
        chameleon_error("CHAMELEON_dpotrf_Tile_Async", "illegal value of uplo");
        return chameleon_request_fail((RUNTIME_sequence_t *) apSequence, (RUNTIME_request_t *) apRequest, -1);
    }
    /* Quick return */
    ExaGeoStatParallelPotrfDiagonal(aUpperLower, apA, aBand, apSequence, apRequest);
    return CHAMELEON_SUCCESS;
}


template<typename T>
void ChameleonImplementationDST<T>::ExaGeoStatParallelPotrfDiagonal(const common::UpperLower &aUpperLower, void *apA,
                                                                    int aBand, void *apSequence, void *apRequest) {
    CHAM_context_t *chameleon_context;
    RUNTIME_option_t options;

    int k, m, n;
    int ldak, ldam, ldan;
    int tempkm, tempmm, tempnn;
    size_t ws_host = 0;

    auto zone = (T) 1.0;
    auto mzone = (T) -1.0;

    chameleon_context = chameleon_context_self();
    if (((RUNTIME_sequence_t *) apSequence)->status != CHAMELEON_SUCCESS)
        return;

    RUNTIME_options_init(&options, chameleon_context, ((RUNTIME_sequence_t *) apSequence),
                         ((RUNTIME_request_t *) apRequest));
    RUNTIME_options_ws_alloc(&options, 0, ws_host);
    auto A = (CHAM_desc_t *) apA;

    /*
     *  ChamLower
     */
    if (aUpperLower == EXAGEOSTAT_LOWER) {
        for (k = 0; k < A->mt; k++) {
            RUNTIME_iteration_push(chameleon_context, k);

            tempkm = k == A->mt - 1 ? A->m - k * A->mb : A->mb;
            ldak = BLKLDD(A, k);

            options.priority = 2 * A->mt - 2 * k;
            INSERT_TASK_dpotrf(&options, ChamLower, tempkm, A->mb, A, k, k, A->nb * k);

            for (m = k + 1; m < A->mt && m < k + aBand; m++) {
                tempmm = m == A->mt - 1 ? A->m - m * A->mb : A->mb;
                ldam = BLKLDD(A, m);

                options.priority = 2 * A->mt - 2 * k - m;
                INSERT_TASK_dtrsm(&options, ChamRight, ChamLower, ChamTrans, ChamNonUnit, tempmm, A->mb, A->mb, zone, A,
                                  k, k, A, m, k);
            }
            RUNTIME_data_flush((RUNTIME_sequence_t *) apSequence, A, k, k);

            for (n = k + 1; n < A->nt && n < k + aBand; n++) {
                tempnn = n == A->nt - 1 ? A->n - n * A->nb : A->nb;
                ldan = BLKLDD(A, n);

                options.priority = 2 * A->mt - 2 * k - n;
                INSERT_TASK_dsyrk(&options, ChamLower, ChamNoTrans, tempnn, A->nb, A->mb, -1.0, A, n, k, 1.0, A, n, n);

                for (m = n + 1; m < A->mt && m < n + aBand; m++) {
                    tempmm = m == A->mt - 1 ? A->m - m * A->mb : A->mb;
                    ldam = BLKLDD(A, m);

                    options.priority = 2 * A->mt - 2 * k - n - m;
                    INSERT_TASK_dgemm(&options, ChamNoTrans, ChamTrans, tempmm, tempnn, A->mb, A->mb, mzone, A, m, k, A,
                                      n, k, zone, A, m, n);
                }
                RUNTIME_data_flush((RUNTIME_sequence_t *) apSequence, A, n, k);
            }
            RUNTIME_iteration_pop(chameleon_context);
        }
    }
        /*
         *  ChamUpper
         */
    else {
        for (k = 0; k < A->nt; k++) {
            RUNTIME_iteration_push(chameleon_context, k);

            tempkm = k == A->nt - 1 ? A->n - k * A->nb : A->nb;
            ldak = BLKLDD(A, k);

            options.priority = 2 * A->nt - 2 * k;
            INSERT_TASK_dpotrf(&options, ChamUpper, tempkm, A->mb, A, k, k, A->nb * k);

            for (n = k + 1; n < A->nt; n++) {
                tempnn = n == A->nt - 1 ? A->n - n * A->nb : A->nb;

                options.priority = 2 * A->nt - 2 * k - n;
                INSERT_TASK_dtrsm(&options, ChamLeft, ChamUpper, ChamTrans, ChamNonUnit, A->mb, tempnn, A->mb, zone, A,
                                  k, k, A, k, n);
            }
            RUNTIME_data_flush((RUNTIME_sequence_t *) apSequence, A, k, k);

            for (m = k + 1; m < A->mt; m++) {
                tempmm = m == A->mt - 1 ? A->m - m * A->mb : A->mb;
                ldam = BLKLDD(A, m);

                options.priority = 2 * A->nt - 2 * k - m;
                INSERT_TASK_dsyrk(&options, ChamUpper, ChamTrans, tempmm, A->mb, A->mb, -1.0, A, k, m, 1.0, A, m, m);

                for (n = m + 1; n < A->nt; n++) {
                    tempnn = n == A->nt - 1 ? A->n - n * A->nb : A->nb;

                    options.priority = 2 * A->nt - 2 * k - n - m;
                    INSERT_TASK_dgemm(&options, ChamTrans, ChamNoTrans, tempmm, tempnn, A->mb, A->mb, mzone, A, k, m, A,
                                      k, n, zone, A, m, n);
                }
                RUNTIME_data_flush((RUNTIME_sequence_t *) apSequence, A, k, m);
            }
            RUNTIME_iteration_pop(chameleon_context);
        }
    }

    RUNTIME_options_ws_free(&options);
    RUNTIME_options_finalize(&options, chameleon_context);
}