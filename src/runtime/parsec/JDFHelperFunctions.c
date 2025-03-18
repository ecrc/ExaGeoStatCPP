
// Copyright (c) 2017-2024 King Abdullah University of Science and Technology,
// All rights reserved.
// ExaGeoStat is a software package, provided by King Abdullah University of Science and Technology (KAUST).

/**
 * @file JDFHelperFunctions.c
 * @brief Implementation of JDF helper functions.
 * @version 2.0.0
 * @author Mahmoud ElKarargy
 * @author Sameh Abdulah
 * @author Qinglei Cao
 * @date 2024-10-20
**/

#include <runtime/parsec/JDFHelperFunctions.h>
#include <gb24/gb24.h>

#ifdef USE_CUDA
#include <runtime/parsec/GPUHelperFunctions.h>
#endif


int CalculateSingleIndex(int aN, int aM) {
    return aN * (aN + 1) / 2 + aM;
}

double SumDoubleData(double *apData, int aColumn, int aRow) {
    double sum = 0.0;
    for (int j = 0; j < aRow; j++) {
        for (int i = 0; i < aColumn; i++) {
            sum += apData[j * aColumn + i];
        }
    }
    return sum;
}

complex double SumComplexData(complex double *apData, int aColumn, int aRow) {
    complex double sum = 0.0;
    for (int j = 0; j < aRow; j++) {
        for (int i = 0; i < aColumn; i++) {
            sum += apData[j * aColumn + i];
        }
    }
    return sum;
}

void ForwardSHTHelper(double *apFlm, complex double *apF_data, int aFDataM, int aFDataN,
                      complex double *apEt1, int aEt1M, complex double *apEt2, int aEt2M,
                      complex double *apEp, int aEpM, int aEpN, complex double *apSlmn,
                      int aSlmnM, int aSlmnN, complex double *apIe, int aIeM, int aIeN,
                      complex double *apIo,  int aIoM, int aIoN, complex double *apP,
                      int aPM, int aPN, complex double *apD, complex double *apGmtheta_r,
                      complex double *apFmnm, complex double *apTmp1,
                      complex double *apTmp2, int aL){

    complex double alpha_complex, beta_complex;
    double alpha_double, beta_double;

    assert(aFDataN == aEpM);
    alpha_complex = (complex double) 1.0;
    beta_complex = (complex double) 0.0;

    int gmtheta_r_M = aFDataM;
    int gmtheta_r_N = aEpN;
    int gmtheta_r_K = aFDataN;

    cblas_zgemm(CblasColMajor, CblasNoTrans, CblasNoTrans,
                gmtheta_r_M, gmtheta_r_N, gmtheta_r_K,
                &alpha_complex, apF_data, gmtheta_r_M,
                apEp, gmtheta_r_K, &beta_complex,
                apGmtheta_r, gmtheta_r_M);

    int fmnm_M = aEt1M;
    int fmnm_N = gmtheta_r_N;
    int fmnm_K = gmtheta_r_M;

    cblas_zgemm(CblasColMajor, CblasNoTrans, CblasNoTrans,
                fmnm_M, fmnm_N, fmnm_K,
                &alpha_complex,apEt1, fmnm_M,
                apGmtheta_r, fmnm_K,
                &beta_complex, apFmnm, fmnm_M);

    int tmp1_M = aEt2M;
    int tmp1_N = aPN;
    int tmp1_K = aPM;
    cblas_zgemm(CblasColMajor, CblasNoTrans, CblasNoTrans,
                tmp1_M, tmp1_N, tmp1_K,
                &alpha_complex, apEt2, tmp1_M,
                apP, tmp1_K,
                &beta_complex, apTmp1, tmp1_M);

    assert(tmp1_N == gmtheta_r_M);
    int tmp2_M = tmp1_M;
    int tmp2_N = gmtheta_r_N;
    int tmp2_K = gmtheta_r_M;

    cblas_zgemm(CblasColMajor, CblasNoTrans, CblasNoTrans,
                tmp2_M, tmp2_N, tmp2_K,
                &alpha_complex, apTmp1, tmp2_M,
                apGmtheta_r, tmp2_K,
                &beta_complex, apTmp2, tmp2_M);

    assert(fmnm_M == tmp2_M);
    fmnm_K = tmp2_N;
    beta_complex = (complex double) 1.0;
    cblas_zgemm(CblasColMajor, CblasNoTrans, CblasNoTrans,
                fmnm_M, fmnm_N, fmnm_K,
                &alpha_complex, apTmp2, fmnm_M,
                apD, fmnm_K,
                &beta_complex, apFmnm, fmnm_M);

    assert(aSlmnN == aIeM);
    assert(aIeN == fmnm_M);

    int flmn_matrix_M = aL;

    complex double *pFlmn_matrix = apTmp1;
    complex double *pFmnm_tmp;
    complex double *pSlmn_tmp;
    complex double *multipy_tmp = apTmp2 + fmnm_M + aSlmnN;

    for (int m = 0; m < aL; m++) {
        pFmnm_tmp = apFmnm + (aL + m - 1) * fmnm_M;
        if (0 == m % 2) {
            for (int n = m; n < aL; n++) {
                pSlmn_tmp = apSlmn + CalculateSingleIndex(n, m);

                alpha_complex = (complex double) 1.0;
                beta_complex = (complex double) 0.0;
                cblas_zgemv(CblasColMajor, CblasNoTrans,
                            aIeM, aIeN,
                            &alpha_complex, apIe, aIeM,
                            pFmnm_tmp, 1,
                            &beta_complex, multipy_tmp, 1);
                cblas_zdotu_sub(aSlmnN, pSlmn_tmp, aSlmnM, multipy_tmp, 1, &pFlmn_matrix[m * flmn_matrix_M + n]);
            }
        } else {
            for (int n = m; n < aL; n++) {
                pSlmn_tmp = apSlmn + CalculateSingleIndex(n, m);

                alpha_complex = (complex double) 1.0;
                beta_complex = (complex double) 0.0;
                cblas_zgemv(CblasColMajor, CblasNoTrans,
                            aIoM, aIoN,
                            &alpha_complex, apIo, aIoM,
                            pFmnm_tmp, 1,
                            &beta_complex, multipy_tmp, 1);

                cblas_zdotu_sub(aSlmnN, pSlmn_tmp, aSlmnM, multipy_tmp, 1, &pFlmn_matrix[m * flmn_matrix_M + n]);
            }
        }

    }
    for (int n = 0; n < aL; n++) {
        for (int m = 0; m <= n; m++) {
            apFlm[n * n + n + m] = creal(pFlmn_matrix[m * flmn_matrix_M + n]);
            if (m != 0) {
                apFlm[n * n + n - m] = cimag(pFlmn_matrix[m * flmn_matrix_M + n]);
            }
        }
    }
}

void InverseSHTHelper(double *apFlm, double *apFspatial, double *apZlm, double *apSC, double *apSmt, int aL) {

    int index_Zlm, index_flm;
    int Smt_M = aL + 1;
    int Smt_N = 2 * aL - 1;

    memset(apSmt, 0, Smt_M * Smt_N * sizeof(double));

    for (int m = -(aL - 1); m < aL; m++) {
        for (int n = abs(m); n < aL; n++) {
            index_Zlm = CalculateSingleIndex(n, abs(m));
            index_flm = n * n + n + m;
            cblas_daxpy(Smt_M, apFlm[index_flm], apZlm + index_Zlm * aL + 1, 1, apSmt + (m + aL - 1) * Smt_M, 1);
        }
    }

    int f_spatial_M = Smt_M;
    int f_spatial_N = 2 * aL;
    int f_spatial_K = Smt_N;

    cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans,
                f_spatial_M, f_spatial_N, f_spatial_K,
                (double) 1.0, apSmt, f_spatial_M,
                apSC, f_spatial_K,
                (double) 0.0, apFspatial, f_spatial_M);

}

double ParsecMatrixSumCore(double *apA, int aMb, int aNb, int aLda) {
    double res = 0.0;
    int indicator = 0;
    for (int j = 0; j < aNb; j++) {
        for (int i = 0; i < aMb; i++) {
            if (isnan(apA[j * aLda + i])) {
                indicator = 1;
            }
            res += apA[j * aLda + i];
        }
    }
    return res;
}

#ifdef USE_CUDA

void ForwardSHTGPUCore(double *apFlm, cuDoubleComplex *apF_data, cuDoubleComplex *apEt1,
                       cuDoubleComplex *apEt2, cuDoubleComplex *apEp, cuDoubleComplex *apSlmn,
                       cuDoubleComplex *apIe, cuDoubleComplex *apIo, cuDoubleComplex *apP,
                       cuDoubleComplex *apD, parsec_device_cuda_module_t *apCudaDevice,
                       parsec_gpu_task_t *apGpuTask, parsec_cuda_exec_stream_t *apCudaStream,
                       void* apWorkSpace, int aL) {

    int f_data_M = ((parsec_tiled_matrix_t*)apF_data)->mb;
    int Ep_N = ((parsec_tiled_matrix_t*)apEp)->nb;
    int f_data_N = ((parsec_tiled_matrix_t*)apF_data)->mb;
    int Et1_M = ((parsec_tiled_matrix_t*)apEt1)->mb;
    int Et2_M = ((parsec_tiled_matrix_t*)apEt2)->mb;
    int P_N = ((parsec_tiled_matrix_t*)apP)->nb;
    int P_M = ((parsec_tiled_matrix_t*)apP)->mb;

    int Slmn_N = ((parsec_tiled_matrix_t*)apSlmn)->nb;
    int Ie_M = ((parsec_tiled_matrix_t*)apIe)->mb;
    int Ie_N = ((parsec_tiled_matrix_t*)apIe)->nb;
    int Slmn_M = ((parsec_tiled_matrix_t*)apSlmn)->mb;
    int Io_M = ((parsec_tiled_matrix_t*)apIo)->mb;
    int Io_N = ((parsec_tiled_matrix_t*)apIo)->mb;

    /* Create handle_cublas_gpu */
	cuDoubleComplex alpha_complex, beta_complex;
    cublasStatus_t status;

    /* Find workspace */
    WorkSpace *gpu_work_space = (WorkSpace *) apWorkSpace;
    StreamWorkSpace *stream_found = LookupGPUWorkspace(apCudaDevice, apCudaStream, gpu_work_space);

    /* Get handle_cublas */
    cublasHandle_t handle = stream_found->handle;
    cublasSetStream( handle, apCudaStream->cuda_stream );

    // Temp buffer
    cuDoubleComplex *Gmtheta_r = (cuDoubleComplex *)stream_found->Gmtheta_r;
    cuDoubleComplex *Fmnm = (cuDoubleComplex *)stream_found->Fmnm;
    cuDoubleComplex *tmp1 = (cuDoubleComplex *)stream_found->tmp1;
    cuDoubleComplex *tmp2 = (cuDoubleComplex *)stream_found->tmp2;

    assert(f_data_N == Ep_N);
    alpha_complex = make_cuDoubleComplex(1.0, 0);
    beta_complex = make_cuDoubleComplex(0.0, 0);
    int Gmtheta_r_M = f_data_M;
    int Gmtheta_r_N = Ep_N;
    int Gmtheta_r_K = f_data_N;
    cublasZgemm(handle, CUBLAS_OP_N, CUBLAS_OP_N,
            Gmtheta_r_M, Gmtheta_r_N, Gmtheta_r_K,
            &alpha_complex, apF_data, Gmtheta_r_M,
                            apEp, Gmtheta_r_K,
            &beta_complex, Gmtheta_r, Gmtheta_r_M);
    assert(Et1_M == Gmtheta_r_M);
    int Fmnm_M = Et1_M;
    int Fmnm_N = Gmtheta_r_N;
    int Fmnm_K = Gmtheta_r_M;

    cublasZgemm(handle, CUBLAS_OP_N, CUBLAS_OP_N,
            Fmnm_M, Fmnm_N, Fmnm_K,
            &alpha_complex, apEt1, Fmnm_M,
                            Gmtheta_r, Fmnm_K,
            &beta_complex, Fmnm, Fmnm_M);

	assert(Et2_M == P_M);
    int tmp1_M = Et2_M;
    int tmp1_N = P_N;
    int tmp1_K = P_M;

    cublasZgemm(handle, CUBLAS_OP_N, CUBLAS_OP_N,
            tmp1_M, tmp1_N, tmp1_K,
            &alpha_complex, apEt2, tmp1_M,
                           apP, tmp1_K,
            &beta_complex, tmp1, tmp1_M);

	assert(tmp1_N == Gmtheta_r_M);
    int tmp2_M = tmp1_M;
    int tmp2_N = Gmtheta_r_N;
    int tmp2_K = Gmtheta_r_M;

    cublasZgemm(handle, CUBLAS_OP_N, CUBLAS_OP_N,
            tmp2_M, tmp2_N, tmp2_K,
            &alpha_complex, tmp1, tmp2_M,
                            Gmtheta_r, tmp2_K,
            &beta_complex,  tmp2, tmp2_M);

	// TODO: diagonal matrix
	assert(Fmnm_M == tmp2_M);
	assert(Fmnm_N == P_N);
	assert(tmp2_N == P_M);
	Fmnm_K = tmp2_N;
	beta_complex = make_cuDoubleComplex(1.0, 0);
    cublasZgemm(handle, CUBLAS_OP_N, CUBLAS_OP_N,
            Fmnm_M, Fmnm_N, Fmnm_K,
            &alpha_complex, tmp2, Fmnm_M,
                           apD, Fmnm_K,
            &beta_complex, Fmnm, Fmnm_M);

    assert(Slmn_N == Ie_M);
    assert(Ie_N == Fmnm_M);

	int flmn_matrix_M = aL;
	int flmn_matrix_N = aL;

	cuDoubleComplex *flmn_matrix = tmp1;
	cuDoubleComplex *Fmnm_tmp;
	cuDoubleComplex *Slmn_tmp;
	cuDoubleComplex *multipy_tmp = tmp2 + Fmnm_M + Slmn_N;

	for(int m = 0; m < aL; m++) {
        Fmnm_tmp = Fmnm + (aL+m-1)*Fmnm_M;
		if( 0 == m % 2) {
			for(int n = m; n < aL; n++) {
#if 0
				index = CalculateSingleIndex(n, m);
				for(int j = 0; j < Slmn_N; j++) {
					Slmn_tmp[j] = apSlmn[j*Slmn_M+index];
				}
#endif
                Slmn_tmp = apSlmn + CalculateSingleIndex(n, m);
				alpha_complex = make_cuDoubleComplex(1.0, 0);
				beta_complex = make_cuDoubleComplex(0.0, 0);
				cublasZgemv(handle, CUBLAS_OP_N,
						Ie_M, Ie_N,
						&alpha_complex, apIe, Ie_M,
						Fmnm_tmp, 1,
						&beta_complex, multipy_tmp, 1);

				cublasZdotu(handle, Slmn_N, Slmn_tmp, Slmn_M, multipy_tmp, 1, &flmn_matrix[m*flmn_matrix_M+n]);
			}
		} else {
			for(int n = m; n < aL; n++) {
#if 0
                index = climate_emulator_getSingleIndex(n, m);
                for(int j = 0; j < Slmn_N; j++) {
                    Slmn_tmp[j] = apSlmn[j*Slmn_M+index];
                }
#endif
                Slmn_tmp = apSlmn + CalculateSingleIndex(n, m);

                alpha_complex = make_cuDoubleComplex(1.0, 0);
                beta_complex = make_cuDoubleComplex(0.0, 0);
                cublasZgemv(handle, CUBLAS_OP_N,
                        Io_M, Io_N,
                        &alpha_complex, apIo, Io_M,
                        Fmnm_tmp, 1,
                        &beta_complex, multipy_tmp, 1);
                cublasZdotu(handle, Slmn_N, Slmn_tmp, Slmn_M, multipy_tmp, 1, &flmn_matrix[m*flmn_matrix_M+n]);
			}
		}

	}

#if 0
	for(int n = 0; n < aL; n++) {
		for(int m = 0; m <= n; m++) {
			apFlm[n*n+n+m] = cuCreal(flmn_matrix[m*flmn_matrix_M+n]);
            if( m != 0) {
				apFlm[n*n+n-m] = cuCimag(flmn_matrix[m*flmn_matrix_M+n]);
			}
		}
	}
#endif
    gb24_reshape_GPU(apFlm, flmn_matrix, aL, flmn_matrix_M, apCudaStream->cuda_stream);
}
#endif