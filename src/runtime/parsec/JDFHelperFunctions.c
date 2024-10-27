
// Copyright (c) 2017-2024 King Abdullah University of Science and Technology,
// All rights reserved.
// ExaGeoStat is a software package, provided by King Abdullah University of Science and Technology (KAUST).

/**
 * @file JDFHelperFunctions.c
 * @brief Implementation of JDF helper functions.
 * @version 2.0.0
 * @author Mahmoud ElKarargy
 * @date 2024-10-20
**/

#include <runtime/parsec/jdf/JDFHelperFunctions.h>

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

    int Zlm_M = ((parsec_matrix_block_cyclic_t *) apZlm)->super.mb;

    int f_spatial_M = ((parsec_matrix_block_cyclic_t *) apZlm)->super.mb;
    int f_spatial_N = ((parsec_matrix_block_cyclic_t *) apZlm)->super.nb;

    int index_Zlm, index_flm;

    int Smt_M = aL + 1;
    int Smt_N = 2 * aL - 1;
    memset(apSmt, 0, Smt_M * Smt_N * sizeof(double));

    for (int m = -(aL - 1); m < aL; m++) {
        for (int n = abs(m); n < aL; n++) {
            index_Zlm = CalculateSingleIndex(n, abs(m));
            index_flm = n * n + n + m;
            cblas_daxpy(Smt_M, apFlm[index_flm], apZlm + index_Zlm * Zlm_M, 1, apSmt + (m + aL - 1) * Smt_M, 1);
        }
    }

    int SC_M = ((parsec_matrix_block_cyclic_t *) apZlm)->super.mb;
    int SC_N = ((parsec_matrix_block_cyclic_t *) apZlm)->super.nb;
    assert(Smt_N == SC_M);
    assert(f_spatial_M == Smt_M);
    assert(f_spatial_N == SC_N);
    cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans,
                f_spatial_M, f_spatial_N, Smt_N,
                (double) 1.0, apSmt, f_spatial_M,
                apSC, Smt_N,
                (double) 0.0, apFspatial, f_spatial_M);

}

void WriteMatrixToBinaryFile(const char *apFilename, size_t aM, size_t aN, double *apDataMatrix) {
    FILE *file = fopen(apFilename, "wb"); // Open file for writing in binary mode
    if (file == NULL) {
        printf("Unable to open file for writing.\n");
        return;
    }

    // Write matrix data column by column
    for (size_t j = 0; j < aN; j++) {
        fwrite(apDataMatrix + j * aM, sizeof(double), aM, file);
    }

    fclose(file);
    printf("Matrix written to binary file successfully.\n");
}

void
WriteData(parsec_context_t *apContext, const char *apFilename, parsec_tiled_matrix_t *apParsecMatrix,
          hicma_parsec_data_t *apHicmaParsecData, int aSymm, int aProcessRank, int aNodes) {
    int P = ((parsec_matrix_block_cyclic_t *) apParsecMatrix)->grid.rows;
    int Q = ((parsec_matrix_block_cyclic_t *) apParsecMatrix)->grid.cols;
    int MB = apParsecMatrix->mb;
    int NB = apParsecMatrix->nb;
    int M = apParsecMatrix->lm;
    int N = apParsecMatrix->ln;

    if (aProcessRank == 0) {
        fprintf(stderr, RED "\nWrite matrix to file:\n" RESET);
        fprintf(stderr, RED "\nP %d Q %d MB %d NB %d M %d N %d:\n" RESET, P, Q, MB, NB, M, N);
    }

#if MPIIO
    MPI_Datatype darrayA;

    CreateMPIDatatype2D(aNodes, aProcessRank,
            //N, M, MPI_DISTRIBUTE_CYCLIC, MPI_DISTRIBUTE_CYCLIC,
                        M, N, MPI_DISTRIBUTE_BLOCK, MPI_DISTRIBUTE_BLOCK,
                        MB, NB, P, Q,
                        MPI_DOUBLE, &darrayA);

    double *DA = (double *) parsec_data_allocate((size_t) apParsecMatrix->nb_local_tiles *
                                                 (size_t) apParsecMatrix->bsiz *
                                                 (size_t) parsec_datadist_getsizeoftype(apParsecMatrix->mtype));

//        if(symm)
//            tile_to_lapack_sym( parsec, (parsec_tiled_matrix_t *)A, DA, P, nodes/P);
//        else
//            tile_to_lapack( parsec, (parsec_tiled_matrix_t *)A, DA, P, nodes/P);

    MPIWriteDFile(apFilename, aProcessRank, DA, apParsecMatrix->llm * apParsecMatrix->lln, MPI_DOUBLE, darrayA);
    parsec_data_free(DA);
#else

    double *DA = (double *) malloc((size_t) apParsecMatrix->llm * (size_t) apParsecMatrix->lln *
                                   (size_t) parsec_datadist_getsizeoftype(apParsecMatrix->mtype));

//    if(symm)
//        tile_to_lapack_sym( parsec, (parsec_tiled_matrix_t *)A, DA, P, nodes/P);
//    else
//        tile_to_lapack( parsec, (parsec_tiled_matrix_t *)A, DA, P, nodes/P);

    WriteMatrixToBinaryFile(apFilename, N, N, DA);

    free(DA);
#endif
}

void ReadData(parsec_context_t *apContext,
              const char *aFileName,
              parsec_tiled_matrix_t *apParsecMatrix,
              hicma_parsec_data_t *apHicmaParsecData,
              int aP, int aQ,
              int aSymm, int aRank, int aNodes) {
    int MB = apParsecMatrix->mb;
    int NB = apParsecMatrix->nb;
    int M = apParsecMatrix->lm;
    int N = apParsecMatrix->ln;

    if (aRank == 0) fprintf(stderr, RED "\nRead matrix from file:\n" RESET);

    MPI_Datatype darrayA;

    CreateMPIDatatype2D(aNodes, aRank,
                        M, N, MPI_DISTRIBUTE_CYCLIC, MPI_DISTRIBUTE_CYCLIC,
            //M, N, MPI_DISTRIBUTE_BLOCK, MPI_DISTRIBUTE_BLOCK,
                        MB, NB, aP, aQ,
                        MPI_DOUBLE, &darrayA);

    double *DA = (double *) parsec_data_allocate((size_t) apParsecMatrix->nb_local_tiles *
                                                 (size_t) apParsecMatrix->bsiz *
                                                 (size_t) parsec_datadist_getsizeoftype(apParsecMatrix->mtype));

    MPIReadDFile(aFileName, aRank, DA, apParsecMatrix->llm * apParsecMatrix->lln, MPI_DOUBLE, darrayA);
//    lapack_to_tile( apContext, apParsecMatrix, DA, aP, aQ);
    parsec_data_free(DA);
}

size_t PARSECGetAddressCM(parsec_tiled_matrix_t *apDescA, int aM, int aN, int aP, int aQ) {
    size_t mm = aM + apDescA->i / apDescA->mb;
    size_t nn = aN + apDescA->j / apDescA->nb;
    size_t eltsize = (size_t) parsec_datadist_getsizeoftype(apDescA->mtype);
    size_t offset = 0;

    mm = mm / aP;
    nn = nn / aQ;

    offset = (size_t) (apDescA->llm * apDescA->nb) * nn + (size_t) (apDescA->mb) * mm;
    return (offset * eltsize);
}

static void CreateMPIDatatype2D(int aProcesses, int aMPIRank,
                                int M, int N, int aDistM, int aDistN,
                                int mb, int nb, int aPRow, int aPCol,
                                MPI_Datatype aOldType, MPI_Datatype *apNewType) {
    int pdims[2], distribs[2], dims[2], dargs[2];
    int ierr;

    pdims[0] = aPRow;
    pdims[1] = aPCol;

    dims[0] = M;
    dims[1] = N;

    distribs[0] = aDistM;
    distribs[1] = aDistN;

    dargs[0] = mb;
    dargs[1] = nb;

    ierr = MPI_Type_create_darray(aProcesses, aMPIRank, 2,
                                  dims, distribs, dargs, pdims, MPI_ORDER_FORTRAN,
                                  aOldType, apNewType);

    if (ierr != 0) {
        printf("\n ====> ierr MPI_Type_create_darray:%d\n", ierr);
        fflush(stdout);
    }

    ierr = MPI_Type_commit(apNewType);
    if (ierr != 0) {
        printf("\n ====> ierr MPI_Type_create_darray:%d\n", ierr);
        fflush(stdout);
    }
}

void MPIReadDFile(const char *apFilename, int aMPIRank, double *apDataMatrix, MPI_Count aSize, MPI_Datatype aOldType,
                  MPI_Datatype aNewType) {
    int ierr;
    MPI_File infile;
    MPI_Status mpistatus;

    MPI_Aint file_type_extent, lb;
    MPI_Count file_type_size;
    MPI_Type_get_extent(aNewType, &lb, &file_type_extent);
    MPI_Type_size_x(aNewType, &file_type_size);
    uint64_t buffer_size = file_type_size / sizeof(double);

    if (buffer_size != aSize) {
        printf("\n On myrank:%d reading data on: %s parsec size:%lld not matching mpi local:%lu\n", aMPIRank,
               apFilename, aSize, buffer_size);
        fflush(stdout);
    }


    ierr = MPI_File_open(MPI_COMM_WORLD, apFilename, MPI_MODE_RDONLY,
                         MPI_INFO_NULL, &infile);
    if (ierr != 0) {
        printf("\n ====> MPI_File_open ierr:%d\n", ierr);
        fflush(stdout);
    }
    ierr = MPI_File_set_view(infile, 0, aOldType, aNewType, "native", MPI_INFO_NULL);

    if (ierr != 0) {
        printf("\n ====> MPI_File_set_view ierr:%d\n", ierr);
        fflush(stdout);
    }

    ierr = MPI_File_read_all(infile, apDataMatrix, aSize, aOldType, MPI_STATUS_IGNORE);

    if (ierr != 0) {
        printf("\n ====> MPI_File_read_all ierr:%d\n", ierr);
        fflush(stdout);
    }

    ierr = MPI_File_close(&infile);
    if (ierr != 0) {
        printf("\n ====> MPI_File_close ierr:%d\n", ierr);
        fflush(stdout);
    }

}

void MPIWriteDFile(const char *apFilename, int aMPIRank, double *apDataMatrix, MPI_Count aSize, MPI_Datatype aOldType,
                   MPI_Datatype aNewType) {
    int ierr;
    MPI_File infile;
    MPI_Status mpistatus;

    MPI_Aint file_type_extent, lb;
    MPI_Count file_type_size;
    MPI_Type_get_extent(aNewType, &lb, &file_type_extent);
    MPI_Type_size_x(aNewType, &file_type_size);
    int64_t buffer_size = file_type_size / sizeof(double);

    if (buffer_size != aSize) {
        printf("\n On myrank:%d reading data on: %s parsec size:%lld not matching mpi local:%ld\n", aMPIRank, apFilename,
               aSize, buffer_size);

        fflush(stdout);
    }

    ierr = MPI_File_open(MPI_COMM_WORLD, apFilename, MPI_MODE_CREATE | MPI_MODE_WRONLY,
                         MPI_INFO_NULL, &infile);
    if (ierr != 0) {
        printf("\n ====> MPI_File_open ierr:%d\n", ierr);
        fflush(stdout);
    }

    ierr = MPI_File_set_view(infile, 0, aOldType, aNewType, "native", MPI_INFO_NULL);
    if (ierr != 0) {
        printf("\n ====> MPI_File_set_view ierr:%d\n", ierr);
        fflush(stdout);
    }

    ierr = MPI_File_write_all(infile, apDataMatrix, aSize, aOldType, &mpistatus);
    if (ierr != 0) {
        printf("\n ====> MPI_File_write_all ierr:%d\n", ierr);
        fflush(stdout);
    }

    ierr = MPI_File_close(&infile);
    if (ierr != 0) {
        printf("\n ====> MPI_File_close ierr:%d\n", ierr);
        fflush(stdout);
    }
}
