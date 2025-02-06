
// Copyright (c) 2017-2024 King Abdullah University of Science and Technology,
// All rights reserved.
// ExaGeoStat is a software package, provided by King Abdullah University of Science and Technology (KAUST).

/**
 * @file GPUHelperFunctions.c
 * @brief Implementation of GPU helper functions.
 * @version 2.0.0
 * @author Mahmoud ElKarargy
 * @author Sameh Abdulah
 * @author Qinglei Cao
 * @date 2024-10-20
**/

#include <runtime/parsec/GPUHelperFunctions.h>

StreamWorkSpace *LookupGPUWorkspace(parsec_device_cuda_module_t *apCudaDevice,
                                    parsec_cuda_exec_stream_t *apCudaSream,
                                    WorkSpace *apWorkSPace){
    int i, j;
    /* Look for device */
    for(i = 0; i < (int)parsec_nb_devices; i++) {
        parsec_device_module_t *device = parsec_mca_device_get(i);
        if( NULL == device ) continue;
        if( device->type != PARSEC_DEV_CUDA && device->type != PARSEC_DEV_HIP ) continue;
        parsec_device_cuda_module_t *cuda_device_compare = (parsec_device_cuda_module_t*)device;

        if(apCudaDevice->cuda_index == cuda_device_compare->cuda_index)
            break;
    }

    /* Look for stream; 0, h2d; 1 d2h*/
    for(j = 2; j < apCudaDevice->super.max_exec_streams; j++) {
        if( apCudaSream == (parsec_cuda_exec_stream_t *)apCudaDevice->super.exec_stream[j])
            break;
    }

    return &apWorkSPace->gpu_workspace[i].stream_workspace[j];
}

/* Allocate memory for workspace */
static void WorkspaceMemoryAllocate(WorkSpace **apWorkSpace){
    *apWorkSpace = (WorkSpace *)malloc( sizeof(WorkSpace) );
    (*apWorkSpace)->gpu_workspace = (GPUWorkSpace *)malloc( parsec_nb_devices * sizeof(GPUWorkSpace));

    for( int i = 0; i < parsec_nb_devices; i++ ) {
        (*apWorkSpace)->gpu_workspace[i].stream_workspace = (StreamWorkSpace *)malloc( PARSEC_GPU_MAX_STREAMS * sizeof(StreamWorkSpace));
        (*apWorkSpace)->gpu_workspace[i].cuda_device = (parsec_device_cuda_module_t *)malloc( sizeof(parsec_device_cuda_module_t) );
    }
}

WorkSpace *InitializeWorkSpace(parsec_context_t *apContext, int aFDatataM, int aEpN, int aEt1M, int aEt2M, int aPN){
    /* Allocate memory */
    WorkSpace * pWork_space;
    WorkspaceMemoryAllocate(&pWork_space);

    /* Traverse all gpu device */
    for(int i = 0; i < (int)parsec_nb_devices; i++) {
        parsec_device_module_t *device = parsec_mca_device_get(i);
        if( NULL == device ) continue;
        if( device->type != PARSEC_DEV_CUDA && device->type != PARSEC_DEV_HIP ) continue;

        /* Set cuda_device */
        parsec_device_cuda_module_t *aCuda_device = (parsec_device_cuda_module_t*)device;
        cudaSetDevice(aCuda_device->cuda_index);
        pWork_space->gpu_workspace[i].cuda_device = aCuda_device;

        /* Traverse all streams */
        for(int j = 0; j < aCuda_device->super.max_exec_streams; j++) {
            /* j 0, h2d; j 1, d2h */
            if( j <= 1 ) continue;

            parsec_cuda_exec_stream_t* cuda_stream = (parsec_cuda_exec_stream_t*)aCuda_device->super.exec_stream[j];

            cublasStatus_t status;
            cublasHandle_t handle;

            /* Create handle_cublas_gpu */
            status = cublasCreate(&handle);
            cublasSetStream( handle, cuda_stream->cuda_stream );
            printf("statues of the assert is %d\n", status);
//            assert(CUBLAS_STATUS_SUCCESS == status);
            pWork_space->gpu_workspace[i].stream_workspace[j].handle = handle;

            /* Buffer */
            pWork_space->gpu_workspace[i].stream_workspace[j].Gmtheta_r = zone_malloc( aCuda_device->super.memory, aFDatataM * aEpN * sizeof(cuDoubleComplex) );
            assert(NULL != pWork_space->gpu_workspace[i].stream_workspace[j].Gmtheta_r);

            pWork_space->gpu_workspace[i].stream_workspace[j].Fmnm = zone_malloc( aCuda_device->super.memory, aEt1M * aEpN * sizeof(cuDoubleComplex) );
            assert(NULL != pWork_space->gpu_workspace[i].stream_workspace[j].Fmnm);

            pWork_space->gpu_workspace[i].stream_workspace[j].tmp1 = zone_malloc( aCuda_device->super.memory, aEt2M * aPN * sizeof(cuDoubleComplex) );
            assert(NULL != pWork_space->gpu_workspace[i].stream_workspace[j].tmp1);

            pWork_space->gpu_workspace[i].stream_workspace[j].tmp2 = zone_malloc( aCuda_device->super.memory, aEt2M * aEpN * sizeof(cuDoubleComplex) );
            assert(NULL != pWork_space->gpu_workspace[i].stream_workspace[j].tmp2);

        }
    }

    return pWork_space;

}
