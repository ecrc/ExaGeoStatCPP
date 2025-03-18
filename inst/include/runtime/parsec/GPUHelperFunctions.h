/**
 * @file GPUHelperFunctions.h
 * @brief Header file for GPU helper function declarations and related data structures.
 * @details Contains definitions of structures and declarations of functions used in
 *          GPU-accelerated versions of computations. These functions manage and allocate
 *          GPU memory resources, handle CUDA streams, and set up GPU workspaces.
 * @version 2.0.0
 * @date 2024-10-20
 * @authod Mahmoud ElKarargy
 */

#ifndef GPUHelperFunctions_H
#define GPUHelperFunctions_H

#ifdef __cplusplus
extern "C" {
#endif

#include <runtime/parsec/ParsecHeader.h>
#include <complex.h>

/**
 * @struct StreamWorkSpace
 * @brief Manages CUDA resources and intermediate matrices for GPU computations.
 */
typedef struct StreamWorkSpace {
    ///Handle for CUBLAS operations on GPU
    cublasHandle_t handle;
    ///Pointer to matrix Gmtheta_r, size: gb->f_data_M * gb->Ep_N
    cuDoubleComplex *Gmtheta_r;
    ///Pointer to matrix Fmnm, size: gb->Et1_M * gb->Ep_N
    cuDoubleComplex *Fmnm;
    ///Temporary matrix tmp1, size: gb->Et2_M * gb->P_N
    cuDoubleComplex *tmp1;
    ///Temporary matrix tmp2, size: gb->Et2_M * gb->Ep_N
    cuDoubleComplex *tmp2;
} StreamWorkSpace;

/**
 * @struct GPUWorkSpace
 * @brief Encapsulates GPU workspace resources, including stream workspace and device information.
 */
typedef struct GPUWorkSpace {
    ///Pointer to the stream workspace
    StreamWorkSpace *stream_workspace;
    ///CUDA device module for GPU tasks
    parsec_device_cuda_module_t *cuda_device;
} GPUWorkSpace;

/**
 * @struct WorkSpace
 * @brief Provides a unified workspace structure including GPU resources and status info.
 */
typedef struct WorkSpace {
    ///Pointer to GPU workspace structure */
    GPUWorkSpace *gpu_workspace;
    ///Status information or error code */
    int info;
} WorkSpace;

/**
 * @brief Retrieves a GPU StreamWorkSpace for the specified CUDA device and execution stream.
 * @param[in] apCudaDevice Pointer to the CUDA device module.
 * @param[in] apCudaSream Pointer to the CUDA execution stream.
 * @param[in] apWorkSPace Pointer to the main workspace structure.
 * @return A pointer to the corresponding StreamWorkSpace for the specified device and stream.
 *
 */
StreamWorkSpace *LookupGPUWorkspace(parsec_device_cuda_module_t *apCudaDevice,
                                    parsec_cuda_exec_stream_t *apCudaSream,
                                    WorkSpace *apWorkSPace);

/**
 * @brief Allocates memory for the WorkSpace structure and its substructures.
 * @param[out] apWorkSpace Address of the pointer to the workspace structure to allocate.
 *
 */
static void WorkspaceMemoryAllocate(WorkSpace **apWorkSpace);

/**
 * @brief Initializes a WorkSpace structure with specified dimensions and parameters.
 * @param[in] parsec Pointer to the main Parsec context.
 * @param[in] f_data_M Integer representing the aFDatataM dimension.
 * @param[in] Ep_N Integer representing the aEpN dimension.
 * @param[in] Et1_M Integer representing the aEt1M dimension.
 * @param[in] Et2_M Integer representing the aEt2M dimension.
 * @param[in] P_N Integer representing the aPN dimension.
 * @return Pointer to the initialized WorkSpace structure.
 *
 */
WorkSpace *InitializeWorkSpace(parsec_context_t *apContext, int aFDatataM, int aEpN, int aEt1M, int aEt2M, int aPN);

#ifdef __cplusplus
}
#endif

#endif // GPUHelperFunctions_H
