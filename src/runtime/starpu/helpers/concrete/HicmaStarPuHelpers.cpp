
// Copyright (c) 2017-2023 King Abdullah University of Science and Technology,
// All rights reserved.
// ExaGeoStat is a software package, provided by King Abdullah University of Science and Technology (KAUST).

/**
 * @file HicmaStarPuHelpers.cpp
 * @brief A class for Hicma implementation of StarPu helpers.
 * @version 1.1.0
 * @author Mahmoud ElKarargy
 * @date 2024-02-25
**/

#include <linear-algebra-solvers/concrete/HicmaHeaders.hpp>
#include <runtime/starpu/helpers/concrete/HicmaStarPuHelpers.hpp>

using namespace exageostat::runtime;

void HicmaStarPuHelpers::ExaGeoStatOptionsInit(void *apOptions, void *apSequence, void *apRequest) {
    HICMA_RUNTIME_options_init((HICMA_option_t *) apOptions, (HICMA_context_t *) ExaGeoStatHardware::GetHicmaContext(),
                               (HICMA_sequence_t *) apSequence, (HICMA_request_t *) apRequest);
}

void HicmaStarPuHelpers::ExaGeoStatOptionsFree(void *apOptions) {
    HICMA_RUNTIME_options_ws_free((HICMA_option_t *) apOptions);
}

void HicmaStarPuHelpers::ExaGeoStatOptionsFinalize(void *apOptions) {
    auto *options = (HICMA_option_t *) apOptions;
    HICMA_RUNTIME_options_finalize(options, (HICMA_context_t *) ExaGeoStatHardware::GetHicmaContext());
}

void *HicmaStarPuHelpers::ExaGeoStatDataGetAddr(void *apDescriptor, const int &aDescRow, const int &aDescCol) {
    return HICMA_RUNTIME_data_getaddr((HICMA_desc_t *) apDescriptor, aDescRow, aDescCol);
}

int HicmaStarPuHelpers::GetMT(void *apDescriptor) {
    auto descriptor = (HICMA_desc_t *) apDescriptor;
    return descriptor->mt;
}

int HicmaStarPuHelpers::GetM(void *apDescriptor) {
    auto descriptor = (HICMA_desc_t *) apDescriptor;
    return descriptor->m;
}

int HicmaStarPuHelpers::GetMB(void *apDescriptor) {
    auto descriptor = (HICMA_desc_t *) apDescriptor;
    return descriptor->mb;
}

void *HicmaStarPuHelpers::GetOptions() {
    return new HICMA_option_t;
}

void HicmaStarPuHelpers::DeleteOptions(void *apOptions) {
    delete (HICMA_option_t *) apOptions;
}
