
// Copyright (c) 2017-2024 King Abdullah University of Science and Technology,
// All rights reserved.
// ExaGeoStat is a software package, provided by King Abdullah University of Science and Technology (KAUST).

/**
 * @file ChameleonStarPuHelpers.cpp
 * @brief A class for Chameleon implementation of StarPu helpers.
 * @version 1.1.0
 * @author Mahmoud ElKarargy
 * @date 2024-02-25
**/

#include <linear-algebra-solvers/concrete/ChameleonHeaders.hpp>
#include <runtime/starpu/helpers/concrete/ChameleonStarPuHelpers.hpp>

using namespace exageostat::runtime;

void
ChameleonStarPuHelpers::ExaGeoStatOptionsInit(void *apOptions, void *apSequence, void *apRequest) {

    RUNTIME_options_init((RUNTIME_option_t *) apOptions, (CHAM_context_t *) ExaGeoStatHardware::GetChameleonContext(),
                         (RUNTIME_sequence_t *) apSequence, (RUNTIME_request_t *) apRequest);
}

void ChameleonStarPuHelpers::ExaGeoStatOptionsFree(void *apOptions) {
    RUNTIME_options_ws_free((RUNTIME_option_t *) apOptions);
}


void ChameleonStarPuHelpers::ExaGeoStatOptionsFinalize(void *apOptions) {
    auto *options = (RUNTIME_option_t *) apOptions;
    RUNTIME_options_finalize(options, (CHAM_context_t *) ExaGeoStatHardware::GetChameleonContext());
}

void *ChameleonStarPuHelpers::ExaGeoStatDataGetAddr(void *apDescriptor, const int &aDescRow, const int &aDescCol) {
    return RUNTIME_data_getaddr((CHAM_desc_t *) apDescriptor, aDescRow, aDescCol);
}

int ChameleonStarPuHelpers::GetMT(void *apDescriptor) {
    auto descriptor = (CHAM_desc_t *) apDescriptor;
    return descriptor->mt;
}

int ChameleonStarPuHelpers::GetM(void *apDescriptor) {
    auto descriptor = (CHAM_desc_t *) apDescriptor;
    return descriptor->m;
}

int ChameleonStarPuHelpers::GetMB(void *apDescriptor) {
    auto descriptor = (CHAM_desc_t *) apDescriptor;
    return descriptor->mb;
}

void *ChameleonStarPuHelpers::GetOptions() {
    return new RUNTIME_option_t;
}

void ChameleonStarPuHelpers::DeleteOptions(void *apOptions) {
    delete (RUNTIME_option_t *) apOptions;
}
