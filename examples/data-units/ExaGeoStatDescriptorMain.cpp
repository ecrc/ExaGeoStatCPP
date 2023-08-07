//// Copyright (c) 2017-2023 King Abdullah University of Science and Technology,
//// All rights reserved.
//// ExaGeoStat is a software package, provided by King Abdullah University of Science and Technology (KAUST).
//
///**
// * @file ExaGeoStatDescriptorMain.cpp
// * @brief
// * @version 1.0.0
// * @author Sameh Abdulah
// * @date 2023-07-17
//**/
//
//#include <data-units/ExaGeoStatDescriptor.hpp>
//#include <api/ExaGeoStat.hpp>
//#ifdef EXAGEOSTAT_USE_CHAMELEON
//#include <chameleon/struct.h>
//#endif
//
//#ifdef EXAGEOSTAT_USE_HiCMA
//#include <hicma_struct.h>
//#endif
//
//using namespace exageostat::api;
//using namespace exageostat::common;
//using namespace std;
//
int main(){
//
//#ifdef EXAGEOSTAT_USE_CHAMELEON
//    exageostat::configurations::Configurations::GetConfigurations()->SetComputation(exageostat::common::EXACT_DENSE);
//#endif
//
//#ifdef EXAGEOSTAT_USE_HiCMA
//    exageostat::configurations::Configurations::GetConfigurations()->SetComputation(exageostat::common::TILE_LOW_RANK);
//#endif
//    cout << "** Initialise ExaGeoStat hardware ** " << endl;
//    // Initialize ExaGeoStat hardware
//    ExaGeoStat<double>::ExaGeoStatInitializeHardware();
//
//    cout << "** Initialise ExaGeoStat Descriptor ** " << endl;
//    // Create stack variable of ExaGeoStatDescriptor
//    exageostat::dataunits::ExaGeoStatDescriptor<double> exaGeoStatDescriptor{};
//    // Set values of the descriptor.
//
//    exaGeoStatDescriptor.mIsOOC = false;
//    auto *x = new double[5];
//    for(int i = 0; i < 5; i ++){
//        x[i] = i * i;
//    }
//    exaGeoStatDescriptor.mat = x;
//    exaGeoStatDescriptor.dtyp = EXAGEOSTAT_REAL_DOUBLE;
//    exaGeoStatDescriptor.mb = 5;
//    exaGeoStatDescriptor.nb = 5;
//    exaGeoStatDescriptor.bsiz = 5 * 5;
//    exaGeoStatDescriptor.lm = 10;
//    exaGeoStatDescriptor.ln = 10;
//    exaGeoStatDescriptor.i = 0;
//    exaGeoStatDescriptor.j = 0;
//    exaGeoStatDescriptor.m = 10;
//    exaGeoStatDescriptor.n = 10;
//    exaGeoStatDescriptor.p = 1;
//    exaGeoStatDescriptor.q = 1;
//
//#ifdef EXAGEOSTAT_USE_CHAMELEON
//    cout << "** Creating Chameleon Descriptor ** " << endl;
//    CHAM_desc_t *desc = exaGeoStatDescriptor.CreateChameleonDescriptor(exageostat::common::EXACT_DENSE);
//#endif
//#ifdef EXAGEOSTAT_USE_HiCMA
//    cout << "** Creating Hicma Descriptor ** " << endl;
//    HICMA_desc_t *desc = exaGeoStatDescriptor.CreateHicmaDescriptor(exageostat::common::TILE_LOW_RANK);
//#endif
//
//    cout << "** ExaGeoStat Descriptor Content** " << endl;
//    cout << "\tmb: " << desc->mb << "\tnb: " << desc->nb << endl;
//    cout << "\tlm: " << desc->lm << "\tln: " << desc->ln << endl;
//    cout << "\ti: " << desc->i << "\tj: " << desc->j << endl;
//    cout << "\tm: " << desc->m << "\tn: " << desc->n << endl;
//    cout << "\tp: " << desc->p << "\tq: " << desc->q << endl;
//    cout << "\tbsize: " << desc->bsiz << endl;
//    cout << "\tmatrix: ";
//    for(int i = 0; i < 5; i ++){
//        cout << ((double*) (desc->mat))[i] << " ";
//    }
//    cout << endl;
//    cout << "** Finalize ExaGeoStat Hardware** " << endl;
////    ExaGeoStat<double>::ExaGeoStatFinalizeHardware();
}

