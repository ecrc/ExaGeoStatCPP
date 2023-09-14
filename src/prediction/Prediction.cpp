
// Copyright (c) 2017-2023 King Abdullah University of Science and Technology,
// All rights reserved.
// ExaGeoStat is a software package, provided by King Abdullah University of Science and Technology (KAUST).

/**
 * @file Prediction.cpp
 * @brief Contains the implementation of the Prediction class.
 * @version 1.0.0
 * @author Mahmoud ElKarargy
 * @author Sameh Abdulah
 * @date 2023-06-08
**/

#include <prediction/Prediction.hpp>
#include <prediction/PredictionHelpers.hpp>
#include <prediction/PredictionAuxiliaryFunctions.hpp>
#include <linear-algebra-solvers/LinearAlgebraFactory.hpp>

using namespace exageostat::prediction;
using namespace exageostat::hardware;
using namespace exageostat::configurations;
using namespace exageostat::dataunits;
using namespace exageostat::linearAlgebra;

template<typename T>
void
Prediction<T>::PredictMissingData(const ExaGeoStatHardware &aHardware, ExaGeoStatData<T> &aData, Configurations &aConfigurations) {

    //If the number of missed locations isn't a positive value, No prediction needed.
    if(aConfigurations.GetUnknownObservationsNb() <= 0){
        return;
    }

    int i;
    int number_of_mspe = 3;
    int p = aConfigurations.GetP();
    int z_miss_number = aConfigurations.GetUnknownObservationsNb();
    int n_z_obs = aConfigurations.CalculateZObsNumber();

    printf("nZobs = %d\n", aConfigurations.GetP() * n_z_obs);

    T *z_obs = new T[p * n_z_obs];
    T *z_miss = new T[p * z_miss_number];
    T *z_actual = new T[p * z_miss_number];

    std::vector<T> avg_pred_value(number_of_mspe);
    auto miss_locations = new Locations<T>(z_miss_number, aData.GetLocations()->GetDimension());
    auto obs_locations = new Locations<T>(n_z_obs, aData.GetLocations()->GetDimension());
    auto linear_algebra_solver = linearAlgebra::LinearAlgebraFactory<T>::CreateLinearAlgebraSolver(
            aConfigurations.GetComputation());

    InitializePredictionArguments(aConfigurations, aData, linear_algebra_solver, z_obs, z_actual, *miss_locations,
                                  *obs_locations);
    PredictionAuxiliaryFunctions<T> prediction_auxiliary_functions;

    if (aConfigurations.GetIsIDW()) {
        T *mspe = new T[number_of_mspe];
        prediction_auxiliary_functions.PredictTileIDW(z_miss, z_actual, z_obs, z_miss_number, n_z_obs, *miss_locations, *obs_locations, mspe);

        fprintf(stderr, "\nAverage prediction Error (IDW): ");
        for (i = 0; i < number_of_mspe; i++) {
            avg_pred_value[i] += mspe[i];
            fprintf(stderr, "   %3.9f   ", avg_pred_value[i]);
        }
        fprintf(stderr, "\n");
        delete[] mspe;
    }

    if (aConfigurations.GetIsMSPE()) {
        T *prediction_error_mspe = linear_algebra_solver->ExaGeoStatMLEPredictTILE(aData,
                                                                                   (T *) aConfigurations.GetStartingTheta().data(),
                                                                                   z_miss_number, n_z_obs, z_obs,
                                                                                   z_actual, z_miss, aHardware,
                                                                                   aConfigurations, *miss_locations,
                                                                                   *obs_locations);
        for (i = 0; i < number_of_mspe; i++) {
            avg_pred_value[i] += prediction_error_mspe[i];
        }

        fprintf(stderr, "\nAverage prediction Error (mspe): ");
        for (i = 0; i < number_of_mspe; i++) {
            fprintf(stderr, "   %3.9f   ", avg_pred_value[i]);
        }
        fprintf(stderr, "\n");
        delete[] prediction_error_mspe;
    }

    delete[] z_obs;
    delete[] z_miss;
    delete[] z_actual;
    delete miss_locations;
    delete obs_locations;
}

template<typename T>
void Prediction<T>::InitializePredictionArguments(Configurations &aConfigurations, ExaGeoStatData<T> &aData,
                                                  std::unique_ptr<LinearAlgebraMethods<T>> &aLinearAlgebraSolver,
                                                  T *apZObs, T *apZActual, Locations<T> &aMissLocation,
                                                  Locations<T> &aObsLocation) {

    int N = aConfigurations.GetProblemSize();
    T *z = new T[N];
    aLinearAlgebraSolver->GetZObs(z, N, *aData.GetDescriptorData());
    PredictionHelpers<T>::PickRandomPoints(aConfigurations, aData, apZObs, apZActual, z, aMissLocation, aObsLocation);

    delete[] z;
}