
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
using namespace exageostat::configurations;
using namespace exageostat::dataunits;

template<typename T>
void Prediction<T>::PredictMissingData(const hardware::ExaGeoStatHardware &aHardware, ExaGeoStatData<T> &aData,
                                       Configurations &aConfigurations, T *apMeasurementsMatrix,
                                       const kernels::Kernel<T> &aKernel) {

    bool can_predict = true;
    int num_params = aKernel.GetParametersNumbers();
    for (int i = 0; i < num_params; i++) {
        if (aConfigurations.GetEstimatedTheta()[i] == -1) {
            can_predict = false;
            break;
        }
    }
    if (!can_predict && (aConfigurations.GetIsMLOEMMOM() || aConfigurations.GetIsMSPE())) {
        throw std::runtime_error(
                "Can't predict without an estimated theta, please either pass --etheta or run the modeling module before prediction");
    }

    //If the number of missed locations isn't a positive value, No prediction needed.
    if (aConfigurations.GetUnknownObservationsNb() <= 0) {
        return;
    }

    int i;
    int number_of_mspe = 3;
    int p = aConfigurations.GetP();
    int z_miss_number = aConfigurations.GetUnknownObservationsNb();
    int n_z_obs = aConfigurations.CalculateZObsNumber();

    LOGGER("- Number of Z observations: " << aConfigurations.GetP() * n_z_obs)
    results::Results::GetInstance()->SetZMiss(aConfigurations.GetP() * n_z_obs);
    T *z_obs = new T[p * n_z_obs];
    T *z_miss = new T[p * z_miss_number];
    T *z_actual = new T[p * z_miss_number];
    std::vector<T> avg_pred_value(number_of_mspe);
    auto miss_locations = new Locations<T>(z_miss_number, aData.GetLocations()->GetDimension());
    auto obs_locations = new Locations<T>(n_z_obs, aData.GetLocations()->GetDimension());
    // We Predict date with only Exact computation. This is a pre-request.
    auto linear_algebra_solver = linearAlgebra::LinearAlgebraFactory<T>::CreateLinearAlgebraSolver(common::EXACT_DENSE);
    InitializePredictionArguments(aConfigurations, aData, linear_algebra_solver, z_obs, z_actual, *miss_locations,
                                  *obs_locations, apMeasurementsMatrix);

    // MLOE MMOM Auxiliary Function Call
    if (aConfigurations.GetIsMLOEMMOM()) {
        LOGGER("---- Using Auxiliary Function MLOE MMOM ----")
        linear_algebra_solver->ExaGeoStatMLETileMLOEMMOM(aConfigurations, aData, aHardware,
                                                         (T *) aConfigurations.GetInitialTheta().data(),
                                                         (T *) aConfigurations.GetEstimatedTheta().data(),
                                                         *miss_locations, *obs_locations, aKernel);
        LOGGER(" ---- mloe_mmom Time(main): %6.2f seconds")
    }

    // IDW Auxiliary Function Call
    if (aConfigurations.GetIsIDW()) {
        LOGGER("---- Using Auxiliary Function IDW ----")
        T *mspe = new T[number_of_mspe];
        PredictionAuxiliaryFunctions<T>::PredictIDW(z_miss, z_actual, z_obs, z_miss_number, n_z_obs, *miss_locations,
                                                    *obs_locations, mspe);

        LOGGER("- Average prediction Error (IDW): ", true)
        for (i = 0; i < number_of_mspe; i++) {
            avg_pred_value[i] += mspe[i];
            LOGGER_PRECISION(avg_pred_value[i], 9)
            if (i != number_of_mspe - 1) {
                LOGGER_PRECISION(", ")
            }
        }
        results::Results::GetInstance()->SetIDWError(
                std::vector<double>(avg_pred_value.data(), avg_pred_value.data() + number_of_mspe));
        LOGGER("")
        delete[] mspe;
    }

    // MSPE Prediction Function Call
    if (aConfigurations.GetIsMSPE()) {
        LOGGER("---- Using MSPE ----")
        T *prediction_error_mspe = linear_algebra_solver->ExaGeoStatMLEPredictTile(aData,
                                                                                   (T *) aConfigurations.GetEstimatedTheta().data(),
                                                                                   z_miss_number, n_z_obs, z_obs,
                                                                                   z_actual, z_miss, aHardware,
                                                                                   aConfigurations, *miss_locations,
                                                                                   *obs_locations, aKernel);
        for (i = 0; i < number_of_mspe; i++) {
            avg_pred_value[i] += prediction_error_mspe[i];
        }
        LOGGER("- Average prediction Error (mspe): ", true)
        for (i = 0; i < number_of_mspe; i++) {
            LOGGER_PRECISION(avg_pred_value[i] << "\t", 8)
        }
        LOGGER("")
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
                                                  std::unique_ptr<linearAlgebra::LinearAlgebraMethods<T>> &aLinearAlgebraSolver,
                                                  T *apZObs, T *apZActual, Locations<T> &aMissLocation,
                                                  Locations<T> &aObsLocation, T *apMeasurementsMatrix) {

    int N = aConfigurations.GetProblemSize();
    T *z = new T[N];

    aLinearAlgebraSolver->ExaGeoStatGetZObs(aConfigurations, z, N, *aData.GetDescriptorData(), apMeasurementsMatrix);
    PredictionHelpers<T>::PickRandomPoints(aConfigurations, aData, apZObs, apZActual, z, aMissLocation, aObsLocation);
    delete[] z;
}