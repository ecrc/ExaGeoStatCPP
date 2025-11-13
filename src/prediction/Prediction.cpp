
// Copyright (c) 2017-2024 King Abdullah University of Science and Technology,
// All rights reserved.
// ExaGeoStat is a software package, provided by King Abdullah University of Science and Technology (KAUST).

/**
 * @file Prediction.cpp
 * @brief Contains the implementation of the Prediction class.
 * @version 1.1.0
 * @author Mahmoud ElKarargy
 * @author Sameh Abdulah
 * @date 2024-02-04
**/

#include <cstring>

#ifdef USE_MKL
#include <mkl_service.h>
#endif

#include <prediction/Prediction.hpp>
#include <prediction/PredictionHelpers.hpp>
#include <prediction/PredictionAuxiliaryFunctions.hpp>
#include <linear-algebra-solvers/LinearAlgebraFactory.hpp>

using namespace std;

using namespace exageostat::prediction;
using namespace exageostat::dataunits;
using namespace exageostat::results;
using namespace exageostat::configurations;

template<typename T>
void Prediction<T>::PredictMissingData(unique_ptr<ExaGeoStatData<T>> &aData, Configurations &aConfigurations,
                                       T *apMeasurementsMatrix, const kernels::Kernel<T> &aKernel,
                                       Locations<T> *apTrainLocations, Locations<T> *apTestLocations) {

#if DEFAULT_RUNTIME
    int i, j;
    bool can_predict = true;
    int num_params = aKernel.GetParametersNumbers();
    for (i = 0; i < num_params; i++) {
        if (aConfigurations.GetEstimatedTheta()[i] == -1) {
            can_predict = false;
            break;
        }
    }

    if (!can_predict &&
        (aConfigurations.GetIsMLOEMMOM() || aConfigurations.GetIsMSPE() || aConfigurations.GetIsFisher())) {
        throw runtime_error(
                "Can't predict without an estimated theta, please either pass --etheta or run the modeling module before prediction");
    }

    int number_of_mspe = 3;
    int p = aKernel.GetVariablesNumber();
    int z_miss_number, n_z_obs;
    T *z_actual;
    if (!apTrainLocations && !apTestLocations) {
        z_miss_number = aConfigurations.GetUnknownObservationsNb();
        n_z_obs = aConfigurations.CalculateZObsNumber();
        z_actual = new T[z_miss_number * p];
    } else {
        z_miss_number = apTestLocations->GetSize();
        n_z_obs = apTrainLocations->GetSize();
        z_actual = nullptr;
    }
    aConfigurations.SetObservationNumber(n_z_obs);
    auto linear_algebra_solver = linearAlgebra::LinearAlgebraFactory<T>::CreateLinearAlgebraSolver(common::EXACT_DENSE);

    VERBOSE("\t- Total number of Z: " << aConfigurations.GetProblemSize())
    LOGGER("\t- Number of Z Miss: " << z_miss_number)
    LOGGER("\t- Number of Z observations: " << n_z_obs)

    // FISHER Prediction Function Call
    if (aConfigurations.GetIsFisher()) {
        LOGGER("\t---- Using Prediction Function Fisher ----")
        T *fisher_results;
        fisher_results = linear_algebra_solver->ExaGeoStatFisherTile(aConfigurations, aData,
                                                                     (T *) aConfigurations.GetEstimatedTheta().data(),
                                                                     aKernel);
        vector<double> fisher_vector;
        fisher_vector.reserve(num_params * num_params); // Reserve memory in advance for efficiency
        for (size_t idx = 0; idx < num_params * num_params; ++idx) {
            fisher_vector.push_back(fisher_results[idx]);
        }
        Results::GetInstance()->SetFisherMatrix(fisher_vector);
        VERBOSE("\t\t- Sd of sigma2, alpha, nu: " << sqrt(fisher_results[0]) << " " << sqrt(fisher_results[4]) << " "
                                                  << sqrt(fisher_results[8]))
        VERBOSE("\t\t- CI for sigma2: " << aConfigurations.GetEstimatedTheta()[0] - Q_NORM * sqrt(fisher_results[0])
                                        << " "
                                        << aConfigurations.GetEstimatedTheta()[0] + Q_NORM * sqrt(fisher_results[0]))

        VERBOSE("\t\t- CI for alpha: " << aConfigurations.GetEstimatedTheta()[1] - Q_NORM * sqrt(fisher_results[4])
                                       << " "
                                       << aConfigurations.GetEstimatedTheta()[1] + Q_NORM * sqrt(fisher_results[4]))
        VERBOSE("\t\t- CI for nu: " << aConfigurations.GetEstimatedTheta()[2] - Q_NORM * sqrt(fisher_results[8]) << " "
                                    << aConfigurations.GetEstimatedTheta()[2] + Q_NORM * sqrt(fisher_results[8]))

        VERBOSE("\t\t- Fisher Matrix:")
        for (i = 0; i < num_params; i++) {
            VERBOSE("\t\t ", true)
            for (j = 0; j < num_params; j++) {
                VERBOSE(fisher_results[i * num_params + j],true)
                if (j != num_params - 1) {
                    VERBOSE(", ",true)
                }
            }
            VERBOSE("")
        }
        delete[] fisher_results;
    }

    if (z_miss_number <= 0) {
        return;
    }

    Results::GetInstance()->SetZMiss(z_miss_number);
    aConfigurations.SetUnknownObservationsNb(z_miss_number);
    T *z_obs = new T[n_z_obs * p];
    T *z_miss = new T[z_miss_number];

    vector<T> avg_pred_value(number_of_mspe);
    auto miss_locations = new Locations<T>(z_miss_number, aConfigurations.GetDimension());
    // Prediction is only supported with 2D.
    auto obs_locations = new Locations<T>(n_z_obs, aConfigurations.GetDimension());

    // We Predict date with only Exact computation. This is a pre-request.
    InitializePredictionArguments(aConfigurations, aData, linear_algebra_solver, z_obs, z_actual, *miss_locations,
                                  *obs_locations, apMeasurementsMatrix, p, apTrainLocations, apTestLocations);

    // MLOE MMOM Auxiliary Function Call
    if (aConfigurations.GetIsMLOEMMOM()) {
        LOGGER("---- Using Auxiliary Function MLOE MMOM ----")
        linear_algebra_solver->ExaGeoStatMLETileMLOEMMOM(aConfigurations, aData,
                                                         (T *) aConfigurations.GetInitialTheta().data(),
                                                         (T *) aConfigurations.GetEstimatedTheta().data(),
                                                         *miss_locations, *obs_locations, aKernel);
    }

    // IDW Auxiliary Function Call
    if (aConfigurations.GetIsIDW()) {
        LOGGER("\t---- Using Auxiliary Function IDW ----")
        T *mspe = new T[number_of_mspe];

        if (!z_actual) {
            z_actual = new T[z_miss_number * p];
            memcpy(z_actual, apMeasurementsMatrix + n_z_obs, z_miss_number * sizeof(T));
        }

        PredictionAuxiliaryFunctions<T>::PredictIDW(z_miss, z_actual, z_obs, z_miss_number, n_z_obs, *miss_locations,
                                                    *obs_locations, mspe);

        vector<double> idw_error;
        idw_error.reserve(number_of_mspe); // Reserve memory in advance for efficiency
        for (size_t idx = 0; idx < number_of_mspe; ++idx) {
            idw_error.push_back(mspe[idx]);
        }

        Results::GetInstance()->SetIDWError(idw_error);
        delete[] mspe;
    }

    // MSPE Prediction Function Call
    if (aConfigurations.GetIsMSPE()) {
        LOGGER("\t---- Using Prediction Function MSPE ----")
        T *prediction_error_mspe;
        if (aConfigurations.GetIsNonGaussian()) {
            prediction_error_mspe = linear_algebra_solver->ExaGeoStatMLENonGaussianPredictTile(aData,
                                                                                               (T *) aConfigurations.GetEstimatedTheta().data(),
                                                                                               z_miss_number, n_z_obs,
                                                                                               z_obs, z_actual, z_miss,
                                                                                               aConfigurations,
                                                                                               *miss_locations,
                                                                                               *obs_locations, aKernel);
        } else {
            prediction_error_mspe = linear_algebra_solver->ExaGeoStatMLEPredictTile(aData,
                                                                                    (T *) aConfigurations.GetEstimatedTheta().data(),
                                                                                    z_miss_number, n_z_obs, z_obs,
                                                                                    z_actual, z_miss, aConfigurations,
                                                                                    *miss_locations, *obs_locations,
                                                                                    aKernel);
        }
        for (i = 0; i < number_of_mspe; i++) {
            avg_pred_value[i] += prediction_error_mspe[i];
        }
        vector<double> z_miss_vector;
        z_miss_vector.reserve(z_miss_number); // Reserve memory in advance for efficiency
        for (size_t idx = 0; idx < z_miss_number; ++idx) {
            z_miss_vector.push_back(z_miss[idx]);
        }
        Results::GetInstance()->SetPredictedMissedValues(z_miss_vector);
        if (z_actual) {
            LOGGER("\t\t- Prediction value: " << avg_pred_value[0])
        }
        delete[] prediction_error_mspe;
    }


    // Due to a leak in Chameleon, exactly trsm We had to free the buffer manually.
#ifdef USE_MKL
    mkl_free_buffers();
#endif

    delete[] z_obs;
    delete[] z_miss;
    delete[] z_actual;
    delete miss_locations;
    delete obs_locations;
#endif
}

template<typename T>
void Prediction<T>::InitializePredictionArguments(Configurations &aConfigurations, unique_ptr<ExaGeoStatData<T>> &aData,
                                                  unique_ptr<linearAlgebra::LinearAlgebraMethods<T>> &aLinearAlgebraSolver,
                                                  T *apZObs, T *apZActual, Locations<T> &aMissLocation,
                                                  Locations<T> &aObsLocation, T *apMeasurementsMatrix, const int &aP,
                                                  Locations<T> *apTrainLocations, Locations<T> *apTestLocations) {
#if DEFAULT_RUNTIME
    int full_problem_size = aConfigurations.GetProblemSize() * aP;
    T *z = new T[full_problem_size];

    aLinearAlgebraSolver->ExaGeoStatGetZObs(aConfigurations, z, full_problem_size, *aData->GetDescriptorData(),
                                            apMeasurementsMatrix, aP);

    if (!apTrainLocations && !apTestLocations) {
        PredictionHelpers<T>::PickRandomPoints(aConfigurations, aData, apZObs, apZActual, z, aMissLocation,
                                               aObsLocation, aP);
    } else {
        for (int i = 0; i < apTrainLocations->GetSize(); ++i) {
            aObsLocation.GetLocationX()[i] = apTrainLocations->GetLocationX()[i];
            aObsLocation.GetLocationY()[i] = apTrainLocations->GetLocationY()[i];
        }
        for (int i = 0; i < apTestLocations->GetSize(); ++i) {
            aMissLocation.GetLocationX()[i] = apTestLocations->GetLocationX()[i];
            aMissLocation.GetLocationY()[i] = apTestLocations->GetLocationY()[i];
        }
        memcpy(apZObs, apMeasurementsMatrix, aObsLocation.GetSize() * sizeof(T));
    }
    delete[] z;
#endif
}