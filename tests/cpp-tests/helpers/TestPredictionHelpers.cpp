
// Copyright (c) 2017-2023 King Abdullah University of Science and Technology,
// All rights reserved.
// ExaGeoStat is a software package, provided by King Abdullah University of Science and Technology (KAUST).

/**
 * @file TestPrediction.cpp
 * @brief Unit tests for the TestPrediction class in the ExaGeoStat software package.
 * @details This file contains Catch2 unit tests that validate the functionality of the TestPrediction class
 * @version 1.0.0
 * @author Mahmoud ElKarargy
 * @date 2023-12-08
**/

#include <catch2/catch_all.hpp>

#include <prediction/PredictionHelpers.hpp>
#include <prediction/PredictionAuxiliaryFunctions.hpp>

using namespace exageostat::prediction;
using namespace exageostat::dataunits;
using namespace exageostat::common;

void TEST_SHUFFLE_HELPER_FUNCTIONS() {
    SECTION("Test Shuffle with one array") {
        int size = 5;
        Dimension dimension = Dimension2D;

        //Create array and Locations to be shuffled.
        double array[] = {0, 1, 2, 3, 4};
        double expected_array[] = {0, 1, 2, 3, 4};

        Locations<double> locations(size, dimension);
        for (int i = 0; i < size; i++) {
            locations.GetLocationX()[i] = i;
            locations.GetLocationY()[i] = i;
        }

        PredictionHelpers<double>::Shuffle(array, locations, size);

        //Checks if the shuffle function changed the size.
        REQUIRE(sizeof(array) / sizeof(double) == size);
        REQUIRE(locations.GetSize() == size);

        //Check that the array is indeed shuffled.
        std::vector<int> original_array_elements(expected_array, expected_array + size);
        std::vector<int> shuffled_X_elements(locations.GetLocationX(), locations.GetLocationX() + size);
        std::vector<int> shuffled_Y_elements(locations.GetLocationY(), locations.GetLocationY() + size);
        std::vector<int> shuffled_array_elements(array, array + size);
        REQUIRE(original_array_elements != shuffled_array_elements);
        REQUIRE(shuffled_X_elements != original_array_elements);
        REQUIRE(shuffled_Y_elements != original_array_elements);

        //Check that even shuffled the array still has its original elements.
        std::sort(original_array_elements.begin(), original_array_elements.end());
        std::sort(shuffled_array_elements.begin(), shuffled_array_elements.end());
        std::sort(shuffled_X_elements.begin(), shuffled_X_elements.end());
        std::sort(shuffled_Y_elements.begin(), shuffled_Y_elements.end());
        REQUIRE(original_array_elements == shuffled_array_elements);
        REQUIRE(original_array_elements == shuffled_X_elements);
        REQUIRE(original_array_elements == shuffled_Y_elements);

    }SECTION("Test Shuffle with two arrays") {
        int size = 5;
        Dimension dimension = Dimension2D;
        //Create arrays and Locations to be shuffled.
        double array_1[] = {0, 1, 2, 3, 4};
        double expected_array_1[] = {0, 1, 2, 3, 4};

        double array_2[] = {5, 6, 7, 8, 9};
        double expected_array_2[] = {5, 6, 7, 8, 9};

        Locations<double> locations(size, dimension);
        for (int i = 0; i < size; i++) {
            locations.GetLocationX()[i] = i;
            locations.GetLocationY()[i] = i;
        }

        PredictionHelpers<double>::Shuffle(array_1, array_2, locations, size);

        //Checks if the shuffle function changed the size.
        REQUIRE(sizeof(array_1) / sizeof(double) == size);
        REQUIRE(sizeof(array_2) / sizeof(double) == size);
        REQUIRE(locations.GetSize() == size);

        //Check that the array is indeed shuffled.
        std::vector<int> original_array_elements_1(expected_array_1, expected_array_1 + size);
        std::vector<int> original_array_elements_2(expected_array_2, expected_array_2 + size);

        std::vector<int> shuffled_X_elements(locations.GetLocationX(), locations.GetLocationX() + size);
        std::vector<int> shuffled_Y_elements(locations.GetLocationY(), locations.GetLocationY() + size);

        std::vector<int> shuffled_array_elements_1(array_1, array_1 + size);
        std::vector<int> shuffled_array_elements_2(array_2, array_2 + size);

        REQUIRE(original_array_elements_1 != shuffled_array_elements_1);
        REQUIRE(original_array_elements_2 != shuffled_array_elements_2);

        REQUIRE(shuffled_X_elements != original_array_elements_1);
        REQUIRE(shuffled_Y_elements != original_array_elements_1);

        //Check that even shuffled the array_1 still has its original elements.
        std::sort(original_array_elements_1.begin(), original_array_elements_1.end());
        std::sort(shuffled_array_elements_1.begin(), shuffled_array_elements_1.end());
        std::sort(original_array_elements_2.begin(), original_array_elements_2.end());
        std::sort(shuffled_array_elements_2.begin(), shuffled_array_elements_2.end());
        std::sort(shuffled_X_elements.begin(), shuffled_X_elements.end());
        std::sort(shuffled_Y_elements.begin(), shuffled_Y_elements.end());

        REQUIRE(original_array_elements_1 == shuffled_array_elements_1);
        REQUIRE(original_array_elements_2 == shuffled_array_elements_2);
        REQUIRE(original_array_elements_1 == shuffled_X_elements);
        REQUIRE(original_array_elements_1 == shuffled_Y_elements);

    }SECTION("Test Shuffle with three arrays") {
        int size = 5;
        Dimension dimension = Dimension2D;

        //Create arrays and Locations to be shuffled.
        double array_1[] = {0, 1, 2, 3, 4};
        double expected_array_1[] = {0, 1, 2, 3, 4};

        double array_2[] = {5, 6, 7, 8, 9};
        double expected_array_2[] = {5, 6, 7, 8, 9};

        double array_3[] = {10, 11, 12, 13, 14};
        double expected_array_3[] = {10, 11, 12, 13, 14};

        Locations<double> locations(size, dimension);
        for (int i = 0; i < size; i++) {
            locations.GetLocationX()[i] = i;
            locations.GetLocationY()[i] = i;
        }

        PredictionHelpers<double>::Shuffle(array_1, array_2, array_3, locations, size);

        //Checks if the shuffle function changed the size.
        REQUIRE(sizeof(array_1) / sizeof(double) == size);
        REQUIRE(sizeof(array_2) / sizeof(double) == size);
        REQUIRE(sizeof(array_3) / sizeof(double) == size);

        REQUIRE(locations.GetSize() == size);

        //Check that the array is indeed shuffled.
        std::vector<int> original_array_elements_1(expected_array_1, expected_array_1 + size);
        std::vector<int> original_array_elements_2(expected_array_2, expected_array_2 + size);
        std::vector<int> original_array_elements_3(expected_array_3, expected_array_3 + size);

        std::vector<int> shuffled_X_elements(locations.GetLocationX(), locations.GetLocationX() + size);
        std::vector<int> shuffled_Y_elements(locations.GetLocationY(), locations.GetLocationY() + size);

        std::vector<int> shuffled_array_elements_1(array_1, array_1 + size);
        std::vector<int> shuffled_array_elements_2(array_2, array_2 + size);
        std::vector<int> shuffled_array_elements_3(array_3, array_3 + size);

        REQUIRE(original_array_elements_1 != shuffled_array_elements_1);
        REQUIRE(original_array_elements_2 != shuffled_array_elements_2);
        REQUIRE(original_array_elements_3 != shuffled_array_elements_3);

        REQUIRE(shuffled_X_elements != original_array_elements_1);
        REQUIRE(shuffled_Y_elements != original_array_elements_1);

        //Check that even shuffled, the array still has its original elements.
        std::sort(original_array_elements_1.begin(), original_array_elements_1.end());
        std::sort(shuffled_array_elements_1.begin(), shuffled_array_elements_1.end());

        std::sort(original_array_elements_2.begin(), original_array_elements_2.end());
        std::sort(shuffled_array_elements_2.begin(), shuffled_array_elements_2.end());

        std::sort(original_array_elements_3.begin(), original_array_elements_3.end());
        std::sort(shuffled_array_elements_3.begin(), shuffled_array_elements_3.end());

        std::sort(shuffled_X_elements.begin(), shuffled_X_elements.end());
        std::sort(shuffled_Y_elements.begin(), shuffled_Y_elements.end());

        REQUIRE(original_array_elements_1 == shuffled_array_elements_1);
        REQUIRE(original_array_elements_2 == shuffled_array_elements_2);
        REQUIRE(original_array_elements_3 == shuffled_array_elements_3);

        REQUIRE(original_array_elements_1 == shuffled_X_elements);
        REQUIRE(original_array_elements_1 == shuffled_Y_elements);

    }
}

void TEST_SORT_HELPER_FUNCTION() {

    SECTION("Test SortArray - 1D") {
        const int size = 8;
        //Create array and Locations to be sorted.
        uint32_t data[size] = {170, 45, 75, 90, 802, 24, 2, 66};
        uint32_t expected_sorted_data[size] = {2, 24, 45, 66, 75, 90, 170, 802};

        // Sort the data using RadixSort
        PredictionHelpers<double>::SortArray(data, size);

        // Verify that the data is sorted correctly
        for (int i = 0; i < size; i++) {
            REQUIRE(data[i] == expected_sorted_data[i]);
        }
    }

    SECTION("Test SortArray for duplicated values - 1D") {
        const int size = 10;
        //Create array and Locations to be sorted.
        uint32_t data[size] = {5, 4, 3, 2, 1, 5, 4, 3, 3, 3};
        uint32_t expected_sorted_data[size] = {1, 2, 3, 3, 3, 3, 4, 4, 5, 5};

        // Sort the data using RadixSort
        PredictionHelpers<double>::SortArray(data, size);

        // Verify that the data is sorted correctly
        for (int i = 0; i < size; i++) {
            REQUIRE(data[i] == expected_sorted_data[i]);
        }
    }


    SECTION("Test SortArray - 2D") {
        const int size = 10;
        const int dim = 2;
        uint32_t data[size][dim] = {
                {170, 75},
                {802, 2},
                {60,  34},
                {74,  15},
                {74,  15},
                {15,  74},
                {45,  90},
                {24,  66},
                {22,  33},
                {45,  5}
        };
        uint32_t expected_sorted_data[
                size * dim] = {2, 5, 15, 15, 15, 22, 24, 33, 34, 45, 45, 60, 66, 74, 74, 74, 75, 90, 170, 802};

        // Sort the 2D array using 'RadixSort'
        PredictionHelpers<double>::SortArray(*data, size * dim);

        for (int i = 0; i < size * dim; i++) {
            REQUIRE(*(*data + i) == expected_sorted_data[i]);
        }
    }
}

TEST_CASE("Prediction tests") {
    TEST_SHUFFLE_HELPER_FUNCTIONS();

    TEST_SORT_HELPER_FUNCTION();
}

