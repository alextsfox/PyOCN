#include <stdint.h>
#include <math.h>
#include <stdlib.h>
#include <assert.h>
#include <string.h>
#include "array.h"

static void assert_inbounds(Arr2D arr, int32_t i, int32_t j){
    assert(
        arr.data != NULL 
        && i >= 0 
        && j >= 0
        && i < arr.rows 
        && j < arr.cols 
    );
}

Arr2D empty_Arr2D(int32_t rows, int32_t cols) {
    Arr2D arr = {
        .rows = rows,
        .cols = cols,
        .data = (double *)malloc(rows * cols * sizeof(double))
    };
    assert(arr.data != NULL && rows > 0 && cols > 0);
    return arr;
}

Arr2D zeros_Arr2D(int32_t rows, int32_t cols) {
    Arr2D arr = {
        .rows = rows,
        .cols = cols,
        .data = (double *)calloc(rows * cols, sizeof(double))
    };
    assert(arr.data != NULL && rows > 0 && cols > 0);
    return arr;
}

Arr2D get_slice_Arr2D(Arr2D arr, int32_t i0, int32_t i1, int32_t j0, int32_t j1) {
    assert_inbounds(arr, i0, j0);
    assert_inbounds(arr, i1 - 1, j1 - 1);

    Arr2D subarr = {
        .rows = i1 - i0,
        .cols = j1 - j0,
        .data = arr.data + (i0 * arr.cols) + j0
    };
    return subarr;
}

Arr2D copy_Arr2D(Arr2D arr) {
    Arr2D copy = empty_Arr2D(arr.rows, arr.cols);
    memcpy(copy.data, arr.data, arr.rows * arr.cols * sizeof(double));
    return copy;
}

double get_Arr2D(Arr2D arr, int32_t i, int32_t j) {
    assert_inbounds(arr, i, j);
    return arr.data[i * arr.cols + j];
}

void set_Arr2D(Arr2D *arr, int32_t i, int32_t j, double value) {
    assert_inbounds(*arr, i, j);
    arr->data[i * arr->cols + j] = value;
}

void free_Arr2D(Arr2D *arr) {
    if (arr->data != NULL) free(arr->data);
    arr->data = NULL;
    arr->rows = 0;
    arr->cols = 0;
}

