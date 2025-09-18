#ifndef BOOLARRAY_H 
#define BOOLARRAY_H

#include <stdint.h>
#include <stdlib.h>
#include <stdbool.h>
#include "returncodes.h"
#include "idx.c"
typedef struct {Idx n; bool *data;} BoolArray;

Status false_BoolArray_data(BoolArray *arr, Idx n) {
    if (n < 0) return OOB_ERROR;
    arr->n = n;
    arr->data = calloc(n*sizeof(bool));  // each index is a pair (row, col)
    if (arr->data == NULL) {
        return NULL_POINTER_ERROR;
    }
    return SUCCESS;
}

Status set_BoolArray_to_zero(BoolArray *arr) {
    if (arr->data == NULL) return NULL_POINTER_ERROR;
    for (Idx i = 0; i < arr->n; i++) {
        arr->data[i] = false;
    }
    return SUCCESS;
}

Status get_BoolArray_bc(BoolArray arr, bool *out, Idx i) {
    if (i < 0 || i >= arr.n) return OOB_ERROR;
    if (out == NULL) return NULL_POINTER_ERROR;
    *out = arr.data[i];
    return SUCCESS;
}
Status set_BoolArray_bc(BoolArray *arr, Idx i, bool val) {
    if (i < 0 || i >= arr->n) return OOB_ERROR;
    if (arr->data == NULL) return NULL_POINTER_ERROR;
    arr->data[i] = val;
    return SUCCESS;
}

void free_BoolArray_data(BoolArray *arr) {
    if (arr->data != NULL) {
        free(arr->data);
        arr->data = NULL;
    }
}

#endif // BOOLARRAY_H