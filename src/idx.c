#ifndef IDX_C
#define IDX_C

#include <stdint.h>
#include <stdlib.h>
#include "returncodes.h"

typedef uint16_t Idx; // row or column index
typedef struct {Idx row; Idx col;} IdxPair;
typedef struct {Idx n; IdxPair *pairs;} IdxArray;
Status empty_IdxArray_data(IdxArray *arr, Idx n) {
    if (n < 0) return OOB_ERROR;
    arr->n = n;
    arr->pairs = malloc(n*sizeof(IdxPair));  // each index is a pair (row, col)
    if (arr->pairs == NULL) {
        return NULL_POINTER_ERROR;
    }
    return SUCCESS;
}
void free_IdxArray_data(IdxArray *arr) {
    if (arr->pairs != NULL) {
        free(arr->pairs);
        arr->pairs = NULL;
    }
}

Status get_IdxArray_bc(IdxArray arr, IdxPair out, Idx i) {
    if (i < 0 || i >= arr.n) return OOB_ERROR;
    out.row = arr.pairs[i].row;  // row
    out.col = arr.pairs[i].col;  // col
    return SUCCESS;
}
Status set_IdxArray_bc(IdxArray *arr, Idx i, Idx row, Idx col) {
    if (i < 0 || i >= arr->n) return OOB_ERROR;
    if (arr->pairs == NULL) return NULL_POINTER_ERROR;
    arr->pairs[i].row = row;
    arr->pairs[i].col = col;
    return SUCCESS;
}

#endif // IDX_C