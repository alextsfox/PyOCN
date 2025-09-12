#ifndef ARRAY_H
#define ARRAY_H

#include <stdint.h>

/**
 * @file array.h
 * @brief Definition and basic operations for 2D arrays (Arr2D).
 */

/**
 * @struct Arr2D
 * @brief A simple 2D array structure.
 * @param rows Number of rows
 * @param cols Number of columns
 * @param data Pointer to the array data (flattened in row-major order)
 */
typedef struct {
    int32_t rows;
    int32_t cols;
    double *data;
} Arr2D;

/**
 * @brief Creates an empty 2D array with the specified dimensions.
 * @param rows Number of rows
 * @param cols Number of columns
 * @return Arr2D struct with allocated data (uninitialized)
 */
Arr2D empty_Arr2D(int32_t rows, int32_t cols);

/**
 * @brief Creates a 2D array filled with zeros.
 * @param rows Number of rows
 * @param cols Number of columns
 * @return Arr2D struct with allocated data (zero-initialized)
 */
Arr2D zeros_Arr2D(int32_t rows, int32_t cols);

/**
 * @brief Creates a subarray (slice) from an existing 2D array as a view into the original array.
 *        Ordered (increasing) inbounds ranges are required.
 * @param arr The original Arr2D array
 * @param i0 Starting row index for the slice (inclusive)
 * @param i1 Ending row index for the slice (exclusive)
 * @param j0 Starting column index for the slice (inclusive)
 * @param j1 Ending column index for the slice (exclusive)
 * @return Arr2D struct representing the subarray view
 */
Arr2D get_slice_Arr2D(Arr2D arr, int32_t i0, int32_t i1, int32_t j0, int32_t j1);

/**
 * @brief Creates a copy of a 2D array.
 * @param arr The array to copy
 * @return Arr2D struct with copied data
 */
Arr2D copy_Arr2D(Arr2D arr);

/**
 * @brief Gets the value at the specified row and column in a 2D array.
 * @param arr The array
 * @param i Row index
 * @param j Column index
 * @return The value at (i, j)
 */
double get_Arr2D(Arr2D arr, int32_t i, int32_t j);

/**
 * @brief Sets the value at the specified row and column in a 2D array.
 * @param arr Pointer to the array
 * @param i Row index
 * @param j Column index
 * @param value Value to set
 */
void set_Arr2D(Arr2D *arr, int32_t i, int32_t j, double value);

/**
 * @brief Frees the memory associated with a 2D array and resets its fields.
 * @param arr Pointer to the array to free
 */
void free_Arr2D(Arr2D *arr);

#endif