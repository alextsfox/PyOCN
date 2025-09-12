#ifndef ARRAYMATH_H
#define ARRAYMATH_H

#include "array.h"

/**
 * @file arraymath.h
 * @brief Element-wise math operations for Arr2D arrays.
 */

/**
 * @brief Element-wise addition: result = a + b
 * @param a First input array
 * @param b Second input array
 * @param result Output array (must be allocated, same shape as a and b)
 */
void add_Arr2D(Arr2D a, Arr2D b, Arr2D *result);

/**
 * @brief Element-wise subtraction: result = a - b
 */
void sub_Arr2D(Arr2D a, Arr2D b, Arr2D *result);

/**
 * @brief Element-wise multiplication: result = a * b
 */
void mul_Arr2D(Arr2D a, Arr2D b, Arr2D *result);

/**
 * @brief Element-wise division: result = a / b (asserts b.data[i] != 0)
 */
void div_Arr2D(Arr2D a, Arr2D b, Arr2D *result);

/**
 * @brief Element-wise power: result = pow(a, b)
 */
void pow_Arr2D(Arr2D a, Arr2D b, Arr2D *result);

/**
 * @brief Element-wise negation: result = -a
 */
void negate_Arr2D(Arr2D a, Arr2D *result);

/**
 * @brief Element-wise addition with scalar: result = a + scalar
 */
void scalar_add_Arr2D(Arr2D a, double scalar, Arr2D *result);

/**
 * @brief Element-wise subtraction with scalar: result = a - scalar
 */
void scalar_sub_Arr2D(Arr2D a, double scalar, Arr2D *result);

/**
 * @brief Element-wise multiplication with scalar: result = a * scalar
 */
void scalar_mul_Arr2D(Arr2D a, double scalar, Arr2D *result);

/**
 * @brief Element-wise division by scalar: result = a / scalar (asserts scalar != 0)
 */
void scalar_div_Arr2D(Arr2D a, double scalar, Arr2D *result);

/**
 * @brief Element-wise power with scalar: result = pow(a, scalar)
 */
void scalar_pow_Arr2D(Arr2D a, double scalar, Arr2D *result);

/**
 * @brief Element-wise exponentiation: result = pow(scalar, b)
 */
void scalar_exp_Arr2D(double scalar, Arr2D b, Arr2D *result);

#endif // ARRAYMATH_H
