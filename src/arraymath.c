// #include <stdint.h>
// #include <math.h>
// #include <stdlib.h>
// #include <assert.h>
// #include "array.h"

// static void assert_binary_op_compatible(Arr2D a, Arr2D b, Arr2D *result) {
//     assert(
//         a.rows == b.rows 
//         && a.cols == b.cols 
//         && result->rows == a.rows 
//         && result->cols == a.cols
//         && a.data != NULL
//         && b.data != NULL
//         && result->data != NULL
//     );
// }


// // Macros for element-wise operations
// #define ARR2D_BINARY_OP(name, op) \
// void name(Arr2D a, Arr2D b, Arr2D *result) { \
//     assert_binary_op_compatible(a, b, result); \
//     for (int32_t i = 0; i < a.rows * a.cols; i++) { \
//         result->data[i] = a.data[i] op b.data[i]; \
//     } \
// }

// #define ARR2D_UNARY_OP(name, expr) \
// void name(Arr2D a, Arr2D *result) { \
//     assert_binary_op_compatible(a, a, result); \
//     for (int32_t i = 0; i < a.rows * a.cols; i++) { \
//         result->data[i] = (expr); \
//     } \
// }

// #define ARR2D_SCALAR_OP(name, op) \
// void name(Arr2D a, double scalar, Arr2D *result) { \
//     assert_binary_op_compatible(a, a, result); \
//     for (int32_t i = 0; i < a.rows * a.cols; i++) { \
//         result->data[i] = a.data[i] op scalar; \
//     } \
// }

// // Binary operations //
// ARR2D_BINARY_OP(add_Arr2D, +)
// ARR2D_BINARY_OP(sub_Arr2D, -)
// ARR2D_BINARY_OP(mul_Arr2D, *)
// ARR2D_BINARY_OP(div_Arr2D, /)
// ARR2D_UNARY_OP(negate_Arr2D, -a.data[i])

// // power and division need special handling //
// void pow_Arr2D(Arr2D a, Arr2D b, Arr2D *result) {
//     assert_binary_op_compatible(a, b, result);
//     for (int32_t i = 0; i < a.rows * a.cols; i++) {
//         result->data[i] = pow(a.data[i], b.data[i]);
//     }
// }

// void div_Arr2D(Arr2D a, Arr2D b, Arr2D *result) {
//     assert_binary_op_compatible(a, b, result);
//     for (int32_t i = 0; i < a.rows * a.cols; i++) {
//         assert(b.data[i] != 0.0);
//         result->data[i] = a.data[i] / b.data[i];
//     }
// }

// ARR2D_SCALAR_OP(scalar_add_Arr2D, +)
// ARR2D_SCALAR_OP(scalar_sub_Arr2D, -)
// ARR2D_SCALAR_OP(scalar_mul_Arr2D, *)

// void scalar_div_Arr2D(Arr2D a, double scalar, Arr2D *result) {
//     assert(scalar != 0.0);
//     assert_binary_op_compatible(a, a, result);
//     for (int32_t i = 0; i < a.rows * a.cols; i++) {
//         result->data[i] = a.data[i] / scalar;
//     }
// }

// void scalar_pow_Arr2D(Arr2D a, double scalar, Arr2D *result) {
//     assert_binary_op_compatible(a, a, result);
//     for (int32_t i = 0; i < a.rows * a.cols; i++) {
//         result->data[i] = pow(a.data[i], scalar);
//     }
// }

// void scalar_exp_Arr2D(double a, Arr2D b, Arr2D *result) {
//     assert_binary_op_compatible(b, b, result);
//     for (int32_t i = 0; i < b.rows * b.cols; i++) {
//         result->data[i] = pow(a, b.data[i]);    
//     } 
// }

#include <stdint.h>
#include <math.h>
#include <stdlib.h>
#include <assert.h>
#include "array.h"
#include "arraymath.h"

static void assert_binary_op_compatible(Arr2D a, Arr2D b, Arr2D *result) {
    assert(
        a.rows == b.rows 
        && a.cols == b.cols 
        && result->rows == a.rows 
        && result->cols == a.cols
        && a.data != NULL
        && b.data != NULL
        && result->data != NULL
    );
}

void add_Arr2D(Arr2D a, Arr2D b, Arr2D *result) {
    assert_binary_op_compatible(a, b, result);
    for (int32_t i = 0; i < a.rows * a.cols; i++) {
        result->data[i] = a.data[i] + b.data[i];
    }
}
void add_Arr2D_by_ref(Arr2D *a, Arr2D *b, Arr2D *result) {
    assert_binary_op_compatible(*a, *b, result);
    for (int32_t i = 0; i < a->rows * a->cols; i++) {
        result->data[i] = a->data[i] + b->data[i];
    }
}
void sub_Arr2D(Arr2D a, Arr2D b, Arr2D *result) {
    assert_binary_op_compatible(a, b, result);
    for (int32_t i = 0; i < a.rows * a.cols; i++) {
        result->data[i] = a.data[i] - b.data[i];
    }
}
void negate_Arr2D(Arr2D a, Arr2D *result) {
    assert_binary_op_compatible(a, a, result);
    for (int32_t i = 0; i < a.rows * a.cols; i++) {
        result->data[i] = -a.data[i];
    }
}
void pow_Arr2D(Arr2D a, Arr2D b, Arr2D *result) {
    assert_binary_op_compatible(a, b, result);
    for (int32_t i = 0; i < a.rows * a.cols; i++) {
        result->data[i] = pow(a.data[i], b.data[i]);
    }
}
void mul_Arr2D(Arr2D a, Arr2D b, Arr2D *result) {
    assert_binary_op_compatible(a, b, result);
    for (int32_t i = 0; i < a.rows * a.cols; i++) {
        result->data[i] = a.data[i] * b.data[i];
    }
}
void div_Arr2D(Arr2D a, Arr2D b, Arr2D *result) {
    assert_binary_op_compatible(a, b, result);
    for (int32_t i = 0; i < a.rows * a.cols; i++) {
        assert(b.data[i] != 0.0);
        result->data[i] = a.data[i] / b.data[i];
    }
}

void scalar_add_Arr2D(Arr2D a, double scalar, Arr2D *result) {
    assert_binary_op_compatible(a, a, result);
    for (int32_t i = 0; i < a.rows * a.cols; i++) {
        result->data[i] = a.data[i] + scalar;
    }
}
void scalar_sub_Arr2D(Arr2D a, double scalar, Arr2D *result) {
    assert_binary_op_compatible(a, a, result);
    for (int32_t i = 0; i < a.rows * a.cols; i++) {
        result->data[i] = a.data[i] - scalar;
    }
}
void scalar_mul_Arr2D(Arr2D a, double scalar, Arr2D *result) {
    assert_binary_op_compatible(a, a, result);
    for (int32_t i = 0; i < a.rows * a.cols; i++) {
        result->data[i] = a.data[i] * scalar;
    }
}
void scalar_div_Arr2D(Arr2D a, double scalar, Arr2D *result) {
    assert(scalar != 0.0);
    assert_binary_op_compatible(a, a, result);
    for (int32_t i = 0; i < a.rows * a.cols; i++) {
        result->data[i] = a.data[i] / scalar;
    }
}
void scalar_pow_Arr2D(Arr2D a, double scalar, Arr2D *result) {
    assert_binary_op_compatible(a, a, result);
    for (int32_t i = 0; i < a.rows * a.cols; i++) {
        result->data[i] = pow(a.data[i], scalar);
    }
}
void scalar_exp_Arr2D(double a, Arr2D b, Arr2D *result) {
    assert_binary_op_compatible(b, b, result);
    for (int32_t i = 0; i < b.rows * b.cols; i++) {
        result->data[i] = pow(a, b.data[i]);    
    } 
}


