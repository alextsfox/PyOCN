#include <stdint.h>
#include <stdlib.h>
#include "returncodes.h"
#include "neighbormatrix.c"

/**
 * @brief computes the drained area of each vertex in the neighbor matrix A, given a root vertex.
 * @param A The neighbor matrix
 * @param root The root vertex (e.g. the outlet)
 * @return The drained area of each vertex in the neighbor matrix A.
 */

/*
Thinking out loud:
* we don't want to have to visit each vertex multiple times
* each time a vertex is visited, we 
*/
int32_t compute_drained_area(NeighborMat A, Vertex root, int32_t *drained_area) {
    return 0;
}