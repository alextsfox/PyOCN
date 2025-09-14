#ifndef NEIGHBORMATRIX_C
#define NEIGHBORMATRIX_C

#include <stdint.h>
#include <stdlib.h>
#include <stdio.h>
#include "returncodes.h"

#define edge_exists(vert, hand) (vert.edges & (1 << hand))
#define downstream_edge_exists(vert, hand) (edge_exists(vert, hand) & (hand == vert.downstream))
#define upstream_edge_exists(vert, hand) (edge_exists(vert, hand) & (hand != vert.downstream))

/**
 * @brief Represents a row or column index in the neighbor matrix.
 */
typedef int16_t Idx;

/**
 * @brief Represents the presence or absence of all edges connecting to a vertex as an 8-bit unsigned int.
 */
typedef uint8_t LocalEdges;

/**
 * @brief Represents a cardinal direction (0-7). 0 = top, 1 = top right, 2 = right, etc.
 */
typedef uint8_t Clockhand; // 0-7. 0 = top, 1 = top right, 2 = right, etc.

/**
 * @brief Represents a vertex in the neighbor matrix.
 * @param edges 8 bits representing edges to 8 neighbors
 * @param downstream The downstream direction
 */
typedef struct{
    LocalEdges edges; // 8 bits representing edges to 8 neighbors
    Clockhand downstream;
} Vertex;

/**
 * @brief Neighbor matrix structure.
 * 
 * The neighbor matrix is a 2D array where each element indicates the presence
 * (1) or absence (0) of an edge between vertices in a graph. Each element N[i,j] 
 * represents the edges connecting the vertex at row i, column j to its
 * 8 immediate neighbors as an 8 bit unsigned int. The least significant bit 
 * represents the topmost neighbor, and the bits proceed clockwise from there.
 * 
 * Example
 * -------
 * Consider the following spanning tree of a 3x3 matrix
 * o o o
 *  \|/
 * o-o-o
 *    /|
 * o-o o
 * 
 * Vertex 1,1 (in the middle) is connected to neighbors 0 (top), 1 (top right),
 * 2 (right), 6 (left), and 7 (top left), so A[1,1] = 0b11000111 = 0xc7
 * 
 * Additionally, NeighborMat contains a "downstream" array. For each element, it contains the downstream direction.
 */
typedef struct {
    Idx m;
    Idx n;
    Vertex *vertices;
} NeighborMat;

/**
 * @brief Retrieves the vertex stored at (i,j) in the NeighborMat A and stores it in vert. Performs bounds checking.
 * @param A The neighbor matrix
 * @param vert Pointer to store the retrieved vertex
 * @param i Row index
 * @param j Column index
 * @return SUCCESS or an error code
 */
Status get_NeighborMat_bc(NeighborMat A, Vertex *vert, Idx i, Idx j) {
    if (i < 0 || i >= A.m || j < 0 || j >= A.n) return OOB_ERROR;
    if (vert == NULL) return NULL_POINTER_ERROR;
    *vert = A.vertices[i * A.n + j];
    return SUCCESS;
}

/**
 * @brief Sets the vertex stored at (i,j) in the NeighborMat A to the given value. Performs bounds checking.
 * @param A Pointer to the neighbor matrix
 * @param i Row index
 * @param j Column index
 * @param vert The new vertex value to set
 * @return SUCCESS or an error code
 */
Status set_NeighborMat_bc(NeighborMat *A, Idx i, Idx j, Vertex vert) {
    if (i < 0 || i >= A->m || j < 0 || j >= A->n) return OOB_ERROR;
    if (A->vertices == NULL) return NULL_POINTER_ERROR;
    A->vertices[i * A->n + j] = vert;
    return SUCCESS;
}

static const Idx row_neighbors[8] = {-1, -1, 0, 1, 1, 1, 0, -1};
static const Idx col_neighbors[8] = {0, 1, 1, 1, 0, -1, -1, -1};

/**
 * @brief Get the row and column indices of neighboring vertices for a given cell.
 * @param A The adjacency matrix
 * @param out Array of size 2 to store the row and column indices of the neighbor
 * @param i Row index of the vertex in question
 * @param j Column index of the vertex in question
 * @return SUCCESS or an error code
 */
Status get_clockhand_neighbor(NeighborMat A, Idx out[2], Idx i, Idx j, Clockhand hand) {
    // returns the index of the neighbor pointed to by the clockhand (0-7) of the vertex at (i,j).
    Vertex vert;
    Status code = get_NeighborMat_bc(A, &vert, i, j);
    if (code != SUCCESS) return code;
    if (edge_exists(vert, hand) == 0) return NO_EDGE_WARNING; // no neighbor in that direction

    out[0] = i + row_neighbors[hand];
    out[1] = j + col_neighbors[hand];
    return SUCCESS;
}

/**
 * @brief Move the outgoing edge from direction handa to handb for the vertex at (i,j) in the neighbor matrix A.
 * @param A Pointer to the neighbor matrix
 * @param i Row index of the vertex
 * @param j Column index of the vertex
 * @param newdownstream The new downstream direction
 * @return SUCCESS or an error code. 
 */
/*
1. Extract the source vertex at (i,j).
2. Extract the 
*/
Status swap_edges(NeighborMat *A, Idx i, Idx j, Clockhand newdownstream) {
    NeighborMat Acopy = *A ; // avoid dereferencing A multiple times. Used only for lookup, not modification.
    Vertex sourcenode, downnode, newdownnode;
    Status code;
    
    // confirm that the source node exists
    code = get_NeighborMat_bc(Acopy, &sourcenode, i, j);
    if (code != SUCCESS) return code;
    if (newdownstream == sourcenode.downstream) return SUCCESS; // nothing to do
    if (edge_exists(sourcenode, newdownstream)) return SWAP_WARNING; // can't swap to an upstream edge
    // So long as the initial array was constructed correctly, this should not be necessary.
    // could remove later.
    if (!downstream_edge_exists(sourcenode, sourcenode.downstream)) return NO_EDGE_WARNING;  

    LocalEdges sourcemask = (1 << sourcenode.downstream) | (1 << newdownstream);
    LocalEdges downmask = (1 << sourcenode.downstream) + 4;
    LocalEdges newdownmask = (1 << newdownstream) + 4;

    // we don't use get/set_NeighborMat_bc here to avoid redundant bounds checking
    Idx downidx[2], newdownidx[2];
    code = get_clockhand_neighbor(Acopy, downidx, i, j, sourcenode.downstream);
    if (code != SUCCESS) return code;
    downnode = Acopy.vertices[downidx[0] * Acopy.n + downidx[1]];
    code = get_clockhand_neighbor(Acopy, newdownidx, i, j, newdownstream);
    if (code != SUCCESS) return code;
    newdownnode = Acopy.vertices[newdownidx[0] * Acopy.n + newdownidx[1]];

    A->vertices[i * Acopy.n + j] = (Vertex){sourcenode.edges ^ sourcemask, newdownstream};
    A->vertices[downidx[0] * Acopy.n + downidx[1]] = (Vertex){downnode.edges ^ downmask, downnode.downstream};
    A->vertices[newdownidx[0] * Acopy.n + newdownidx[1]] = (Vertex){newdownnode.edges ^ newdownmask, newdownnode.downstream};

    return SUCCESS;
}



void print_vert(Vertex vert){
    
    // top row of Os
    printf("O O O\n");

    // top row of connections
    printf(" ");
    if (edge_exists(vert, 7)) printf("\\"); else printf(" ");
    if (edge_exists(vert, 0)) printf("|"); else printf(" ");
    if (edge_exists(vert, 1)) printf("/"); else printf(" ");
    printf("\n");
    // middle row of Os
    printf("O");
    if (edge_exists(vert, 6)) printf("-"); else printf(" ");
    printf("O");
    if (edge_exists(vert, 2)) printf("-"); else printf(" ");
    printf("O");
    printf("\n");
    
    // bottom row of connections
    printf(" ");
    if (edge_exists(vert, 5)) printf("/"); else printf(" ");
    if (edge_exists(vert, 4)) printf("|"); else printf(" ");
    if (edge_exists(vert, 3)) printf("\\"); else printf(" ");

    // bottom row of O's
    printf("O O O\n");
}

void free_NeighborMat(NeighborMat *A) {
    if (A->vertices != NULL) {
        free(A->vertices);
        A->vertices = NULL;
    }
}

#endif // NEIGHBORMATRIX_C