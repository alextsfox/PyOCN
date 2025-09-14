#include <stdint.h>
#include <stdlib.h>
#include <assert.h>
#include <stdbool.h>
#include "adjacencymatrix.h"

#define edge_exists(vert, hand) ((vert) & (1 << (hand)))
#define SUCCESS 256
#define NO_EDGE 257
#define OUT_OF_BOUNDS 258

typedef int16_t Idx; // for indexing rows and columns
typedef uint16_t Vertex; // 16-bit unsigned int representing edges to 8 neighbors. Values over 255 are used as return codes.
typedef uint8_t Clockhand; // 0-7

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
 */
typedef struct {
    Idx rows;
    Idx cols;
    Vertex *data;
} NeighborMat;

/**
 * @brief Retrieves the vertex stored at (i,j) in the NeighborMat A. Performs bounds checking.
 * @param A The neighbor matrix
 * @param i Row index
 * @param j Column index
 * @return The vertex at (i,j) or an error code
 */
Vertex get_NeighborMat_bc(NeighborMat A, Idx i, Idx j) {
    if (i < 0 || i >= A.rows || j < 0 || j >= A.cols) return OUT_OF_BOUNDS;
    return A.data[i * A.cols + j];
}

/**
 * @brief Sets the vertex stored at (i,j) in the NeighborMat A to the given value. Performs bounds checking.
 * @param A Pointer to the neighbor matrix
 * @param i Row index
 * @param j Column index
 * @param vert The new vertex value to set
 * @return SUCCESS or an error code
 */
Vertex set_NeighborMat_bc(NeighborMat *A, Idx i, Idx j, Vertex vert) {
    if (i < 0 || i >= A->rows || j < 0 || j >= A->cols) return OUT_OF_BOUNDS;
    A->data[i * A->cols + j] = vert;
    return SUCCESS;
}

static const Idx row_neighbors[8] = {-1, -1, 0, 1, 1, 1, 0, -1};
static const Idx col_neighbors[8] = {0, 1, 1, 1, 0, -1, -1, -1};

/**
 * @brief Get the row and column indices of neighboring vertices for a given cell.
 * @param A The adjacency matrix
 * @param i Row index of the cell
 * @param j Column index of the cell
 * @return A NeighborVertices struct containing the number of neighbors and their indices
 */
Vertex get_clockhand_neighbor(NeighborMat A, Idx i, Idx j, Clockhand hand, Idx out[2]) {
    // returns the index of the neighbor pointed to by the clockhand (0-7) of the vertex at (i,j).
    Vertex vert = get_NeighborMat_bc(A, i, j);
    if (vert == OUT_OF_BOUNDS) return OUT_OF_BOUNDS;
    if (!edge_exists(vert, hand)) return NO_EDGE; // no neighbor in that direction

    out[0] = i + row_neighbors[hand];
    out[1] = j + col_neighbors[hand];
    return SUCCESS;
}

/**
 * @brief swap edges in the adjacency matrix according to the given permutation.
 * @param A Pointer to the neighbor matrix
 * @param perm The permutation to apply
 * @return SUCCESS or an error code
 */
Vertex swap_edges(NeighborMat *A, Permutation perm) {
    // permutes a single edge in NeighborMat according to the permutation mapping perm.
   
    NeighborMat Acopy = *A;  // avoid dereferencing A multiple times. Used only for lookup, not modification.

    Vertex source = get_NeighborMat_bc(Acopy, perm.i, perm.j);
    if (source == OUT_OF_BOUNDS) return OUT_OF_BOUNDS;
    
    if (!edge_exists(source, perm.oldhand)) return NO_EDGE;
    if (edge_exists(source, perm.newhand)) return NO_EDGE;

    // 1. The source vertex (remove old edge to target, add new edge to new target)
    // 2. The old target vertex (remove old edge to source)
    // 3. The new target vertex (add new edge to source)

    // bitmasks for toggling edges. Will be XORed with existing vert.
    Vertex sourcemask = (1 << perm.oldhand) | (1 << perm.newhand);
    Vertex oldtgtmask = (1 << perm.oldhand) + 4;
    Vertex newtgtmask = (1 << perm.newhand) + 4;

    Idx oldtgtidx[2], newtgtidx[2];
    Vertex oldretval = get_clockhand_neighbor(Acopy, perm.i, perm.j, perm.oldhand, oldtgtidx);
    if (oldretval != SUCCESS) return oldretval;
    Vertex newretval = get_clockhand_neighbor(Acopy, perm.i, perm.j, perm.newhand, newtgtidx);
    if (newretval != SUCCESS) return newretval;

    Idx ioldtgt = oldtgtidx[0];
    Idx inewtgt = newtgtidx[0];

    Idx joldtgt = oldtgtidx[1];
    Idx jnewtgt = newtgtidx[1];

    // we don't use get/set_NeighborMat_bc here to avoid redundant bounds checking
    Vertex oldtarget = Acopy.data[ioldtgt * Acopy.cols + joldtgt];
    Vertex newtarget = Acopy.data[inewtgt * Acopy.cols + jnewtgt];
    A->data[perm.i * Acopy.cols + perm.j] = source ^ sourcemask;
    A->data[ioldtgt * Acopy.cols + joldtgt] = oldtarget ^ oldtgtmask;
    A->data[inewtgt * Acopy.cols + jnewtgt] = newtarget ^ newtgtmask;

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
    if (A->data != NULL) {
        free(A->data);
        A->data = NULL;
    }
}







