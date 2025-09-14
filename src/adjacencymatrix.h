#ifndef ADJACENCYMATRIX_H
#define ADJACENCYMATRIX_H

#include <stdint.h>
#include <stdbool.h>

/**
 * @file adjacencymatrix.h
 * @brief Definition of adjacency matrix and related operations.
 */

typedef int16_t Idx; // for indexing rows and columns
typedef uint16_t Vertex; // 16-bit unsigned int representing edges to 8 neighbors. Values over 255 are used as return codes.
typedef uint8_t Clockhand; // 0-7

/**
 * @brief Neighbor matrix structure.
 * 
 * The neighbor matrix is a 2D array where each element indicates the presence
 * (1) or absence (0) of an edge between nodes in a graph. Each element N[i,j] 
 * represents the edges connecting the node at row i, column j to its
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
 * Node 1,1 (in the middle) is connected to neighbors 0 (top), 1 (top right),
 * 2 (right), 6 (left), and 7 (top left), so A[1,1] = 0b11000111 = 0xc7
 */
typedef struct {
    Idx rows;
    Idx cols;
    Vertex *data;
} NeighborMat;

/**
 * @brief Contains the row and column indices of neighboring vertices.
 * @param n the number of neighbors
 * @param row array of row indices of neighbors
 * @param col array of column indices of neighbors
 */
typedef struct{
    Vertex node;
    Idx row[8];
    Idx col[8];
} NeighborVertices;

Vertex get_NeighborMat_bc(NeighborMat A, Idx i, Idx j);
void set_NeighborMat_bc(NeighborMat *A, Idx i, Idx j, Vertex node);

/**
 * @brief Get the row and column indices of neighboring vertices for a given cell.
 * @param A The adjacency matrix
 * @param i Row index of the cell
 * @param j Column index of the cell
 * @return A NeighborVertices struct containing the number of neighbors and their indices
 */
NeighborVertices get_neighbor_vertices_bc(NeighborMat A, Idx i, Idx j);

/**
 * @brief Get the indices of the neighbor in the direction indicated by the clockhand.
 * @param A The adjacency matrix
 * @param i Row index of the cell
 * @param j Column index of the cell
 * @param clockhand An integer from 0 to 7 indicating the direction (0=top, 1=top-right, ..., 7=top-left)
 * @return An array of two integers [i, j] representing the indices of the neighbor, or NULL if no neighbor exists in that direction
 */
int32_t *get_clockhand_neighbor_bc(NeighborMat A, Idx i, Idx j, Clockhand hand);

/**
 * @brief Represents a permutation of a vertex in the adjacency matrix.
 * Used for reordering vertices.
 * @param i Row index of the vertex
 * @param j Column index of the vertex
 * @param oldclockhand The original neighbor index of the vertex as the bit number (0=top, 1=top-right, 2=right,etc.)
 * @param newclockhand The new neighbor value (bitmask) of the vertex after permutation
 */
typedef struct {
    Idx i;
    Idx j;
    Clockhand oldhand;
    Clockhand newhand;
} Permutation;

/**
 * @brief Permute a vertex in the adjacency matrix according to the given permutation mapping.
 * @param A Pointer to the adjacency matrix to be modified
 * @param perm The permutation mapping
 * @return 0 if the permutation was successful, 1 if the old edge does not exist, 2 if the new edge already exists, or 3 if both errors occurred
 */
int permute_vertices_bc(NeighborMat *A, Permutation perm);

#endif // ADJACENCYMATRIX_H