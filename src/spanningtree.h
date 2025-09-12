#ifndef SPANNINGTREE_H
#define SPANNINGTREE_H


#include <stdint.h>
#include <stdbool.h>
#include "array.h"

/**
 * @file spanningtree.h
 * @brief Definition and basic operations for spanning trees.
 */

/**
 * @struct SpatialEdge
 * @brief Represents an edge in a spatially embedded graph with source and destination coordinates.
 * @param src Coordinates of the source node, 2-element array (row, column)
 * @param dest Coordinates of the destination, 2-element array (row, column)
 */
typedef struct{
    int32_t src[2];
    int32_t dest[2];
} SpatialEdge;

/**
 * @struct SpanningTree
 * @brief A simple spanning tree structure.
 * @param n_nodes Number of nodes in the tree
 * @param n_edges Number of edges in the tree (should be n_nodes - 1)
 * @param edges Pointer to the array of edges (each edge is a pair of node indices)
 */
typedef struct {
    int32_t n_edges;
    SpatialEdge *edges;
} SpanningTree;

// Function prototypes for spanning tree operations
SpanningTree create_spanning_tree(int32_t n_edges);
void free_spanning_tree(SpanningTree *tree);
void swap_edges(SpanningTree *tree, int32_t idx1, int32_t idx2);
bool has_cycle(SpanningTree tree);
bool edges_cross(SpatialEdge e1, SpatialEdge e2);
SpatialEdge *parents(SpatialEdge edge, SpanningTree tree);
SpatialEdge *children(SpatialEdge edge, SpanningTree tree);

#endif // SPANNINGTREE_H