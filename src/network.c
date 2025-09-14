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
algorithm:
1. Starting from root, follow upstream edges up the tree. Each time we visit a vertex or a predecessor vertex of that vertex, increment its drained area by 1.
2. Stop when we reach a vertex with no upstream edges.
N.B. we only need to compute this once. After a swap, we can do the following
Case 1 (circular river): Swap created a cycle
Case 2 (shortcut): Swap moved outflow downstream from where if was before
Case 3 (longcut?): Swap moved outflow upstream from where it was before
Case 4 (stream capture): Swap moved outflow to a separate branch from where it was before.

Case 1:
    Do nothing, this should have been handled by an error code in the swap_nodes function
Case 2:
    Starting from the previous downnode, subtract sourcenode's upstream area from that node and all downstream nodes. 
    Stop when we hit the new downnode (do not modify drained area of new downnode)
Case 3:
    Starting from the new downnode, add sourcenode's upstream area to that node and all downstream nodes.
    Stop when we hit the old downnode (do not modify drained area of old downnode)
Case 4 (hard):
    Find the least common acnestor (or descendant, I suppose, since the DAG points downstream) of the old and new downnodes.
    Starting from the old downnode, subtract the sourcenode's from that node and all downstream nodes.
    Stop when we hit the common ancestor (do not modify drained area of common ancestor)
    Starting from the new downnode, add the sourcenode's upstream area to that node and all downstream nodes.
    Stop when we hit the common ancestor (do not modify drained area of common ancestor)
*/
// In python, this would be something like:
// def set_node_predecessors(DAG, node):
//     predecessors = {node}
//     for n in DAG.predecessors(node):
//         predecessors = predecessors | DAG.nodes[n]['predecessors']
//     DAG.nodes[node]['predecessors'] = predecessors
int32_t compute_drained_area(NeighborMat A, Vertex root, int32_t *drained_area) {
    return SUCCESS;
}