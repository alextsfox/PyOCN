#include <stdint.h>
#include <stdlib.h>
#include <stdbool.h>

#include "returncodes.h"
#include "neighbormatrix.c"
#include "idx.c"

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




/*
TODO create functions to
    (1) return all downstream vertices from a given vertex, in order
    (2) increment/decrement the drained area of all vertices downstream of a given vertex.

    These could all be combined into a single function.
    
    # this function could be run when we instantiate the streamgraph on all nodes to check for cycles
    # it can be run again after a swap to check if the swap created a cycle
    Python pseudocode:
    def flow_downstream(S:StreamGraph, i:Idx, j:Idx):
        visited = [0]*S.m*S.n
        downstream_nodes_i = [i]
        downstream_nodes_j = [j]

        vert = S[i][j]
        visited[i*S.n + j] += 1
        while vert.downstream != NO_CLOCKHAND:  # while vert is not root
            i, j = get_clockhand_neighbor(S, i, j)
            if visited[i*S.n + j] >= 1:
                raise CYCLE_ERROR
            visited[i*S.n + j] += 1
            downstream_nodes_i.append(i)
            downstream_nodes_j.append(j)
            vert = S[i][j]
        return visited, downstream_nodes_i, downstream_nodes_j
    
    # this code flows downstream like the one above, but does not need to keep track of visited nodes, and possibly doesn't even need to track cycles.
    def modify_downstream_drained_area(S:StreamGraph, drainedarea:DrainedAreaGraph, i:Idx, j:Idx, delta:int):
        visited = [0]*S.m*S.n  # can skip if we are confidant that the graph is not malformed/does not contain a cycle. Could be a big time saver.
        vert = S[i][j]
        visited[i*S.n + j] += 1  # again, could skip
        while vert.downstream != NO_CLOCKHAND:  # while vert is not root
            i, j = get_clockhand_neighbor(S, i, j)
            if visited[i*S.n + j] >= 1:  # again, could skip
                raise CYCLE_ERROR  # to handle this, we will have to revert all changes with another call to this function, but with a negative delta.
            drainedarea[i*S.n + j] += delta
            visited[i*S.n + j] += 1  # again, could skip
            vert = S[i][j]
        return drainedarea

    # right now the swap function checks for cycles by itself. 
    # we could instead, in our main function, call swap, then check for cycles manually when we modify the streamgraph, which might save some time.
    # however, this would require use to re-run modify_downstream_drained_area with a negative delta, losing us time.

    # so we have two options:
    # (1) keep the swap function as is, and have it check for cycles itself
    # (2) have the swap function not check for cycles, and instead check for cycles in the main function after every swap, and simultaneously compute the drained area, retracing our steps if need be.
    # (3) have the swap function check for cycles and return the list of affected nodes (in reality this list would be pre-allocated and only modified by the swap function), and then modify the drained area in the main function using that list if no cycles were detected.
    
    # I think (1) is best for now, since it is simplest to implement and understand. We can optimize later if needed.
    # (2) Saves on memory at the cost of additional CPU cycles, since we need to traverse the tree three times: once to check for cycles, once to decrement drained area downstream of the old downnode, and once to increment the drained area downstream of the new downnode.
    # (3) Saves on CPU cycles at the cost of additional memory, since we have to store the list of affected nodes, but would only have to traverse part of the tree each time to modify the drained area.

    # I think the additional memory cost (on the order of ~1MB) is worth cutting our CPU cycles by 75%, so I will implement (3).
*/

typedef int32_t DrainedArea;

/**
 * @brief Rotate the downstream edge on vertex idx to newdownstream for the vertex at (i,j) in the stream graph S, and update the drained area accordingly.
 * Insures that any input array remains valid after the function returns.
 * @param S Pointer to the stream graph
 * @param drainedarea Pre-allocated drained area array to be updated
 * @param olddownnodes Pre-allocated array to store indices of vertices downstream of olddownnode that may be affected by the swap
 * @param newdownnodes Pre-allocated array to store indices of vertices downstream of newdownnode that may be affected by the swap
 * @param visited Pre-allocated FALSE boolean array to track visited nodes during cycle detection. MUST BE INITIALIZED TO FALSE BEFORE CALLING THIS FUNCTION.
 * @param idx The index of the vertex to modify
 * @return SUCCESS or an error code. 
 */
Status swap_edges_and_compute_drained_area(
    StreamGraph *S, 
    DrainedArea *drainedarea,
    IdxArray *olddownnodes,
    IdxArray *newdownnodes,
    bool *visited,
    IdxPair idx, 
    Clockhand newdownstream, 
){
    Status code;
    code = swap_edges(S, idx, newdownstream);
    if (code == SUCCESS || code == NULL_POINTER_ERROR || code == OOB_ERROR) return code;

    // check for cycles and return the list of affected nodes 
    // (in reality this list would be pre-allocated and only modified by the swap function), 
    // and then modify the drained area in the main function using that list if no cycles were detected.

    // visited = [0]*S.m*S.n  # can skip if we are confidant that the graph is not malformed/does not contain a cycle. Could be a big time saver.
    //     vert = S[i][j]
    //     visited[i*S.n + j] += 1  # again, could skip
    //     while vert.downstream != NO_CLOCKHAND:  # while vert is not root
    //         i, j = get_clockhand_neighbor(S, i, j)
    //         if visited[i*S.n + j] >= 1:  # again, could skip
    //             raise CYCLE_ERROR  # to handle this, we will have to revert all changes with another call to this function, but with a negative delta.
    //         drainedarea[i*S.n + j] += delta
    //         visited[i*S.n + j] += 1  # again, could skip
    //         vert = S[i][j]
    //     return drainedarea
}