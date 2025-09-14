/*
TODO implement the NO_CLOCKHAND value.
TODO instill some consistency in our indexing into S. Currently we use a mixture of size-2 arrays and single indices.
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
#ifndef STREAMGRAPH_C
#define STREAMGRAPH_C

#include <stdint.h>
#include <stdlib.h>
#include <stdio.h>
#include <stdbool.h>
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
typedef uint8_t Clockhand; // 0-7. 0 = top, 1 = top right, 2 = right, etc. 255 means no neighbor.
#define NO_CLOCKHAND (Clockhand)255

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
 * O O O
 *  \|/
 * O-O-O
 *    /|
 * O-O R
 * 
 * Vertex 1,1 (in the middle) is connected to neighbors 0 (top), 1 (top right),
 * 2 (right), 6 (left), and 7 (top left), so S[1,1] = 0b11000111 = 0xc7
 * 
 * Additionally, StreamGraph contains a "downstream" array. For each element, it contains the downstream direction.
 */
typedef struct {
    Idx m;
    Idx n;
    Vertex *vertices;
} StreamGraph;

/**
 * @brief Retrieves the vertex stored at (i,j) in the StreamGraph S and stores it in vert. Performs bounds checking.
 * @param S The neighbor matrix
 * @param vert Pointer to store the retrieved vertex
 * @param i Row index
 * @param j Column index
 * @return SUCCESS or an error code
 */
Status get_StreamGraph_bc(StreamGraph S, Vertex *vert, Idx i, Idx j) {
    if (i < 0 || i >= S.m || j < 0 || j >= S.n) return OOB_ERROR;
    if (vert == NULL) return NULL_POINTER_ERROR;
    *vert = S.vertices[i * S.n + j];
    return SUCCESS;
}

/**
 * @brief Sets the vertex stored at (i,j) in the StreamGraph S to the given value. Performs bounds checking.
 * @param S Pointer to the neighbor matrix
 * @param i Row index
 * @param j Column index
 * @param vert The new vertex value to set
 * @return SUCCESS or an error code
 */
Status set_StreamGraph_bc(StreamGraph *S, Idx i, Idx j, Vertex vert) {
    if (i < 0 || i >= S->m || j < 0 || j >= S->n) return OOB_ERROR;
    if (S->vertices == NULL) return NULL_POINTER_ERROR;
    S->vertices[i * S->n + j] = vert;
    return SUCCESS;
}

static const Idx row_neighbors[8] = {-1, -1, 0, 1, 1, 1, 0, -1};
static const Idx col_neighbors[8] = {0, 1, 1, 1, 0, -1, -1, -1};

/**
 * @brief Get the row and column indices of neighboring vertices for a given cell.
 * @param S The adjacency matrix
 * @param out Array of size 2 to store the row and column indices of the neighbor
 * @param i Row index of the vertex in question
 * @param j Column index of the vertex in question
 * @return SUCCESS or an error code
 */
Status get_clockhand_neighbor(StreamGraph S, Idx out[2], Idx i, Idx j, Clockhand hand) {
    // returns the index of the neighbor pointed to by the clockhand (0-7) of the vertex at (i,j).
    Vertex vert;
    Status code = get_StreamGraph_bc(S, &vert, i, j);
    if (code != SUCCESS) return code;
    if (edge_exists(vert, hand) == 0) return NO_EDGE_WARNING; // no neighbor in that direction

    out[0] = i + row_neighbors[hand];
    out[1] = j + col_neighbors[hand];
    return SUCCESS;
}

/**
 * @brief Check if there is a cycle in the graph containing the start vertex.
 * @param S The neighbor matrix
 * @param start The starting vertex
 * @return true if there is a cycle, false otherwise
 */
Status in_cycle(StreamGraph S, bool *out, Idx i, Idx j){
    bool *visited = calloc(S.m * S.n, sizeof(bool));  // 0 means unvisited
    if (visited == NULL) return NULL_POINTER_ERROR;

    Vertex vert;
    Status code;
    code = get_StreamGraph_bc(S, &vert, i, j);
    if (code != SUCCESS) goto cleanup;
    visited[i*S.n + j] = true;

    Idx next_idx[2];
    while (vert.downstream != NO_CLOCKHAND){
        code = get_clockhand_neighbor(S, next_idx, i, j, vert.downstream);
        i = next_idx[0]; j = next_idx[1];
        if (code != SUCCESS) goto cleanup;
        
        if (visited[i*S.n + j] == true){
            *out = true;
            goto cleanup;
        }
        
        code = get_StreamGraph_bc(S, &vert, i, j);
        if (code != SUCCESS) goto cleanup;
        visited[i*S.n + j] = true;
    }
    *out = false;

    cleanup:
        free(visited);
        return code;
}

/**
 * @brief Move the outgoing edge from direction handa to handb for the vertex at (i,j) in the neighbor matrix S.
 * @param S Pointer to the neighbor matrix
 * @param i Row index of the vertex
 * @param j Column index of the vertex
 * @param newdownstream The new downstream direction
 * @return SUCCESS or an error code. 
 */
Status swap_edges(StreamGraph *S, Idx i, Idx j, Clockhand newdownstream) {
    StreamGraph Acopy = *S ; // avoid dereferencing S multiple times. Used only for lookup, not modification.
    Vertex sourcenode, downnode, newdownnode;
    Status code;
    
    // confirm that the source node exists
    code = get_StreamGraph_bc(Acopy, &sourcenode, i, j);
    if (code != SUCCESS) return code;
    if (newdownstream == sourcenode.downstream) return SWAP_WARNING; // nothing to do
    if (edge_exists(sourcenode, newdownstream)) return SWAP_WARNING; // can't swap to an upstream edge
    // So long as the initial array was constructed correctly, this should not be necessary.
    // could remove later.
    if (!downstream_edge_exists(sourcenode, sourcenode.downstream)) return NO_EDGE_WARNING;  

    LocalEdges sourcemask = (1 << sourcenode.downstream) | (1 << newdownstream);
    LocalEdges downmask = (1 << sourcenode.downstream) + 4;
    LocalEdges newdownmask = (1 << newdownstream) + 4;

    // we don't use get/set_StreamGraph_bc here to avoid redundant bounds checking
    Idx downidx[2], newdownidx[2];
    code = get_clockhand_neighbor(Acopy, downidx, i, j, sourcenode.downstream);
    if (code != SUCCESS) return code;
    downnode = Acopy.vertices[downidx[0] * Acopy.n + downidx[1]];
    code = get_clockhand_neighbor(Acopy, newdownidx, i, j, newdownstream);
    if (code != SUCCESS) return code;
    newdownnode = Acopy.vertices[newdownidx[0] * Acopy.n + newdownidx[1]];

    S->vertices[i * Acopy.n + j] = (Vertex){sourcenode.edges ^ sourcemask, newdownstream};

    // check to see if this swap created a cycle
    bool created_cycle;
    code = in_cycle(Acopy, &created_cycle, i, j);
    if (code != SUCCESS){
        S->vertices[i * Acopy.n + j] = (Vertex){sourcenode.edges, sourcenode.downstream};
        return code;
    }
    else if (created_cycle){
        // revert the change
        S->vertices[i * Acopy.n + j] = (Vertex){sourcenode.edges, sourcenode.downstream};
        return SWAP_WARNING;
    }

    S->vertices[downidx[0] * Acopy.n + downidx[1]] = (Vertex){
        downnode.edges ^ downmask, 
        downnode.downstream
    };
    S->vertices[newdownidx[0] * Acopy.n + newdownidx[1]] = (Vertex){
        newdownnode.edges ^ newdownmask, 
        newdownnode.downstream
    };

    return SUCCESS;
}



/**
 * @brief Print a visual representation of the vertex to stdout.
 * @param vert The vertex to print
 */
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

void free_StreamGraph(StreamGraph *S) {
    if (S->vertices != NULL) {
        free(S->vertices);
        S->vertices = NULL;
    }
}

#endif // STREAMGRAPH_C