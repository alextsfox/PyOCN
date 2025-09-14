#ifndef STREAMGRAPH_C
#define STREAMGRAPH_C

#include <stdint.h>
#include <stdlib.h>
#include <stdio.h>
#include <stdbool.h>
#include "returncodes.h"
#include "idx.c"


typedef uint8_t LocalEdges; // 8 bits representing edges to 8 neighbors
typedef uint8_t Clockhand; // 0-7. 0 = top, 1 = top right, 2 = right, etc. 255 means no neighbor.
#define NO_CLOCKHAND (Clockhand)255
typedef struct{LocalEdges edges; Clockhand downstream;} Vertex;  // downstream indicates the direction of the downstream vertex
typedef struct {Idx m; Idx n; Vertex *vertices;} StreamGraph;

#define edge_exists(vert, hand) (vert.edges & (1 << hand))
#define downstream_edge_exists(vert, hand) (edge_exists(vert, hand) & (hand == vert.downstream))
#define upstream_edge_exists(vert, hand) (edge_exists(vert, hand) & (hand != vert.downstream))

Status get_StreamGraph_bc(StreamGraph S, Vertex *vert, IdxPair idx) {
    if (idx[0] < 0 || idx[0] >= S.m || idx[1] < 0 || idx[1] >= S.n) return OOB_ERROR;
    if (vert == NULL) return NULL_POINTER_ERROR;
    *vert = S.vertices[idx[0] * S.n + idx[1]];
    return SUCCESS;
}
Status set_StreamGraph_bc(StreamGraph *S, IdxPair idx, Vertex vert) {
    if (idx[0] < 0 || idx[0] >= S->m || idx[1] < 0 || idx[1] >= S->n) return OOB_ERROR;
    if (S->vertices == NULL) return NULL_POINTER_ERROR;
    S->vertices[idx[0] * S->n + idx[1]] = vert;
    return SUCCESS;
}

static const IdxPair neighbor_offsets[8] = {
    (IdxPair){-1, -1}, (IdxPair){-1, 0}, (IdxPair){-1, 1}, 
    (IdxPair){0, -1}, (IdxPair){0, 1},   
    (IdxPair){1, -1}, (IdxPair){1, 0}, (IdxPair){1, 1}
};

/**
 * @brief Get the row and column indices of neighboring vertices for a given cell.
 * @param S The adjacency matrix
 * @param out Array of size 2 to store the row and column indices of the neighbor
 * @param i Row index of the vertex in question
 * @param j Column index of the vertex in question
 * @return SUCCESS or an error code
 */
Status get_clockhand_neighbor(StreamGraph S, IdxPair *out, IdxPair idx, Clockhand hand) {
    // returns the index of the neighbor pointed to by the clockhand (0-7) of the vertex at (i,j).
    Vertex vert;
    Status code = get_StreamGraph_bc(S, &vert, idx);
    if (code != SUCCESS) return code;
    if (vert.downstream == NO_CLOCKHAND) return NO_EDGE_WARNING; // root node
    if (edge_exists(vert, hand) == 0) return NO_EDGE_WARNING; // no neighbor in that direction

    out->row = idx.row + neighbor_offsets[hand].row;
    out->col = idx.col + neighbor_offsets[hand].col;
    return SUCCESS;
}

/**
 * @brief Move the outgoing edge from direction handa to handb for the vertex at (i,j) in the neighbor matrix S.
 * Will swap edges only if they are are not immediately obviously invalid (ie will not check for large cycles).
 * @param S Pointer to the neighbor matrix
 * @param i Row index of the vertex
 * @param j Column index of the vertex
 * @param newdownstream The new downstream direction
 * @return SUCCESS or an error code. 
 */
Status swap_edges(StreamGraph *S, IdxPair idx, Clockhand newdownstream) {
    StreamGraph Acopy = *S ; // avoid dereferencing S multiple times. Used only for lookup, not modification.
    Vertex sourcenode, downnode, newdownnode;
    Status code;
    
    // confirm that the source node exists
    code = get_StreamGraph_bc(Acopy, &sourcenode, idx);
    if (code != SUCCESS) return code;
    
    // once again, we could remove this check later once we work out a good method to prevent it from happening.
    if (sourcenode.downstream == NO_CLOCKHAND) return SWAP_WARNING; // can't swap if there is no downstream edge
    if (newdownstream == sourcenode.downstream) return SWAP_WARNING; // nothing to do
    if (edge_exists(sourcenode, newdownstream)) return SWAP_WARNING; // can't swap to an upstream edge
    // So long as the initial array was constructed correctly, this should not be necessary.
    // could remove later.
    if (!downstream_edge_exists(sourcenode, sourcenode.downstream)) return NO_EDGE_WARNING;  

    LocalEdges sourcemask = (1 << sourcenode.downstream) | (1 << newdownstream);
    LocalEdges downmask = (1 << sourcenode.downstream) + 4;
    LocalEdges newdownmask = (1 << newdownstream) + 4;

    // we don't use get/set_StreamGraph_bc here to avoid redundant bounds checking
    IdxPair downidx, newdownidx;
    code = get_clockhand_neighbor(Acopy, &downidx, idx, sourcenode.downstream);
    if (code != SUCCESS) return code;
    downnode = Acopy.vertices[downidx.row * Acopy.n + downidx.col];
    code = get_clockhand_neighbor(Acopy, &newdownidx, idx, newdownstream);
    if (code != SUCCESS) return code;
    newdownnode = Acopy.vertices[newdownidx.row * Acopy.n + newdownidx.col];

    S->vertices[idx.row * Acopy.n + idx.col] = (Vertex){sourcenode.edges ^ sourcemask, newdownstream};

    S->vertices[downidx.row * Acopy.n + downidx.col] = (Vertex){
        downnode.edges ^ downmask, 
        downnode.downstream
    };
    S->vertices[newdownidx.row * Acopy.n + newdownidx.col] = (Vertex){
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