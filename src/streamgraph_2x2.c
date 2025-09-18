#include <stdint.h>
#include <stdlib.h>
#include <stdbool.h>

#include "returncodes.h"

typedef int32_t drainedarea_t;
typedef uint16_t cartidx_t;
typedef uint32_t linidx_t;
typedef uint8_t localedges_t; // 8 bits representing edges to 8 neighbors
typedef uint8_t clockhand_t; // 0-7. 0 = top, 1 = top right, 2 = right, etc. 255 means no neighbor.
#define IS_ROOT (clockhand_t)255

// [1111|2222|345-] 12B
typedef struct{
    drainedarea_t da;
    linidx_t adown; // linear index of the downstream vertex
    localedges_t edges;
    clockhand_t downstream;
    uint8_t visited;
} Vertex;

typedef struct {cartidx_t m; cartidx_t n; Vertex *vertices;} StreamGraph;

// ##############################
// # Helpers
// ##############################
static const cartidx_t row_offsets[8] = {-1, -1, 0, 1, 1, 1, 0, -1};
static const cartidx_t col_offsets[8] = {0, 1, 1, 1, 0, -1, -1, -1};

// ##############################
// # Create/destroy streamgraph #
// ##############################
Status new_StreamGraph(StreamGraph *G, cartidx_t m, cartidx_t n){
    if (m % 2 != 0 || n % 2 != 0) return MALFORMED_GRAPH_ERROR;
    if (G == NULL) return NULL_POINTER_ERROR;
    Vertex *vertices = calloc(m*n, sizeof(Vertex));
    if (vertices == NULL) return NULL_POINTER_ERROR;
    *G = (StreamGraph){m, n, vertices};
    return SUCCESS;
}
Status free_StreamGraph(StreamGraph *G){
    if (G != NULL && G->vertices != NULL){
        free(G->vertices);
        G->vertices = NULL;
    }
    G = NULL;
    return SUCCESS;
}


// ##############################
// # Getters + Setters          #
// ##############################

linidx_t cart_to_lin(row, col, m, n){
    div_t rdiv2 = div(row, 2);
    div_t cdiv2 = div(col, 2);
    linidx_t k = rdiv2.quot * n/2 * cdiv2.quot;
    linidx_t a = 4 * k + 2 * rdiv2.rem + cdiv2.rem;
    return a;
}
cartidx_t lin_to_cart_row(a, n){
    // c % 2 = a % 2
    // r % 2 = a % 2 + 1
    // k = ....?
    return 0;
}

// Gets the value at [row, col] on the map and loads into *out, with safeguards.
Status get_StreamGraph_cart_safe(StreamGraph G, Vertex *out, cartidx_t row, cartidx_t col){
    if (out == NULL || G.vertices == NULL) return NULL_POINTER_ERROR;
    if (row < 0 || row >= G.m || col < 0 || col >= G.n) return OOB_ERROR;
    div_t rdiv2 = div(row, 2);
    div_t cdiv2 = div(col, 2);
    linidx_t a = G.n/2 * rdiv2.quot + 2*rdiv2.rem + cdiv2.quot + cdiv2.rem;
    *out = G.vertices[a];
    return SUCCESS;
}
// Gets the value at [row, col] on the map and loads into *out, without safeguards.
Status get_StreamGraph_cart(StreamGraph G, Vertex *out, cartidx_t row, cartidx_t col){
    div_t rdiv2 = div(row, 2);
    div_t cdiv2 = div(col, 2);
    linidx_t a = G.n/2 * rdiv2.quot + 2*rdiv2.rem + cdiv2.quot + cdiv2.rem;
    *out = G.vertices[a];
    return SUCCESS;
}
// Sets the value at [row, col] on the map from *vert, with safeguards.
Status set_StreamGraph_cart_safe(Vertex vert, StreamGraph *G, cartidx_t row, cartidx_t col){
    if (G == NULL || G->vertices == NULL) return NULL_POINTER_ERROR;
    if (row < 0 || row >= G->m || col < 0 || col >= G->n) return OOB_ERROR;
    div_t rdiv2 = div(row, 2);
    div_t cdiv2 = div(col, 2);
    linidx_t a = G->n/2 * rdiv2.quot + 2*rdiv2.rem + cdiv2.quot + cdiv2.rem;
    G->vertices[a] = vert;
    return SUCCESS;
}
// Sets the value at [row, col] on the map from *vert, without safeguards.
Status set_StreamGraph_cart(Vertex vert, StreamGraph *G, cartidx_t row, cartidx_t col){
    div_t rdiv2 = div(row, 2);
    div_t cdiv2 = div(col, 2);
    linidx_t a = G->n/2 * rdiv2.quot + 2*rdiv2.rem + cdiv2.quot + cdiv2.rem;
    G->vertices[a] = vert;
    return SUCCESS;
}
// Gets the value at [a] in memory and loads into *out, with safeguards.
Status get_StreamGraph_lin_safe(StreamGraph G, Vertex *out, linidx_t a){
    if (a < 0 || a >= (G.m * G.n)) return OOB_ERROR;
    if (out == NULL || G.vertices == NULL) return NULL_POINTER_ERROR;
    *out = G.vertices[a];
    return SUCCESS;
}
// Gets the value at [a] in memory and loads into *out, with safeguards.
Status get_StreamGraph_lin(StreamGraph G, Vertex *out, linidx_t a){
    *out = G.vertices[a];
    return SUCCESS;
}
// Sets the value at [a] in memory to vert, with safeguards.
Status set_StreamGraph_lin_safe(Vertex vert, StreamGraph *G, linidx_t a){
    if (G == NULL || G->vertices == NULL) return NULL_POINTER_ERROR;
    if (a < 0 || a >= (G->m * G->n)) return OOB_ERROR;
    G->vertices[a] = vert;
    return SUCCESS;
}
// Gets the value at [a] in memory to vert, with safeguards.
Status set_StreamGraph_lin(Vertex vert, StreamGraph *G, linidx_t a){
    G->vertices[a] = vert;
    return SUCCESS;
}


// ##############################
// # Network traversal          #
// ##############################

// moves a downstream edge to point to a new node. Only checks for immediate issues in swapping
Status change_vertex_outflow(linidx_t a, StreamGraph *G, clockhand_t downnew){
    // 1. Get G[a], G[adownold], G[adownnew] safely
    // 2. Check that everything is valid regarding the swap
    // 3. Make the swap, adjust old and new downstream edges accordingly

    // 1.
    if (G == NULL) return NULL_POINTER_ERROR;
    StreamGraph Gcopy = *G;

    Status code;

    Vertex vert;
    Vertex vdownold;
    Vertex vdownnew;
    
    code = get_StreamGraph_lin_safe(Gcopy, &vert, a);
    if (code != SUCCESS) return code;

    linidx_t adownold = vert.adown;
    code = get_StreamGraph_lin_safe(Gcopy, &vdownold, adownold);
    if (code != SUCCESS) return code;

    linidx_t adownnew = cart_to_lin(row_offsets[downnew], col_offsets[downnew], Gcopy.m, Gcopy.n);
    code = get_StreamGraph_lin_safe(Gcopy, &vdownnew, adownnew);
    if (code != SUCCESS) return code;

    // 2. 
    // (Could move to just below get_StreamGraph_lin_safe(Gcopy, &vert, a))
    if (
        vert.downstream == IS_ROOT  // root node
        || downnew == vert.downstream  // no change
        || 1 << downnew & vert.edges != 0  // new downstream not a valid edge
    ) return SWAP_WARNING;

    // 3.
    vert.adown = adownnew;
    vert.downstream = downnew;
    vert.edges = vert.edges ^ (vert.edges || (1 << downnew));
    set_StreamGraph_lin(vert, G, a);

    vdownold.edges = vdownold.edges ^ (1 << (downnew + 4));
    set_StreamGraph_lin(vdownold, G, adownold);
    vdownnew.edges = vdownnew.edges ^ (1 << (downnew + 4));
    set_StreamGraph_lin(vdownnew, G, adownnew);

    return SUCCESS;
}

// starting from vertex linear (in-memory) index a, flow downstream and increment nodes by inc until reaching the root node
// visitationinc tracks how many times this function has been called since .vistitations was reset to 0. Can be either 1 or 2
static const uint8_t ncallsregistry[2] = {};

Status increment_downstream_safe(drainedarea_t inc, StreamGraph *G, linidx_t a, uint8_t ncalls){
    uint8_t visitationinc = ncallsregistry[ncalls];

    if (G == NULL) return NULL_POINTER_ERROR;
    StreamGraph Gcopy = *G;

    Vertex vert;
    Status code;

    code = get_StreamGraph_lin_safe(Gcopy, &vert, a);
    if (code != SUCCESS) return code;
    while (vert.downstream != IS_ROOT){
        // if we find ourselves in a cycle, exit immediately and signal to the caller
        if (vert.visited == ncalls) return CYCLE_WARNING;
        vert.visited += ncalls;
        vert.da += inc;
        
        // any other errorcode: exit
        code = get_StreamGraph_lin_safe(Gcopy, &vert, vert.adown);
        if (code != SUCCESS) return code;
    }

    return SUCCESS;
}
Status increment_downstream(drainedarea_t inc, StreamGraph *G, linidx_t a, uint8_t ncalls){
    uint8_t visitationinc = ncallsregistry[ncalls];
    StreamGraph Gcopy = *G;
    Vertex vert;
    Status code;

    get_StreamGraph_lin(Gcopy, &vert, a);
    while (vert.downstream != IS_ROOT){
        // if we find ourselves in a cycle, exit immediately and signal to the caller
        if (vert.visited == ncalls) return CYCLE_WARNING;
        vert.visited += ncalls;
        vert.da += inc;
        
        // any other errorcode: exit
        get_StreamGraph_lin(Gcopy, &vert, vert.adown);
    }

    return SUCCESS;
}

Status main_ocn_loop(uint8_t maxtries, StreamGraph *G){
    linidx_t a;
    clockhand_t downnew;
    Status code;
    uint8_t ntries = 0;

    if (G == NULL) return NULL_POINTER_ERROR;
    if (G->vertices == NULL) return NULL_POINTER_ERROR;
    StreamGraph Gcopy = *G;

     do {
        a = rand() % (G->m * G->n);
        downnew = rand();
        
        code = change_vertex_outflow(a, G, downnew);
        if (code == SWAP_WARNING) continue;
        if (code != SUCCESS) return code;  // not always an error: could be just an invalid swap
        
        
        linidx_t adownnew = cart_to_lin(row_offsets[downnew], col_offsets[downnew], Gcopy.m, Gcopy.n);
        Vertex vert, vdownold, vdownnew;
        get_StreamGraph_lin(Gcopy, &vert, a);
        get_StreamGraph_lin(Gcopy, &vdownold, vert.adown);
        get_StreamGraph_lin(Gcopy, &vdownnew, adownnew);
        
        drainedarea_t inc = vert.da;
        // zero the visitation counter before doing anything else
        for (linidx_t i = 0; i < Gcopy.m * Gcopy.n; i++) G->vertices[i].visited = 0;  
        code = increment_downstream_safe(-inc, G, vert.adown, 1);
        if (code == CYCLE_WARNING){  // if we find ourselves in a cycle, retrace our steps and undo the increment
            for (linidx_t i = 0; i < Gcopy.m * Gcopy.n; i++) G->vertices[i].visited = 0;
            increment_downstream(inc, G, vert.adown, 1);
            return SWAP_WARNING;
        }
        code = increment_downstream_safe(inc, G, adownnew, 2);
        if (code == CYCLE_WARNING){  // if we find ourselves in a cycle, retrace our steps and undo the increment
            for (linidx_t i = 0; i < Gcopy.m * Gcopy.n; i++) G->vertices[i].visited = 0;
            increment_downstream(inc, G, vert.adown, 1);
            increment_downstream(-inc, G, adownnew, 2);
            return SWAP_WARNING;
        }







        ntries += 1;
    }while ((code == CYCLE_WARNING || code == OOB_ERROR) && (ntries <= maxtries));
    
    return code;

    code = change_vertex_outflow(a, G, downnew);
    if (code != SUCCESS) return code;  // not always an error: could be just an invalid swap
    
    
    linidx_t adownnew = cart_to_lin(row_offsets[downnew], col_offsets[downnew], Gcopy.m, Gcopy.n);
    Vertex vert, vdownold, vdownnew;
    get_StreamGraph_lin(Gcopy, &vert, a);
    get_StreamGraph_lin(Gcopy, &vdownold, vert.adown);
    get_StreamGraph_lin(Gcopy, &vdownnew, adownnew);
    
    drainedarea_t inc = vert.da;
    // zero the visitation counter before doing anything else
    for (linidx_t i = 0; i < Gcopy.m * Gcopy.n; i++) G->vertices[i].visited = 0;  
    code = increment_downstream_safe(-inc, G, vert.adown, 1);
    if (code == CYCLE_WARNING){  // if we find ourselves in a cycle, retrace our steps and undo the increment
        for (linidx_t i = 0; i < Gcopy.m * Gcopy.n; i++) G->vertices[i].visited = 0;
        increment_downstream(inc, G, vert.adown, 1);
        return SWAP_WARNING;
    }
    code = increment_downstream_safe(inc, G, adownnew, 2);
    if (code == CYCLE_WARNING){  // if we find ourselves in a cycle, retrace our steps and undo the increment
        for (linidx_t i = 0; i < Gcopy.m * Gcopy.n; i++) G->vertices[i].visited = 0;
        increment_downstream(inc, G, vert.adown, 1);
        increment_downstream(-inc, G, adownnew, 2);
        return SWAP_WARNING;
    }

    return SUCCESS;
}

Status ocniterate(){
    linidx_t astart = rand() % (G->m * G->n);
    clockhand_t downnew = rand();
    Status code;
    uint8_t ntries = 0;
   
}