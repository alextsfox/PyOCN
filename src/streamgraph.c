// naming convention
// defines: ALL_CAPS
// struct types: PascalCase
// simple types: flatcase_t
// functions + variables: snake_case

#ifndef STREAMGRAPH_C
#define STREAMGRAPH_C


#include <stdint.h>
#include <stdlib.h>
#include <stdbool.h>
#include <math.h>
#include <stdio.h>

#include "returncodes.h"

const double gamma = 0.5;
const double temperature = 0.75;

/**
 * @file streamgraph.h
 * @brief Data structures and functions for managing stream graphs.
 * A stream graph represents a tree of vertices with directed edges,
 * where each vertex has a drained area and flows downstream to another vertex.
 * The graph is stored in a 2D grid layout for spatial representation, and laid out in memory using block tiling.
 * The edges connected to each vertex are represented using an 8-bit integer, where each bit corresponds to one of the 8 possible neighboring directions.
 * The downstream direction is represented as a clock hand (0-7), with 255 indicating no downstream (root).
 */

#include <stdint.h>
#include <stdbool.h>
#include "returncodes.h"

// Type definitions
typedef int32_t drainedarea_t;
typedef uint16_t cartidx_t;
typedef uint32_t linidx_t;
typedef uint8_t localedges_t;
typedef uint8_t clockhand_t;

#define IS_ROOT (clockhand_t)255
#define TILE_SIZE 2

typedef struct {
    drainedarea_t drained_area;
    linidx_t adown;
    localedges_t edges;
    clockhand_t downstream;
    uint8_t visited;
} Vertex;

typedef struct {
    cartidx_t m, n;
    cartidx_t i_root, j_root;
    double energy;
    Vertex *vertices;
} StreamGraph;

// ##############################
// # Helpers
// ##############################
const int16_t row_offsets[8] = {-1, -1, 0, 1, 1, 1, 0, -1};
const int16_t col_offsets[8] = {0, 1, 1, 1, 0, -1, -1, -1};

// ##############################
// # Create/destroy streamgraph #
// ##############################
Status sg_create(StreamGraph *G, cartidx_t m, cartidx_t n, cartidx_t i_root, cartidx_t j_root){
    if (m % 2 != 0 || n % 2 != 0) return MALFORMED_GRAPH_ERROR;
    if (G == NULL) return NULL_POINTER_ERROR;
    Vertex *vertices = calloc(m*n, sizeof(Vertex));
    if (vertices == NULL) return NULL_POINTER_ERROR;
    *G = (StreamGraph){m, n, i_root, j_root, 0.0, vertices};
    return SUCCESS;
}
Status sg_destroy(StreamGraph *G){
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

// 2d raster is block-tiled in memory to improve cache locality
#define TILE_SIZE 2
/*
 * Converts cartesian (row, col) to linear index in block-tiled memory
 */
// linidx_t sg_cart_to_lin(cartidx_t row, cartidx_t col, cartidx_t m, cartidx_t n){
//     div_t r_divT = div(row, TILE_SIZE);
//     div_t c_divT = div(col, TILE_SIZE);
//     linidx_t k = r_divT.quot * n/TILE_SIZE + c_divT.quot;
//     linidx_t a = (TILE_SIZE * TILE_SIZE * k) + (TILE_SIZE * r_divT.rem) + c_divT.rem;
//     return a;
// }
// for now, we just use normal row-major order until we finish debugging the main algorithm
linidx_t sg_cart_to_lin(cartidx_t row, cartidx_t col, cartidx_t m, cartidx_t n){
    return (row * n + col);
    // return a;
}




// absolute nightmare of a function, but it works for now
// TODO: write a proper linear to cartesian conversion function
cartidx_t sg_lin_to_cart_row(linidx_t a, cartidx_t m, cartidx_t n){
    for (cartidx_t row = 0; row < m; row++){
        for (cartidx_t col = 0; col < n; col++){
            if (a == sg_cart_to_lin(row, col, m, n)) return row;
        }
    }
    return 0;
}

cartidx_t sg_lin_to_cart_col(linidx_t a, cartidx_t m, cartidx_t n){
    for (cartidx_t row = 0; row < m; row++){
        for (cartidx_t col = 0; col < n; col++){
            if (a == sg_cart_to_lin(row, col, m, n)) return col;
        }
    }
    return 0;
}

// Gets the value at [row, col] on the map and loads into *out, with safeguards.
Status sg_get_cart_safe(StreamGraph G, Vertex *out, cartidx_t row, cartidx_t col){
    if (row < 0 || row >= G.m || col < 0 || col >= G.n) return OOB_ERROR;
    linidx_t a = sg_cart_to_lin(row, col, G.m, G.n);
    *out = G.vertices[a];
    return SUCCESS;
}

// Gets the value at [row, col] on the map and loads into *out, without safeguards.
Status sg_get_cart(StreamGraph G, Vertex *out, cartidx_t row, cartidx_t col){
    linidx_t a = sg_cart_to_lin(row, col, G.m, G.n);
    *out = G.vertices[a];
    return SUCCESS;
}

// Sets the value at [row, col] on the map from *vert, with safeguards.
Status sg_set_cart_safe(Vertex vert, StreamGraph *G, cartidx_t row, cartidx_t col){
    if (row < 0 || row >= G->m || col < 0 || col >= G->n) return OOB_ERROR;
    linidx_t a = sg_cart_to_lin(row, col, G->m, G->n);
    G->vertices[a] = vert;
    return SUCCESS;
}
// Sets the value at [row, col] on the map from *vert, without safeguards.
Status sg_set_cart(Vertex vert, StreamGraph *G, cartidx_t row, cartidx_t col){
    linidx_t a = sg_cart_to_lin(row, col, G->m, G->n);
    G->vertices[a] = vert;
    return SUCCESS;
}
// Gets the value at [a] in memory and loads into *out, with safeguards.
Status sg_get_lin_safe(StreamGraph G, Vertex *out, linidx_t a){
    if (a < 0 || a >= (G.m * G.n)) return OOB_ERROR;
    *out = G.vertices[a];
    return SUCCESS;
}
// Gets the value at [a] in memory and loads into *out, with safeguards.
Status sg_get_lin(StreamGraph G, Vertex *out, linidx_t a){
    *out = G.vertices[a];
    return SUCCESS;
}
// Sets the value at [a] in memory to vert, with safeguards.
Status sg_set_lin_safe(Vertex vert, StreamGraph *G, linidx_t a){
    if (a < 0 || a >= (G->m * G->n)) return OOB_ERROR;
    G->vertices[a] = vert;
    return SUCCESS;
}
// Gets the value at [a] in memory to vert, with safeguards.
Status sg_set_lin(Vertex vert, StreamGraph *G, linidx_t a){
    G->vertices[a] = vert;
    return SUCCESS;
}


// ##############################
// # Network traversal          #
// ##############################

/**
 * @brief Change the outflow direction of a vertex, updating the edges of the new and old downstream vertices accordingly.
 * @param a Linear index into the StreamGraph of the vertex to change
 * @param G Pointer to the StreamGraph
 * @param down_new The new downstream direction (0-7).
 * @return SUCCESS, SWAP_WARNING, or OOB_ERROR.
 */
Status sg_change_vertex_outflow(linidx_t a, StreamGraph *G, clockhand_t down_new){
    // 1. Get G[a], G[a_down_old], G[adownnew] safely
    // 2. Check that everything is valid regarding the swap
    // 3. Make the swap, adjust old and new downstream edges accordingly

    // 1.
    StreamGraph Gcopy = *G;

    Status code;

    Vertex vert;
    Vertex vert_down_old;
    Vertex vert_down_new;
    
    code = sg_get_lin_safe(Gcopy, &vert, a);
    if (code == OOB_ERROR) return OOB_ERROR;

    linidx_t a_down_old = vert.adown;
    code = sg_get_lin_safe(Gcopy, &vert_down_old, a_down_old);
    if (code == OOB_ERROR) return OOB_ERROR;

    linidx_t a_down_new = sg_cart_to_lin((cartidx_t)(sg_lin_to_cart_row(a, Gcopy.m, Gcopy.n) + row_offsets[down_new]), (cartidx_t)(sg_lin_to_cart_col(a, Gcopy.m, Gcopy.n) + col_offsets[down_new]), Gcopy.m, Gcopy.n);
    code = sg_get_lin_safe(Gcopy, &vert_down_new, a_down_new);
    if (code == OOB_ERROR) return OOB_ERROR;

    // 2. 
    // (Could move to just below sg_get_lin_safe(Gcopy, &vert, a))
    if (
        (vert.downstream == IS_ROOT)  // root node
        || (down_new == vert.downstream)  // no change
        || ((1 << down_new) & (vert.edges != 0)) // new downstream not a valid edge
    ) return SWAP_WARNING;

    // 3.
    vert.adown = a_down_new;
    clockhand_t down_old = vert.downstream;
    vert.downstream = down_new;
    vert.edges = vert.edges ^ (vert.edges | (1 << down_new));
    sg_set_lin(vert, G, a);

    vert_down_old.edges = vert_down_old.edges ^ (1 << (down_old + 4));
    sg_set_lin(vert_down_old, G, a_down_old);
    vert_down_new.edges = vert_down_new.edges ^ (1 << (down_new + 4));
    sg_set_lin(vert_down_new, G, a_down_new);

    return SUCCESS;
}

// starting from vertex linear (in-memory) index a, flow downstream and increment nodes by inc until reaching the root node
// visitationinc tracks how many times this function has been called since .vistitations was reset to 0. Can be either 1 or 2
Status sg_increment_downstream_safe(drainedarea_t inc, StreamGraph *G, linidx_t a, uint8_t ncalls){
    StreamGraph Gcopy = *G;

    Vertex vert;
    Status code;

    code = sg_get_lin_safe(Gcopy, &vert, a);
    if (code != SUCCESS) return code;
    while (vert.downstream != IS_ROOT){
        // if we find ourselves in a cycle, exit immediately and signal to the caller
        if (vert.visited == ncalls) return SWAP_WARNING;
        vert.visited += ncalls;

        G->energy -= pow((double)vert.drained_area, gamma);  // remove old energy contribution
        vert.drained_area += inc;
        G->energy += pow((double)vert.drained_area, gamma);  // add new energy contribution

        // any other errorcode: exit
        code = sg_get_lin_safe(Gcopy, &vert, vert.adown);
        if (code != SUCCESS) return code;
    }

    return SUCCESS;
}
Status sg_increment_downstream(drainedarea_t inc, StreamGraph *G, linidx_t a, uint8_t ncalls){
    StreamGraph Gcopy = *G;
    Vertex vert;

    sg_get_lin(Gcopy, &vert, a);
    while (vert.downstream != IS_ROOT){
        // if we find ourselves in a cycle, exit immediately and signal to the caller
        if (vert.visited == ncalls) return SWAP_WARNING;
        vert.visited += ncalls;
        
        G->energy -= pow((double)vert.drained_area, gamma);  // remove old energy contribution
        vert.drained_area += inc;  // update drained area
        G->energy += pow((double)vert.drained_area, gamma);  // add new energy contribution
        
        // any other errorcode: exit
        sg_get_lin(Gcopy, &vert, vert.adown);
    }

    return SUCCESS;
}

Status sg_single_erosion_event(StreamGraph *G){
    Status code;
    StreamGraph Gcopy = *G;

    Vertex vert, vert_down_old, vert_down_new;
    clockhand_t down_old, down_new;
    linidx_t a, a_down_old, a_down_new;
    
    drainedarea_t inc;
    double energy_old = G->energy;
    double energy_new;
    bool accept_bad_value = (((double)rand() / (double)RAND_MAX) / temperature) > 1.0;  // Metro-Hastings criterion

    linidx_t i;
    uint8_t ntries;

    while (true){
        // pick a random vertex and a random new downstream direction,
        a = rand() % (G->m * G->n);
        sg_get_lin(Gcopy, &vert, a);
        
        down_old = vert.downstream;
        
        a_down_old = vert.adown;
        
        inc = vert.drained_area;
        
        // try to find a new downstream direction that doesn't break the graph
        for (ntries = 0, down_new = rand() % 8; ntries < 8; ntries++, down_new++){
            code = sg_change_vertex_outflow(a, G, down_new);
            if (code != SUCCESS) continue;
            
            // this repeats code that was done in sg_change_vertex_outflow
            // TODO: refactor to avoid redundant calculation

            a_down_new = sg_cart_to_lin((cartidx_t)(sg_lin_to_cart_row(a_down_new, Gcopy.m, Gcopy.n) + row_offsets[down_new]), (cartidx_t)(sg_lin_to_cart_col(a_down_new, Gcopy.m, Gcopy.n) + col_offsets[down_new]), Gcopy.m, Gcopy.n);
            code = sg_get_lin_safe(Gcopy, &vert_down_old, a_down_old);
            if (code == OOB_ERROR){
                sg_change_vertex_outflow(a, G, down_old);  // undo the outflow change
                continue;
            }
            code = sg_get_lin_safe(Gcopy, &vert_down_new, a_down_new);
            if (code == OOB_ERROR){
                sg_change_vertex_outflow(a, G, down_old);  // undo the outflow change
                continue;
            }

            // zero the cycle tracker before updating drained area
            for (i = 0; i < Gcopy.m * Gcopy.n; i++) G->vertices[i].visited = 0;  
            // this also repeats code that was done in sg_change_vertex_outflow (get vertex)
            // TODO: refactor to avoid redundant calculation
            code = sg_increment_downstream_safe(-inc, G, a_down_old, 1);

            if (code == SWAP_WARNING){  // if we find ourselves in a cycle, retrace our steps and undo the increment
                for (i = 0; i < Gcopy.m * Gcopy.n; i++){
                    G->vertices[i].visited = 0;
                }
                sg_increment_downstream(inc, G, a_down_old, 1);  // undo the decrement
                sg_change_vertex_outflow(a, G, down_old);  // undo the outflow change
                continue;
            }

            // use a new visitation increment to avoid collision with the previous sg_increment_downstream call
            // TODO: find a way to use visitation increments that don't collide to avoid resetting the entire visited array each time
            code = sg_increment_downstream_safe(inc, G, a_down_new, 2);
            if (code == SWAP_WARNING){
                for (i = 0; i < Gcopy.m * Gcopy.n; i++){
                    G->vertices[i].visited = 0;
                }
                sg_increment_downstream(inc, G, a_down_old, 1);  // undo the decrement
                sg_increment_downstream(-inc, G, a_down_new, 2);  // undo the increment
                sg_change_vertex_outflow(a, G, down_old);  // undo the outflow change
                continue;
            }
            
            energy_new = G->energy;
            break;  // if we reached here, the swap was successful
        }
        if (code != SUCCESS) continue;  // if we exhausted all 8 directions without a successful swap, try again

        // if we reach here, the swap was successful, now we choose whether to accept it
        if (energy_new < energy_old || accept_bad_value){  // accept swap
            break;
        }

        // reject swap: undo everything
        for (i = 0; i < Gcopy.m * Gcopy.n; i++){
            G->vertices[i].visited = 0;
        }
        sg_increment_downstream(inc, G, a_down_old, 1);  // undo the decrement
        sg_increment_downstream(-inc, G, a_down_new, 2);  // undo the increment
        sg_change_vertex_outflow(a, G, a_down_old);  // undo the outflow change
    }

    return SUCCESS;
}

Status sg_outer_ocn_loop(uint32_t niterations, StreamGraph *G){
    Status code;
    uint32_t i = 0;
    for (i = 0; i < niterations; i++){
        code = sg_single_erosion_event(G);
        if (code != SUCCESS) exit(code);
    }
    return SUCCESS;
}

const char RIGHT_ARROW = '-';
const char DOWN_ARROW = '|';
const char LEFT_ARROW = '-';
const char UP_ARROW = '|';
const char DOWNRIGHT_ARROW = '\\';
const char DOWNLEFT_ARROW = '/';
const char UPLEFT_ARROW = '\\';
const char UPRIGHT_ARROW = '/';
const char NO_ARROW = ' ';
const char NODE = 'O';
const char ROOT_NODE = 'X';
void sg_display(StreamGraph G){
    Vertex vert;
    clockhand_t direction;
    for (cartidx_t printed_row = 0; printed_row < G.m*2; printed_row++){
        cartidx_t row = printed_row / 2;
        for (cartidx_t col = 0; col < G.n; col++){
            sg_get_lin(G, &vert, sg_cart_to_lin(row, col, G.m, G.n));
            // same row as vertex: horizontal edges
            if (printed_row % 2 == 0){
                if (vert.downstream == IS_ROOT) printf("%c", ROOT_NODE);
                else printf("%c", NODE);
                direction = (1 << 2);
                if (vert.edges & direction) printf("%c", RIGHT_ARROW);
                else printf("%c", NO_ARROW);
            // row between vertices: vertical and diagonal edges
            } else {
                direction = (1 << 4);
                if (vert.edges & direction){
                    printf("%c", DOWN_ARROW);
                } else {
                    printf("%c", NO_ARROW);
                }
                direction = (1 << 3);
                if (vert.edges & direction){
                    printf("%c", DOWNRIGHT_ARROW);
                } else {
                    printf("%c", NO_ARROW);
                }
            }
        }
        printf("\n");
    }
}

#endif // STREAMGRAPH_C