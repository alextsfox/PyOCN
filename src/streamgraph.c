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
#include <string.h>
#include <stdio.h>
#include <wchar.h>
#include <locale.h>
// #include <locale.h>

typedef uint8_t Status;

const Status SUCCESS = 0;
// vertex return codes
// valid values for .downstream are 0-7
// use remaining values for error codes
const Status OOB_ERROR = 1;
const Status NO_EDGE_WARNING = 2;
const Status NULL_POINTER_ERROR = 3;
const Status SWAP_WARNING = 4;
const Status MALFORMED_GRAPH_WARNING = 5;
const Status CYCLE_WARNING = 6;

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
typedef uint32_t drainedarea_t;
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

struct cartesian_pair{cartidx_t i, j;};


// ##############################
// # Helpers
// ##############################
const int16_t row_offsets[8] = {-1, -1, 0, 1, 1, 1, 0, -1};
const int16_t col_offsets[8] = {0, 1, 1, 1, 0, -1, -1, -1};


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
Status sg_lin_to_cart(linidx_t a, struct cartesian_pair *out, cartidx_t m, cartidx_t n){
    for (cartidx_t row = 0; row < m; row++){
        for (cartidx_t col = 0; col < n; col++){
            if (a == sg_cart_to_lin(row, col, m, n)) {
                out->i = row;
                out->j = col;
                return SUCCESS;
            }
        }
    }
    return OOB_ERROR;
}

cartidx_t sg_lin_to_cart_col(linidx_t a, cartidx_t m, cartidx_t n){
    for (cartidx_t row = 0; row < m; row++){
        for (cartidx_t col = 0; col < n; col++){
            if (a == sg_cart_to_lin(row, col, m, n)) return col;
        }
    }
    return 0;
}
cartidx_t sg_lin_to_cart_row(linidx_t a, cartidx_t m, cartidx_t n){
    for (cartidx_t row = 0; row < m; row++){
        for (cartidx_t col = 0; col < n; col++){
            if (a == sg_cart_to_lin(row, col, m, n)) return row;
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
// # Create/destroy streamgraph #
// ##############################
Status sg_create(StreamGraph *G, cartidx_t m, cartidx_t n, cartidx_t i_root, cartidx_t j_root){
    if (m % 2 != 0 || n % 2 != 0) return MALFORMED_GRAPH_WARNING;
    if (G == NULL) return NULL_POINTER_ERROR;
    Vertex *vertices = calloc(m*n, sizeof(Vertex));
    if (vertices == NULL) return NULL_POINTER_ERROR;
    if (i_root < 0 || i_root >= m || j_root < 0 || j_root >= n) return OOB_ERROR;
    vertices[sg_cart_to_lin(i_root, j_root, m, n)].downstream = IS_ROOT;  // mark root
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

/**
 * @brief Get the downstream neighbor of a vertex in a specified direction, with safeguards.
 * @param m Number of rows in the grid
 * @param n Number of columns in the grid
 * @param vert_down Pointer to store the downstream neighbor vertex
 * @param a Linear index of the current vertex
 * @param down The downstream direction (0-7)
 * @return SUCCESS or OOB_ERROR if out of bounds
 */
Status sg_clockhand_to_lin(cartidx_t m, cartidx_t n, linidx_t *a_down, linidx_t a, clockhand_t down){
    struct cartesian_pair row_col;
    Status code;

    code = sg_lin_to_cart(a, &row_col, m, n);
    if (code != SUCCESS) return code;

    int16_t col_off = col_offsets[down];
    int16_t col = (int16_t)row_col.j + col_off;  // cast to signed int to avoid underflow
    if (col < 0 || col >= (int16_t)n) return OOB_ERROR;

    int16_t row_off = row_offsets[down];
    int16_t row = (int16_t)row_col.i + row_off;
    if (row < 0 || row >= (int16_t)m) return OOB_ERROR;

    *a_down = sg_cart_to_lin((cartidx_t)row, (cartidx_t)col, m, n);  // recast back to unsigned before converting to linear index

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
    cartidx_t m = G->m;
    cartidx_t n = G->n;

    Status code;

    Vertex vert;
    Vertex vert_down_old;
    Vertex vert_down_new;
    
    // retrieve G[a], G[a_down_old], G[adownnew]
    // return OOB_ERROR if any of these fail
    code = sg_get_lin_safe(*G, &vert, a);
    if (code == OOB_ERROR) return OOB_ERROR;

    linidx_t a_down_old = vert.adown;
    code = sg_get_lin_safe(*G, &vert_down_old, a_down_old);
    if (code == OOB_ERROR) return OOB_ERROR;

    
    // int16_t col_off = col_offsets[down_new];
    // int16_t col = (int16_t)row_col.j + col_off;  // cast to signed int to avoid underflow
    // if (col < 0 || col >= (int16_t)Gcopy.n) return OOB_ERROR;
    
    // int16_t row_off = row_offsets[down_new];
    // int16_t row = (int16_t)row_col.i + row_off;
    // if (row < 0 || row >= (int16_t)Gcopy.m) return OOB_ERROR;
    
    // linidx_t a_down_new = sg_cart_to_lin((cartidx_t)row, (cartidx_t)col, Gcopy.m, Gcopy.n);  // recast back to unsigned before converting to linear index
    linidx_t a_down_new;
    code = sg_clockhand_to_lin(m, n, &a_down_new, a, down_new);
    if (code == OOB_ERROR) return OOB_ERROR;
    code = sg_get_lin_safe(*G, &vert_down_new, a_down_new);
    if (code == OOB_ERROR) return OOB_ERROR;
    
    // 2. 
    // check that the new downstream is valid (does not check for large cycles or for root access)
    if (
        (vert.downstream == IS_ROOT)  // root node
        || (down_new == vert.downstream)  // no change
        || ((1u << down_new) & (vert.edges)) // new downstream not a valid edge
    ) return SWAP_WARNING;
    
    // check that we haven't created any crosses
    // how to check: upper right quadrant (new_downstream = 1): vertically upwards vertex cannot have an edge at 2 (edges & (1 << 2) == 0)
    // upper left quadrant (new_downstream = 7): vertically upwards vertex cannot have an edge at 5 (edges & (1 << 5) == 0)
    // lower right quadrant (new_downstream = 3): vertically downwards vertex cannot have an edge at 1 (edges & (1 << 1) == 0)
    // lower left quadrant (new_downstream = 5): vertically downwards vertex cannot have an edge at 6 (edges & (1 << 6) == 0)
    Vertex cross_check_vert;
    struct cartesian_pair row_col;
    code = sg_lin_to_cart(a, &row_col, m, n);
    if (code == OOB_ERROR) return OOB_ERROR;
    cartidx_t check_row = (cartidx_t)row_col.i;
    cartidx_t check_col = (cartidx_t)row_col.j;
    switch (down_new){
        case 1: check_row -= 1; break;
        case 3: check_row += 1; break;
        case 5: check_col -= 1; break;
        case 7: check_col += 1; break;
    }
    sg_get_cart_safe(*G, &cross_check_vert, check_row, check_col);
    switch (down_new){
        case 1: if (cross_check_vert.edges & (1u << 2)) return SWAP_WARNING; break;
        case 3: if (cross_check_vert.edges & (1u << 1)) return SWAP_WARNING; break;
        case 5: if (cross_check_vert.edges & (1u << 6)) return SWAP_WARNING; break;
        case 7: if (cross_check_vert.edges & (1u << 5)) return SWAP_WARNING; break;
    }

    // 3. make the swap
    vert.adown = a_down_new;
    clockhand_t down_old = vert.downstream;
    vert.downstream = down_new;
    vert.edges ^= ((1u << down_old) | (1u << down_new));
    sg_set_lin(vert, G, a);

    vert_down_old.edges ^= (1u << (down_old + 4));
    sg_set_lin(vert_down_old, G, a_down_old);
    vert_down_new.edges ^= (1u << (down_new + 4));
    sg_set_lin(vert_down_new, G, a_down_new);

    return SUCCESS;
}

// simply flows downstream. Does not accumulate any values other than a visit checker
// returns a MALFORMED_GRAPH_WARNING if a cycle is detected, an OOB_ERROR is thrown, or if we cannot find the root node.
Status sg_flow_downstream_safe(StreamGraph *G, linidx_t a, uint8_t ncalls){
    Vertex vert;
    Status code;
    code = sg_get_lin_safe(*G, &vert, a);
    if (code != SUCCESS) return code;

    while (vert.downstream != IS_ROOT){
        // if we find ourselves in a cycle, exit immediately and signal to the caller
        if (vert.visited == ncalls) return CYCLE_WARNING;
        vert.visited += ncalls;
        sg_set_lin(vert, G, a);

        // any other errorcode: exit
        a = vert.adown;
        code = sg_get_lin_safe(*G, &vert, vert.adown);
        if (code != SUCCESS) return code;

    }

    return SUCCESS;  // found root successfully, no cycles found
}

// flows downstream, incrementing drained area of each node, and updating the energy
// UNSAFE; assumes that the graph is well-formed and acyclic. Does no error checking.
Status sg_increment_downstream(drainedarea_t da_inc, StreamGraph *G, linidx_t a){
    StreamGraph Gcopy = *G;
    Vertex vert;
    drainedarea_t da;
    sg_get_lin(Gcopy, &vert, a);
    while (vert.downstream != IS_ROOT){
        da = vert.drained_area;
        G->energy += pow((double)(da + da_inc), gamma) - pow((double)da, gamma);  // update energy
        vert.drained_area += da_inc;  // update drained area of vertex
        sg_set_lin(vert, G, a);
        a = vert.adown;
        sg_get_lin(Gcopy, &vert, vert.adown);  // get next
    }
    return SUCCESS;
}

Status sg_single_erosion_event(StreamGraph *G){
    Status code;

    Vertex vert, vert_down_old, vert_down_new;
    clockhand_t down_old, down_new;
    linidx_t a, a_down_old, a_down_new;
    
    cartidx_t m, n;
    m = G->m;
    n = G->n;

    drainedarea_t da_inc;
    double energy_old = G->energy;
    double energy_new;
    bool accept_bad_value = (((double)rand() / (double)RAND_MAX) / temperature) > 1.0;  // Metro-Hastings criterion

    linidx_t i;
    uint8_t ntries;
    for (uint16_t total_tries = 0; total_tries < 1000; total_tries++){
        // pick a random vertex and a random new downstream direction,
        down_new = rand() % 8;
        a = rand() % (m*n);
        sg_get_lin(*G, &vert, a);
        
        down_old = vert.downstream;
        a_down_old = vert.adown;
        da_inc = vert.drained_area;

        for (ntries = 0; ntries < 8; ntries++){  // try a new direction each time, up to 8 times before giving up and picking a new vertex
            down_new  = (down_new + 1) % 8;

            code = sg_change_vertex_outflow(a, G, down_new);
            if (code != SUCCESS){
                continue;  // try again if the swap is immediately invalid
            }
            
            // retrieve the downstream vertices...shouldn't fail, since we already checked bounds in sg_change_vertex_outflow
            a_down_new = vert.adown;
            sg_get_lin(*G, &vert_down_new, a_down_new);
            sg_get_lin(*G, &vert_down_old, a_down_old);

            // zero the cycle tracker before updating drained area
            for (i = 0; i < m*n; i++) G->vertices[i].visited = 0;
            // this also repeats code that was done in sg_change_vertex_outflow (get vertex)
            // TODO: refactor to avoid redundant calculation
            code = sg_flow_downstream_safe(G, a_down_old, 1);
            if (code != SUCCESS){
                sg_change_vertex_outflow(a, G, down_old);  // undo the swap, try again
                continue;
            }

            // use a new visitation increment to avoid collision with the previous sg_increment_downstream call
            // TODO: find a way to use visitation increments that don't collide to avoid resetting the entire visited array each time
            code = sg_flow_downstream_safe(G, a_down_new, 2);
            if (code != SUCCESS){
                sg_change_vertex_outflow(a, G, down_old);  // undo the swap, try again
                continue;
            }

            break;  // if we reached here, the swap resulted in a well-formed graph, so we can move on the acceptance step
        }
        if (code != SUCCESS) continue;  // if we exhausted all 8 directions without a successful swap, try again. do not go to acceptance step

        // update drained area and energy along both paths
        energy_old = G->energy;
        sg_increment_downstream(-da_inc, G, a_down_old);  // decrement drained area along old path
        sg_increment_downstream(da_inc, G, a_down_new);  // increment drained area along new path
        energy_new = G->energy;

        if (energy_new < energy_old || accept_bad_value){  // accept the swap
            return SUCCESS;
        }

        // reject swap: undo everything and try again
        for (i = 0; i < G->m * G->n; i++){
            G->vertices[i].visited = 0;
        }
        sg_increment_downstream(da_inc, G, a_down_old);  // undo the decrement
        sg_increment_downstream(-da_inc, G, a_down_new);  // undo the increment
        sg_change_vertex_outflow(a, G, a_down_old);  // undo the outflow change
    }

    return MALFORMED_GRAPH_WARNING;  // if we reach here, we failed to find a valid swap in many, many tries
}

Status sg_outer_ocn_loop(uint32_t niterations, StreamGraph *G){
    Status code;
    uint32_t i = 0;
    for (i = 0; i < niterations; i++){
        code = sg_single_erosion_event(G);
        if (code != SUCCESS) return code;  // malformed graph?
    }
    return SUCCESS;
}

// vibe-coded slop, I guess it works.
const char E_ARROW = '-';
const char S_ARROW = '|';
const char W_ARROW = '-';
const char N_ARROW = '|';
const char SE_ARROW = '\\';
const char SW_ARROW = '/';
const char NW_ARROW = '\\';
const char NE_ARROW = '/';
const char NO_ARROW = ' ';
const char NODE = 'O';
const char ROOT_NODE = 'X';


const wchar_t E_ARROW_UTF8 = L'\u2192';
const wchar_t S_ARROW_UTF8 = L'\u2193';
const wchar_t W_ARROW_UTF8 = L'\u2190';
const wchar_t N_ARROW_UTF8 = L'\u2191';
const wchar_t SE_ARROW_UTF8 = L'\u2198';
const wchar_t SW_ARROW_UTF8 = L'\u2199';
const wchar_t NW_ARROW_UTF8 = L'\u2196';
const wchar_t NE_ARROW_UTF8 = L'\u2197';
const wchar_t NO_ARROW_UTF8 = L'\u2002';  // space
const wchar_t NODE_UTF8 = L'\u25EF';  // large empty circle
const wchar_t ROOT_NODE_UTF8 = L'\u25CF';  // circle with x

/**
 * @brief Display the stream graph in a human-readable format.
 * Nodes are represented by 'O' (or 'X' for the root), and edges are represented by arrows.
 * @param G The StreamGraph to display.
 */
void sg_display(StreamGraph G, bool use_utf8){
    if (use_utf8) setlocale(LC_ALL, "");

    if (G.vertices == NULL || G.m == 0 || G.n == 0) return;

    putchar('\n');

    cartidx_t m = G.m;
    cartidx_t n = G.n;

    for (cartidx_t i = 0; i < m; i++){
        // node row
        for (cartidx_t j = 0; j < n; j++){
            Vertex *v = &G.vertices[i*n + j];

            // Draw node
            if (use_utf8) wprintf(L"%lc", v->downstream == IS_ROOT ? ROOT_NODE : NODE);
            else putchar(v->downstream == IS_ROOT ? ROOT_NODE : NODE);
            if (use_utf8) putchar (' ');
            
            if (j < n - 1){
                // Draw horizontal edge going right out of node
                if (v->edges & (1u << 2)){        
                    if (v->downstream == 2) {
                        if (use_utf8) wprintf(L"%lc", E_ARROW_UTF8);        
                        else putchar(E_ARROW);
                    } else {
                        if (use_utf8) wprintf(L"%lc", W_ARROW_UTF8);    
                        else putchar(W_ARROW);    
                    }
                } else {  // no horizontal edge edge
                    if (use_utf8) putchar(NO_ARROW);
                    else wprintf(L"%lc", NO_ARROW_UTF8);
                }
                if (use_utf8) putchar(' ');
            }
        }
        putchar('\n');

        if (i == m - 1) break; // last row: no edges below

        // edge-only row
        for (cartidx_t j = 0; j < n; j++){
            Vertex *v = &G.vertices[i*n + j];

            if (v->edges & (1u << 4)){  // vertical, under node
                if (v->downstream == 4){
                    if (use_utf8) wprintf(L"%lc", S_ARROW_UTF8);    
                    else putchar(S_ARROW);
                } else {
                    if (use_utf8) wprintf(L"%lc", N_ARROW_UTF8);
                    else putchar(N_ARROW);
                }
            } else {  // no vertical under node
                if (use_utf8) wprintf(L"%lc", NO_ARROW_UTF8);
                else putchar(NO_ARROW);
            }
            if (use_utf8) putchar(' ');

            if (j < n - 1){// Diagonal edges occupy the space between two vertical edges
                wchar_t mid_utf8 = NO_ARROW_UTF8;
                char mid = NO_ARROW;  // default to no diagonal edge
                if (v->edges & (1u << 3)) {  // diagonal edge heading SE or NW
                    if (v->downstream == 3){
                        mid_utf8 = SE_ARROW_UTF8;
                        mid = SE_ARROW;
                    } else {
                        mid_utf8 = NW_ARROW_UTF8;
                        mid = NW_ARROW;
                    }
                }
                Vertex *v_right = &G.vertices[i*n + (j+1)];
                if (mid == NO_ARROW && (v_right->edges & (1u << 5))){  // diagonal edge heading NE or SW
                    if (v_right->downstream == 5) {
                        mid_utf8 = SW_ARROW_UTF8;
                        mid = SW_ARROW;
                    } else {
                        mid_utf8 = NE_ARROW_UTF8;
                        mid = NE_ARROW;
                    }   
                } 
                if (use_utf8) wprintf(L"%lc", mid_utf8);
                else putchar(mid);
                if (use_utf8) putchar(' ');
            }
        }
        putchar('\n');
    }
    putchar('\n');
}

StreamGraph sg_make_test_graph(){
    const cartidx_t m = 4;
    const cartidx_t n = 4;


    Vertex cart_vertices[] = {
        (Vertex){1, 0,  0b00010000, 4, 0}, (Vertex){1, 0, 0b00010000, 4, 0}, (Vertex){1, 0, 0b00010000, 4, 0}, (Vertex){1, 0, 0b00010000, 4, 0},
        (Vertex){2, 0, 0b00010001, 4, 0}, (Vertex){2, 0, 0b00010001, 4, 0}, (Vertex){2, 0, 0b00010001, 4, 0}, (Vertex){2, 0, 0b00010001, 4, 0},
        (Vertex){3, 0, 0b00010001, 4, 0}, (Vertex){3, 0, 0b00010001, 4, 0}, (Vertex){3, 0, 0b00010001, 4, 0}, (Vertex){3, 0, 0b00010001, 4, 0},
        (Vertex){4, 0, 0b00000101, 2, 0}, (Vertex){8, 0, 0b01000101, 2, 0}, (Vertex){12, 0,0b01000101, 2, 0}, (Vertex){16, 0,0b01000001, IS_ROOT, 0},
    };
    Vertex lin_vertices[m*n];
    Vertex vert;
    // awful hack until I write a proper conversion function i,j -> a
    // TODO write proper conversion function
    for (linidx_t a = 0; a < m*n; a++){
        for (cartidx_t row = 0; row < m; row++){
            for (cartidx_t col = 0; col < n; col++){
                if (a == sg_cart_to_lin(row, col, m, n)){
                    vert = cart_vertices[row*n + col];
                    vert.adown = sg_cart_to_lin(row + row_offsets[vert.downstream], col + col_offsets[vert.downstream], m, n);
                    lin_vertices[a] = vert;
                }
            }
        }
    }


    StreamGraph G;
    sg_create(&G, m, n, 3, 3);
    G.energy = 0.0;
    memcpy(G.vertices, lin_vertices, G.m * G.n * sizeof(Vertex));
    G.energy = 4 + sqrt(2)*4 + sqrt(3)*4 + sqrt(4)*4;  // pre-calculate initial energy

    return G;
}

#endif // STREAMGRAPH_C