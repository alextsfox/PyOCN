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

#include "ocn.h"
#include "status.h"
#include "streamgraph.h"

clockhand_t IS_ROOT = 255;

// ##############################
// # Helpers Objects
// ##############################
const int16_t row_offsets[8] = {-1, -1, 0, 1, 1, 1, 0, -1};
const int16_t col_offsets[8] = {0, 1, 1, 1, 0, -1, -1, -1};

typedef struct {
    int16_t row;
    int16_t col;
} OffsetPair;
const OffsetPair offsets[8] = {
    {-1,  0}, // N
    {-1,  1}, // NE
    { 0,  1}, // E
    { 1,  1}, // SE
    { 1,  0}, // S
    { 1, -1}, // SW
    { 0, -1}, // W
    {-1, -1}  // NW
};

// 2d raster is block-tiled in memory to improve cache locality
/*
 * Converts cartesian (row, col) to linear index in block-tiled memory
 */
// const cartidx_t TILE_SIZE = 2;
// linidx_t sg_cart_to_lin(CartPair coords, CartPair dims){
//     cartidx_t row = coords.row;
//     cartidx_t col = coords.col;
//     cartidx_t n = dims.col;
//     cartidx_t m = dims.row;
//     div_t r_divT = div(row, TILE_SIZE);
//     div_t c_divT = div(col, TILE_SIZE);
//     linidx_t k = r_divT.quot * n/TILE_SIZE + c_divT.quot;
//     linidx_t a = (TILE_SIZE * TILE_SIZE * k) + (TILE_SIZE * r_divT.rem) + c_divT.rem;
//     return a;
// }
// for now, we just use normal row-major order until we finish debugging the main algorithm
linidx_t sg_cart_to_lin(CartPair coords, CartPair dims){
    return (coords.row * dims.col + coords.col);
}
CartPair sg_lin_to_cart(linidx_t a, CartPair dims){
    div_t adiv = div(a, dims.col);
    return (CartPair){adiv.quot, adiv.rem};
}
// CartPair sg_lin_to_cart(linidx_t a, CartPair dims){
//     cartidx_t p = a % TILE_SIZE;
//     cartidx_t q = (a / TILE_SIZE) % TILE_SIZE;
//     cartidx_t k = a / (TILE_SIZE * TILE_SIZE);

//     cartidx_t j = p + 2 * (k % (dims.col / TILE_SIZE));
//     cartidx_t i = q + 2 * (k / (dims.col / TILE_SIZE));
//     return (CartPair){i, j};
// }
Status sg_clockhand_to_lin_safe(linidx_t *a_down, linidx_t a, clockhand_t down, CartPair dims){
    CartPair row_col = sg_lin_to_cart(a, dims);
    OffsetPair offset = offsets[down];
    OffsetPair cart_down_off = {
        .row = (int16_t)row_col.row + offset.row,
        .col = (int16_t)row_col.col + offset.col
    };
    if (
        cart_down_off.row < 0 
        || cart_down_off.row >= (int16_t)dims.row 
        || cart_down_off.col < 0
        || cart_down_off.col >= (int16_t)dims.col
    ) return OOB_ERROR;
    CartPair cart_down = {
        .row = (cartidx_t)cart_down_off.row,
        .col = (cartidx_t)cart_down_off.col
    };

    *a_down = sg_cart_to_lin(cart_down, dims);

    return SUCCESS;
}

// ##############################
// # Getters + Setters          #
// ##############################
// cartesian
Status sg_get_cart_safe(Vertex *out, StreamGraph *G, CartPair coords){
    if (
        G == NULL 
        || out == NULL 
        || coords.row < 0 
        || coords.row >= G->dims.row 
        || coords.col < 0 
        || coords.col >= G->dims.col
    ) return OOB_ERROR;
    
    linidx_t a = sg_cart_to_lin(coords, G->dims);
    *out = G->vertices[a];
    return SUCCESS;
}
Vertex sg_get_cart(StreamGraph *G, CartPair coords){
    linidx_t a = sg_cart_to_lin(coords, G->dims);
    return G->vertices[a];
}
Status sg_set_cart_safe(StreamGraph *G, Vertex vert, CartPair coords){
    if (
        G == NULL 
        || coords.row < 0 
        || coords.row >= G->dims.row 
        || coords.col < 0 
        || coords.col >= G->dims.col
    ) return OOB_ERROR;
    linidx_t a = sg_cart_to_lin(coords, G->dims);
    G->vertices[a] = vert;
    return SUCCESS;
}
void sg_set_cart(StreamGraph *G, Vertex vert, CartPair coords){
    linidx_t a = sg_cart_to_lin(coords, G->dims);
    G->vertices[a] = vert;
}

// linear
Status sg_get_lin_safe(Vertex *out, StreamGraph *G, linidx_t a){
    if (G == NULL || out == NULL || a < 0 || a >= (G->dims.row * G->dims.col)) return OOB_ERROR;
    *out = G->vertices[a];
    return SUCCESS;
}
Vertex sg_get_lin(StreamGraph *G, linidx_t a){
    return G->vertices[a];
}
Status sg_set_lin_safe(StreamGraph *G, Vertex vert, linidx_t a){
    if (G == NULL || a < 0 || a >= (G->dims.row * G->dims.col)) return OOB_ERROR;
    G->vertices[a] = vert;
    return SUCCESS;
}
// Gets the value at [a] in memory to vert, with safeguards.
void sg_set_lin(StreamGraph *G, Vertex vert, linidx_t a){
    G->vertices[a] = vert;
}

// ##############################
// # Create/destroy streamgraph #
// ##############################
StreamGraph *sg_create_empty_safe(CartPair dims){
    if (dims.row % 2 != 0 || dims.col % 2 != 0) return NULL;  // dimensions must be even

    StreamGraph *G = malloc(sizeof(StreamGraph));
    if (G == NULL) return NULL;

    Vertex *vertices = malloc(dims.row * dims.col * sizeof(Vertex));
    if (vertices == NULL) {
        free(G);
        return NULL;
    }

    G->dims = dims;
    G->vertices = vertices;
    
    return G;
}

Status sg_destroy_safe(StreamGraph *G){
    if (G != NULL){
        if (G->vertices != NULL) free(G->vertices); G->vertices = NULL;
        free(G); 
    }
    return SUCCESS;
}


// ##############################
// # Network manipulation       #
// ##############################
Status sg_change_vertex_outflow(StreamGraph *G, linidx_t a, clockhand_t down_new){
    Status code;
    CartPair dims = G->dims;
    Vertex vert, vert_down_old, vert_down_new;
    linidx_t a_down_old, a_down_new;
    clockhand_t down_old;
    
    // 1. Get G[a], G[a_down_old], G[adownnew] safely
    code = sg_get_lin_safe(&vert, G, a);
    down_old = vert.downstream;
    if (code == OOB_ERROR) return OOB_ERROR;

    a_down_old = vert.adown;
    code = sg_get_lin_safe(&vert_down_old, G, a_down_old);
    if (code == OOB_ERROR) return OOB_ERROR;

    // a_down_new is trickier to get.
    code = sg_clockhand_to_lin_safe(&a_down_new, a, down_new, dims);
    if (code == OOB_ERROR) return OOB_ERROR;
    vert_down_new = sg_get_lin(G, a_down_new);  // we can use unsafe here because we already checked bounds

    // 2. check for any immediate problems with the swap that would malform the graph
    // check that the new downstream is valid (does not check for large cycles or for root access)
    if (
        (down_old == IS_ROOT)  // cannot rerout root node.
        || (down_new == down_old)  // no change
        || ((1u << down_new) & (vert.edges)) // new downstream direction already occupied
    ) return SWAP_WARNING;
    
    // check that we haven't created any crosses
    Vertex cross_check_vert;
    CartPair check_row_col = sg_lin_to_cart(a, dims);
    switch (down_new){
        case 1: check_row_col.row -= 1; break;  // NE flow: check N vertex
        case 7: check_row_col.row -= 1; break;  // NW flow: check N vertex
        case 3: check_row_col.row += 1; break;  // SE flow: check S vertex
        case 5: check_row_col.row += 1; break;  // SW flow: check S vertex
    }
    sg_get_cart_safe(&cross_check_vert, G, check_row_col);
    switch (down_new){
        case 1: if (cross_check_vert.edges & (1u << 3)) return SWAP_WARNING; break;  // NE flow: N vertex cannot have a SE edge
        case 7: if (cross_check_vert.edges & (1u << 5)) return SWAP_WARNING; break;  // NW flow: N vertex cannot have a SW edge
        case 3: if (cross_check_vert.edges & (1u << 1)) return SWAP_WARNING; break;  // SE flow: S vertex cannot have a NE edge
        case 5: if (cross_check_vert.edges & (1u << 7)) return SWAP_WARNING; break;  // SW flow: S vertex cannot have a NW edge
    }

    // 3. make the swap
    vert.adown = a_down_new;
    vert.downstream = down_new;
    vert.edges ^= ((1u << down_old) | (1u << down_new));
    sg_set_lin(G, vert, a);

    vert_down_old.edges ^= (1u << ((down_old + 4)%8));  // old downstream node loses that edge. Note that (down_old + 4) wraps around correctly because down_old is in [0,7]
    sg_set_lin(G, vert_down_old, a_down_old);

    vert_down_new.edges ^= (1u << ((down_new + 4)%8));  // new downstream node gains that edge
    sg_set_lin(G, vert_down_new, a_down_new);

    return SUCCESS;
}



Status sg_flow_downstream_safe(StreamGraph *G, linidx_t a, uint8_t ncalls){
    Vertex vert;
    Status code;
    code = sg_get_lin_safe(&vert, G, a);
    if (code != SUCCESS) return code;

    while (vert.downstream != IS_ROOT){
        // if we find ourselves in a cycle, exit immediately and signal to the caller
        if (vert.visited == ncalls) return MALFORMED_GRAPH_WARNING;
        
        vert.visited += ncalls;
        sg_set_lin(G, vert, a);  // unsafe is ok here because we already checked bounds

        // get next vertex
        a = vert.adown;
        code = sg_get_lin_safe(&vert, G, a);
        if (code != SUCCESS) return code;
    }
    return SUCCESS;  // found root successfully, no cycles found
}

// vibe-coded display function
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

void sg_display_bad(StreamGraph *G, bool use_utf8){
    if (use_utf8) setlocale(LC_ALL, "");

    if (G->vertices == NULL || G->dims.row == 0 || G->dims.col == 0) return;

    putchar('\n');

    cartidx_t m = G->dims.row;
    cartidx_t n = G->dims.col;

    for (cartidx_t i = 0; i < m; i++){
        // node row
        for (cartidx_t j = 0; j < n; j++){
            Vertex *v = &(G->vertices[i*n + j]);

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
            Vertex *v = &(G->vertices[i*n + j]);

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
                Vertex *v_right = &(G->vertices[i*n + (j+1)]);
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

void sg_display(StreamGraph *G, bool use_utf8){
    if (G->vertices == NULL || G->dims.row == 0 || G->dims.col == 0) return;
    
    if (use_utf8) setlocale(LC_ALL, "");
    putchar('\n');

    cartidx_t m = G->dims.row;
    cartidx_t n = G->dims.col;
    linidx_t size = m * n;

    for (cartidx_t i = 0; i < m; i++){
        // Node row
        for (cartidx_t j = 0; j < n; j++){
            linidx_t idx = i*n + j;
            Vertex *v = &(G->vertices[idx]);
            
            // Draw node
            bool is_root = (v->downstream == IS_ROOT || v->downstream == 255);
            if (use_utf8) {
                wprintf(L"%lc", is_root ? ROOT_NODE_UTF8 : NODE_UTF8);
            } else {
                putchar(is_root ? ROOT_NODE : NODE);
            }
            
            // Space after node
            putchar(' ');
            
            // Draw horizontal connector if not at right edge
            if (j < n - 1){
                linidx_t right_idx = i*n + (j+1);
                if (right_idx < size) {
                    Vertex *v_right = &(G->vertices[right_idx]);
                    
                    // Check for horizontal edge in either direction
                    bool has_east_edge = (v->edges & (1u << 2));
                    bool has_west_edge = (v_right->edges & (1u << 6));
                    
                    if (has_east_edge || has_west_edge) {
                        // Determine direction based on downstream
                        if (v->downstream == 2) {
                            if (use_utf8) wprintf(L"%lc", E_ARROW_UTF8);
                            else putchar(E_ARROW);
                        } else if (v_right->downstream == 6) {
                            if (use_utf8) wprintf(L"%lc", W_ARROW_UTF8);
                            else putchar(W_ARROW);
                        } else {
                            // No flow direction specified, use a neutral connector
                            if (use_utf8) wprintf(L"%lc", E_ARROW_UTF8);
                            else putchar(E_ARROW);
                        }
                    } else {
                        // No edge
                        if (use_utf8) wprintf(L"%lc", NO_ARROW_UTF8);
                        else putchar(NO_ARROW);
                    }
                }
                
                // Space after horizontal connector
                putchar(' ');
            }
        }
        putchar('\n');
        
        // Skip vertical connectors on last row
        if (i == m - 1) break;
        
        // Vertical/diagonal connector row
        for (cartidx_t j = 0; j < n; j++){
            linidx_t idx = i*n + j;
            linidx_t below_idx = (i+1)*n + j;
            Vertex *v = &(G->vertices[idx]);
            Vertex *v_below = (below_idx < size) ? &(G->vertices[below_idx]) : NULL;
            
            // Draw vertical connector
            bool has_south_edge = (v->edges & (1u << 4));
            bool has_north_edge = (v_below && (v_below->edges & (1u << 0)));
            
            if (has_south_edge || has_north_edge) {
                if (v->downstream == 4) {
                    if (use_utf8) wprintf(L"%lc", S_ARROW_UTF8);
                    else putchar(S_ARROW);
                } else if (v_below && v_below->downstream == 0) {
                    if (use_utf8) wprintf(L"%lc", N_ARROW_UTF8);
                    else putchar(N_ARROW);
                } else {
                    // Default direction
                    if (use_utf8) wprintf(L"%lc", S_ARROW_UTF8);
                    else putchar(S_ARROW);
                }
            } else {
                if (use_utf8) wprintf(L"%lc", NO_ARROW_UTF8);
                else putchar(NO_ARROW);
            }
            
            putchar(' ');
            
            // Draw diagonal connector if not at right edge
            if (j < n - 1) {
                linidx_t right_idx = i*n + (j+1);
                linidx_t diag_idx = (i+1)*n + (j+1);
                Vertex *v_right = (right_idx < size) ? &(G->vertices[right_idx]) : NULL;
                Vertex *v_diag = (diag_idx < size) ? &(G->vertices[diag_idx]) : NULL;
                
                // SE/NW diagonal
                bool has_se_edge = (v->edges & (1u << 3));
                bool has_nw_edge = (v_diag && (v_diag->edges & (1u << 7)));
                
                // SW/NE diagonal
                bool has_sw_edge = (v_right && (v_right->edges & (1u << 5)));
                bool has_ne_edge = (v_below && (v_below->edges & (1u << 1)));
                
                char mid = NO_ARROW;
                wchar_t mid_utf8 = NO_ARROW_UTF8;
                
                // Prioritize SE/NW diagonal
                if (has_se_edge || has_nw_edge) {
                    if (v->downstream == 3) {
                        mid = SE_ARROW;
                        mid_utf8 = SE_ARROW_UTF8;
                    } else if (v_diag && v_diag->downstream == 7) {
                        mid = NW_ARROW;
                        mid_utf8 = NW_ARROW_UTF8;
                    } else {
                        // Default
                        mid = SE_ARROW;
                        mid_utf8 = SE_ARROW_UTF8;
                    }
                }
                // Then check SW/NE diagonal
                else if (has_sw_edge || has_ne_edge) {
                    if (v_right && v_right->downstream == 5) {
                        mid = SW_ARROW;
                        mid_utf8 = SW_ARROW_UTF8;
                    } else if (v_below && v_below->downstream == 1) {
                        mid = NE_ARROW;
                        mid_utf8 = NE_ARROW_UTF8;
                    } else {
                        // Default
                        mid = SW_ARROW;
                        mid_utf8 = SW_ARROW_UTF8;
                    }
                }
                
                if (use_utf8) wprintf(L"%lc", mid_utf8);
                else putchar(mid);
                
                putchar(' ');
            }
        }
        putchar('\n');
    }
    putchar('\n');
}

#endif // STREAMGRAPH_C