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
// linidx_t sg_cart_to_lin(cartidx_t row, cartidx_t col, cartidx_t m, cartidx_t n){
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
Status sg_clockhand_to_lin_safe(linidx_t *a_down, linidx_t a, clockhand_t down, CartPair dims){
    Status code;

    CartPair row_col = sg_lin_to_cart(a, dims);
    OffsetPair offset = offsets[down];
    OffsetPair cart_down_off = {
        .row = (int16_t)row_col.row + offsets[down].row,
        .col = (int16_t)row_col.col + offsets[down].col
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
Status sg_create_empty_safe(StreamGraph *G, CartPair root, CartPair dims){
    if (dims.row % 2 != 0 || dims.col % 2 != 0) return MALFORMED_GRAPH_WARNING;
    if (G == NULL) return NULL_POINTER_ERROR;

    Vertex *vertices = malloc(dims.row * dims.col * sizeof(Vertex));
    if (vertices == NULL) return NULL_POINTER_ERROR;
    if (root.row < 0 || root.row >= dims.row || root.col < 0 || root.col >= dims.col) return OOB_ERROR;
    vertices[sg_cart_to_lin(root, dims)].downstream = IS_ROOT;

    G->dims = dims;
    G->root = root;
    G->energy = 0.0;
    G->vertices = vertices;

    return SUCCESS;
}
Status sg_destroy_safe(StreamGraph *G){
    if (G != NULL && G->vertices != NULL){
        free(G->vertices);
        G->vertices = NULL;
    }
    G = NULL;
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

StreamGraph *sg_make_test_graph(){
    const cartidx_t m = 4;
    const cartidx_t n = 4;
    const CartPair dims = {m, n};


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
                CartPair coords = {row, col};
                if (a == sg_cart_to_lin(coords, dims)){
                    vert = cart_vertices[row*n + col];
                    linidx_t a_down;
                    sg_clockhand_to_lin_safe(&a_down, a, vert.downstream, dims);
                    vert.adown = a_down;
                    lin_vertices[a] = vert;
                }
            }
        }
    }


    StreamGraph *G;
    sg_create_empty_safe(G, (CartPair){3, 3}, dims);
    memcpy(G->vertices, lin_vertices, (G->dims.row * G->dims.col) * sizeof(Vertex));
    G->energy = 4 + sqrt(2)*4 + sqrt(3)*4 + sqrt(4)*4;  // pre-calculate initial energy

    return G;
}





void sg_display(StreamGraph *G, bool use_utf8){
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
            if (use_utf8) wprintf(L"%lc", v->downstream == IS_ROOT ? ROOT_NODE_UTF8 : NODE_UTF8);
            else putchar(v->downstream == IS_ROOT ? ROOT_NODE : NODE);
            if (use_utf8) putchar(' ');
            
            if (j < n - 1){
                // Draw horizontal edge going right out of node
                // Check if there's an edge from current node to right node (E) or from right node to current node (W)
                Vertex *v_right = &(G->vertices[i*n + (j+1)]);
                bool has_edge = (v->edges & (1u << 2)) || (v_right->edges & (1u << 6));
                
                if (has_edge) {        
                    if (v->downstream == 2) {
                        // Current node flows east
                        if (use_utf8) wprintf(L"%lc", E_ARROW_UTF8);        
                        else putchar(E_ARROW);
                    } else if (v_right->downstream == 6) {
                        // Right node flows west
                        if (use_utf8) wprintf(L"%lc", W_ARROW_UTF8);    
                        else putchar(W_ARROW);    
                    } else {
                        // Non-flow edge, use default direction
                        if (use_utf8) wprintf(L"%lc", E_ARROW_UTF8);
                        else putchar(E_ARROW);
                    }
                } else {  // no horizontal edge
                    if (use_utf8) wprintf(L"%lc", NO_ARROW_UTF8);
                    else putchar(NO_ARROW);
                }
                if (use_utf8) putchar(' ');
            }
        }
        putchar('\n');

        if (i == m - 1) break; // last row: no edges below

        // edge-only row
        for (cartidx_t j = 0; j < n; j++){
            Vertex *v = &(G->vertices[i*n + j]);
            
            // Check both current node's South edge and bottom node's North edge
            bool has_v_edge = false;
            bool flows_south = false;
            bool flows_north = false;
            
            if (i < m - 1) {
                Vertex *v_bottom = &(G->vertices[(i+1)*n + j]);
                has_v_edge = (v->edges & (1u << 4)) || (v_bottom->edges & (1u << 0));
                flows_south = (v->downstream == 4);
                flows_north = (v_bottom->downstream == 0);
            }
            
            if (has_v_edge) {
                if (flows_south) {
                    if (use_utf8) wprintf(L"%lc", S_ARROW_UTF8);    
                    else putchar(S_ARROW);
                } else if (flows_north) {
                    if (use_utf8) wprintf(L"%lc", N_ARROW_UTF8);
                    else putchar(N_ARROW);
                } else {
                    // Non-flow edge, use default direction
                    if (use_utf8) wprintf(L"%lc", S_ARROW_UTF8);    
                    else putchar(S_ARROW);
                }
            } else {  // no vertical under node
                if (use_utf8) wprintf(L"%lc", NO_ARROW_UTF8);
                else putchar(NO_ARROW);
            }
            if (use_utf8) putchar(' ');

            if (j < n - 1) {
                // Diagonal edges occupy the space between two vertical edges
                wchar_t mid_utf8 = NO_ARROW_UTF8;
                char mid = NO_ARROW;  // default to no diagonal edge
                
                bool has_se_edge = false;
                bool has_ne_edge = false;
                bool flows_se = false;
                bool flows_nw = false;
                bool flows_sw = false;
                bool flows_ne = false;
                
                // Check SE/NW edge
                if (i < m - 1 && j < n - 1) {
                    Vertex *v_bottom_right = &(G->vertices[(i+1)*n + (j+1)]);
                    has_se_edge = (v->edges & (1u << 3)) || (v_bottom_right->edges & (1u << 7));
                    flows_se = (v->downstream == 3);
                    flows_nw = (v_bottom_right->downstream == 7);
                }
                
                // Check SW/NE edge
                Vertex *v_right = &(G->vertices[i*n + (j+1)]);
                if (i < m - 1) {
                    Vertex *v_bottom = &(G->vertices[(i+1)*n + j]);
                    has_ne_edge = (v_right->edges & (1u << 5)) || (v_bottom->edges & (1u << 1));
                    flows_sw = (v_right->downstream == 5);
                    flows_ne = (v_bottom->downstream == 1);
                }
                
                // Prioritize SE/NW edge
                if (has_se_edge) {
                    if (flows_se) {
                        mid_utf8 = SE_ARROW_UTF8;
                        mid = SE_ARROW;
                    } else if (flows_nw) {
                        mid_utf8 = NW_ARROW_UTF8;
                        mid = NW_ARROW;
                    } else {
                        // Non-flow edge, use default direction
                        mid_utf8 = SE_ARROW_UTF8;
                        mid = SE_ARROW;
                    }
                } else if (has_ne_edge) {  // If no SE/NW edge, check NE/SW
                    if (flows_sw) {
                        mid_utf8 = SW_ARROW_UTF8;
                        mid = SW_ARROW;
                    } else if (flows_ne) {
                        mid_utf8 = NE_ARROW_UTF8;
                        mid = NE_ARROW;
                    } else {
                        // Non-flow edge, use default direction
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









#endif // STREAMGRAPH_C