#include <stdio.h>
#include "streamgraph.c"
#include "returncodes.h"
#include <stdint.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <stdbool.h>
#include <assert.h>

int main(void){
/*
    O O O O
    | | | |
    O O O O
    | | | |
    O O O O
    | | | |
    O-O-O-X
    */

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
                    // printf("i:%d, j:%d, a:%d, adown:%d, edges:%d, downstream:%d\n", row, col, a, vert.adown, vert.edges, vert.downstream);
                }
            }
        }
    }


    StreamGraph G;
    sg_create(&G, m, n, 3, 3);
    G.energy = 0.0;
    G.vertices = lin_vertices;
    G.energy = 4 + sqrt(2)*4 + sqrt(3)*4 + sqrt(4)*4;  // pre-calculate initial energy

    sg_display(G);

    
    // Test getters
    sg_get_cart_safe(G, &vert, 3, 2);
    assert(vert.drained_area == 12);
    assert(vert.edges == 0b01000101);
    assert(vert.adown == 15);
    assert(vert.downstream == 2);
    assert(vert.visited == 0);
    
    sg_get_lin_safe(G, &vert, 14);
    assert(vert.drained_area == 12);
    assert(vert.edges == 0b01000101);
    assert(vert.adown == 15);
    assert(vert.downstream == 2);
    assert(vert.visited == 0);

    // test vertex outflow switch
    // no change
    assert(sg_change_vertex_outflow(14, &G, 2) == SWAP_WARNING);
    sg_get_lin_safe(G, &vert, 14);
    assert(vert.downstream == 2);
    assert(vert.adown == 15);
    assert(vert.edges == 0b01000101);

    // successful rerout
    // O O O O
    // | | | |
    // O O O O
    // | | | |
    // O O O O
    // | | |/|
    // O-O-O X
    assert(sg_change_vertex_outflow(14, &G, 1) == SUCCESS);
    sg_get_lin_safe(G, &vert, 14);
    assert(vert.downstream == 1);
    assert(vert.adown == 11);
    assert(vert.edges == 0b01000011);
    
    sg_get_lin_safe(G, &vert, 15);
    assert(vert.downstream == IS_ROOT);
    assert(vert.edges == 0b00000001);
    
    sg_get_lin_safe(G, &vert, 11);
    assert(vert.downstream == 4);
    assert(vert.edges == 0b00110001);

    sg_display(G);

    // attempt to create cross-flow
    assert(sg_change_vertex_outflow(10, &G, 3) == SWAP_WARNING);
    sg_get_lin_safe(G, &vert, 10);  // confirm that nothing changed
    assert(vert.downstream == 4);
    assert(vert.adown == 14);
    assert(vert.edges == 0b00010001);
    sg_get_lin_safe(G, &vert, 15);
    assert(vert.downstream == IS_ROOT);
    assert(vert.edges == 0b00000001);
    sg_get_lin_safe(G, &vert, 14);
    assert(vert.downstream == 1);
    assert(vert.edges == 0b01000011);

    // successful rerout
    // O O O O
    // | | | |
    // O O O O
    // | | | |
    // O O O O
    // | | |/|
    // O-O-O X
    sg_display(G);

    // flow downstream from a vertex
    assert(sg_flow_downstream_safe(&G, 0, 1) == SUCCESS);
    // create a cycle and test that it is detected
    assert(sg_change_vertex_outflow(14, &G, 7) == SUCCESS);  // should not detect a cycle yet
    assert(sg_flow_downstream_safe(&G, 9, 2) == CYCLE_WARNING);  // should detect a cycle now

    sg_outer_ocn_loop(10, &G);
    sg_display(G);
}