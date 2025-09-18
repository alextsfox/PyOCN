#include "streamgraph.c"
#include "timeit.c"
#include "returncodes.h"
#include <stdbool.h>
#include <stdint.h>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>


int main(){
    
    /*
    O O O O
    | | | |
    O O O O
    | | | |
    O O O O
    | | | |
    O-O-O-X
    */


    // drainedarea_t drained_area;
    // linidx_t adown;
    // localedges_t edges;
    // clockhand_t downstream;
    // uint8_t visited;

    const cartidx_t m = 4;
    const cartidx_t n = 4;

    Vertex cart_vertices[] = {
        (Vertex){1, 0,  0b00010000, 4, 0}, (Vertex){10, 0, 0b00010000, 4, 0}, (Vertex){10, 0, 0b00010000, 4, 0}, (Vertex){10, 0, 0b00010000, 4, 0},
        (Vertex){20, 0, 0b00010001, 4, 0}, (Vertex){20, 0, 0b00010001, 4, 0}, (Vertex){20, 0, 0b00010001, 4, 0}, (Vertex){20, 0, 0b00010001, 4, 0},
        (Vertex){30, 0, 0b00010001, 4, 0}, (Vertex){30, 0, 0b00010001, 4, 0}, (Vertex){30, 0, 0b00010001, 4, 0}, (Vertex){30, 0, 0b00010001, 4, 0},
        (Vertex){40, 0, 0b00000101, 2, 0}, (Vertex){80, 0, 0b01000101, 2, 0}, (Vertex){120, 0,0b01000101, 2, 0}, (Vertex){160, 0,0b01000001, ROOT_NODE, 0},
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
                    vert.adown = sg_cart_to_lin(row_offsets[vert.downstream], col_offsets[vert.downstream], m, n);
                    lin_vertices[a] = vert;
                    printf("i:%d, j:%d, a:%d, adown:%d, edges:%d\n", row, col, a, vert.adown, vert.edges);
                }
            }
        }
    }

    StreamGraph G;
    sg_create(&G, m, n, 3, 3);
    G.energy = 0.0;
    G.vertices = lin_vertices;

    sg_display(G);
    return SUCCESS;
}