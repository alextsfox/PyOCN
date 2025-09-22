#include <stdio.h>
#include "streamgraph.c"
#include "returncodes.h"
#include <stdint.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <stdbool.h>
#include <assert.h>
#include <locale.h>

#define USE_UTF8 true

int main(void){
    setlocale(LC_ALL, "");
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

    sg_single_erosion_event(&G);
}