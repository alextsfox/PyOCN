#include <stdio.h>
#include <stdint.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <stdbool.h>
#include <assert.h>
#include <locale.h>

#include "status.h"
#include "streamgraph.h"
#include "ocn.h"
#include "rng.h"

int main(void){
    setlocale(LC_ALL, "");
/*
    O O O O
    | | | |
    O O O O
    | | | |
    O O O O
    | | | |
    O-O-O-X but with 36 vertices
    */

    // const cartidx_t m  =  4;
    // const cartidx_t n = 4;

    // Vertex cart_vertices[] = {
    //     (Vertex){1, 0,  0b00010000, 4, 0}, (Vertex){1, 0, 0b00010000, 4, 0}, (Vertex){1, 0, 0b00010000, 4, 0}, (Vertex){1, 0, 0b00010000, 4, 0},
    //     (Vertex){2, 0, 0b00010001, 4, 0}, (Vertex){2, 0, 0b00010001, 4, 0}, (Vertex){2, 0, 0b00010001, 4, 0}, (Vertex){2, 0, 0b00010001, 4, 0},
    //     (Vertex){3, 0, 0b00010001, 4, 0}, (Vertex){3, 0, 0b00010001, 4, 0}, (Vertex){3, 0, 0b00010001, 4, 0}, (Vertex){3, 0, 0b00010001, 4, 0},
    //     (Vertex){4, 0, 0b00000101, 2, 0}, (Vertex){8, 0, 0b01000101, 2, 0}, (Vertex){12, 0,0b01000101, 2, 0}, (Vertex){16, 0,0b01000001, IS_ROOT, 0},
    // };
    // Vertex lin_vertices[m*n];
    // Vertex vert;
    // // awful hack until I write a proper conversion function i,j -> a
    // // TODO write proper conversion function
    // for (linidx_t a = 0; a < m*n; a++){
    //     for (cartidx_t row = 0; row < m; row++){
    //         for (cartidx_t col = 0; col < n; col++){
    //             if (a == sg_cart_to_lin(row, col, m, n)){
    //                 vert = cart_vertices[row*n + col];
    //                 vert.adown = sg_cart_to_lin(row + row_offsets[vert.downstream], col + col_offsets[vert.downstream], m, n);
    //                 lin_vertices[a] = vert;
    //                 // printf("i:%d, j:%d, a:%d, adown:%d, edges:%d, downstream:%d\n", row, col, a, vert.adown, vert.edges, vert.downstream);
    //             }
    //         }
    //     }
    // }


    // StreamGraph G;
    // sg_create(&G, m, n, 3, 3);
    // G.energy = 0.0;
    // G.vertices = lin_vertices;
    // G.energy = 4 + sqrt(2)*4 + sqrt(3)*4 + sqrt(4)*4;  // pre-calculate initial energy

    StreamGraph sg;
    CartPair dims = {6, 6};
    CartPair root = {5, 5};
    sg_create_empty_safe(&sg, root, dims);
    // StreamGraph(dims=(6, 6), root=(5, 5), energy=76.82643119030881, vertices=<36 vertices>)

    Vertex vertices[36] = {
        {.drained_area = 1, .adown = 1, .edges = 4, .downstream = 2, .visited = 0},
        {.drained_area = 1, .adown = 7, .edges = 80, .downstream = 4, .visited = 0},
        {.drained_area = 1, .adown = 8, .edges = 16, .downstream = 4, .visited = 0},
        {.drained_area = 1, .adown = 9, .edges = 16, .downstream = 4, .visited = 0},
        {.drained_area = 1, .adown = 10, .edges = 16, .downstream = 4, .visited = 0},
        {.drained_area = 1, .adown = 10, .edges = 32, .downstream = 5, .visited = 0},
        {.drained_area = 2, .adown = 12, .edges = 16, .downstream = 4, .visited = 0},
        {.drained_area = 2, .adown = 13, .edges = 16, .downstream = 4, .visited = 0},
        {.drained_area = 2, .adown = 15, .edges = 8, .downstream = 3, .visited = 0},
        {.drained_area = 2, .adown = 15, .edges = 16, .downstream = 4, .visited = 0},
        {.drained_area = 2, .adown = 16, .edges = 16, .downstream = 4, .visited = 0},
        {.drained_area = 2, .adown = 17, .edges = 16, .downstream = 4, .visited = 0},
        {.drained_area = 3, .adown = 18, .edges = 16, .downstream = 4, .visited = 0},
        {.drained_area = 3, .adown = 19, .edges = 16, .downstream = 4, .visited = 0},
        {.drained_area = 3, .adown = 20, .edges = 16, .downstream = 4, .visited = 0},
        {.drained_area = 3, .adown = 21, .edges = 144, .downstream = 4, .visited = 0},
        {.drained_area = 3, .adown = 22, .edges = 16, .downstream = 4, .visited = 0},
        {.drained_area = 3, .adown = 23, .edges = 16, .downstream = 4, .visited = 0},
        {.drained_area = 4, .adown = 24, .edges = 16, .downstream = 4, .visited = 0},
        {.drained_area = 4, .adown = 25, .edges = 16, .downstream = 4, .visited = 0},
        {.drained_area = 4, .adown = 26, .edges = 16, .downstream = 4, .visited = 0},
        {.drained_area = 4, .adown = 27, .edges = 16, .downstream = 4, .visited = 0},
        {.drained_area = 4, .adown = 28, .edges = 16, .downstream = 4, .visited = 0},
        {.drained_area = 4, .adown = 29, .edges = 16, .downstream = 4, .visited = 0},
        {.drained_area = 5, .adown = 30, .edges = 16, .downstream = 4, .visited = 0},
        {.drained_area = 5, .adown = 31, .edges = 16, .downstream = 4, .visited = 0},
        {.drained_area = 5, .adown = 31, .edges = 32, .downstream = 5, .visited = 0},
        {.drained_area = 5, .adown = 33, .edges = 16, .downstream = 4, .visited = 0},
        {.drained_area = 5, .adown = 34, .edges = 16, .downstream = 4, .visited = 0},
        {.drained_area = 5, .adown = 35, .edges = 16, .downstream = 4, .visited = 0},
        {.drained_area = 6, .adown = 31, .edges = 4, .downstream = 2, .visited = 0},
        {.drained_area = 12, .adown = 32, .edges = 4, .downstream = 2, .visited = 0},
        {.drained_area = 18, .adown = 33, .edges = 4, .downstream = 2, .visited = 3},
        {.drained_area = 24, .adown = 34, .edges = 4, .downstream = 2, .visited = 3},
        {.drained_area = 30, .adown = 35, .edges = 4, .downstream = 2, .visited = 3},
        {.drained_area = 36, .adown = 0, .edges = 0, .downstream = 255, .visited = 0}
    };
    
    sg.vertices = vertices;
    sg.energy = 76.82643119030881;

    for (int i = 0; i < 100; i++) {
        printf("started iteration %d\n", i);
        ocn_single_erosion_event(&sg, 0.5, 0.01);
    }


    sg_display(&sg, true);
}