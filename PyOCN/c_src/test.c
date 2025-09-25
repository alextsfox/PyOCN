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

    // X shape, 5->10 and 6->9
    StreamGraph sg;
    CartPair dims = {4, 4};
    CartPair root = {3, 3};
    sg_create_empty_safe(&sg, root, dims);
    Vertex vertices[16] = {
        {.drained_area = 1, .adown = 4, .edges = 16, .downstream = 4, .visited = 0},
        {.drained_area = 1, .adown = 5, .edges = 16, .downstream = 4, .visited = 0},
        {.drained_area = 1, .adown = 6, .edges = 16, .downstream = 4, .visited = 0},
        {.drained_area = 1, .adown = 7, .edges = 16, .downstream = 4, .visited = 0},
        {.drained_area = 2, .adown = 8, .edges = 16, .downstream = 4, .visited = 0},
        {.drained_area = 2, .adown = 10, .edges = 8, .downstream = 3, .visited = 0},
        {.drained_area = 2, .adown = 9, .edges = 32, .downstream = 5, .visited = 0},
        {.drained_area = 2, .adown = 11, .edges = 16, .downstream = 4, .visited = 0},
        {.drained_area = 3, .adown = 12, .edges = 16, .downstream = 4, .visited = 0},
        {.drained_area = 3, .adown = 13, .edges = 16, .downstream = 4, .visited = 0},
        {.drained_area = 3, .adown = 14, .edges = 16, .downstream = 4, .visited = 0},
        {.drained_area = 3, .adown = 15, .edges = 16, .downstream = 4, .visited = 0},
        {.drained_area = 4, .adown = 13, .edges = 4, .downstream = 2, .visited = 0},
        {.drained_area = 8, .adown = 14, .edges = 4, .downstream = 2, .visited = 0},
        {.drained_area = 12, .adown = 15, .edges = 4, .downstream = 2, .visited = 0},
        {.drained_area = 16, .adown = 0, .edges = 0, .downstream = 255, .visited = 0}
    };
    sg.vertices = vertices;
    sg.energy = 64.0 ;

    rng_seed(8472);

    sg_display(&sg, true);


    // check that we haven't created any crosses
    linidx_t a = 5;
    Vertex vert = sg.vertices[a];
    clockhand_t down_new = vert.downstream;

    Vertex cross_check_vert;
    CartPair check_row_col = sg_lin_to_cart(a, dims);
    // Convert switch statements to if/else if/else structure
    if (down_new == 1) {
        check_row_col.row -= 1;  // NE flow: check N vertex
    } else if (down_new == 7) {
        check_row_col.row -= 1;  // NW flow: check N vertex
    } else if (down_new == 3) {
        check_row_col.row += 1;  // SE flow: check S vertex
    } else if (down_new == 5) {
        check_row_col.row += 1;  // SW flow: check S vertex
    }

    sg_get_cart_safe(&cross_check_vert, &sg, check_row_col);

    printf("CHECK: row %d col %d lin %d\n", check_row_col.row, check_row_col.col, sg_cart_to_lin(check_row_col, dims));
    printf("CROSS CHECK VERTEX EDGES: %u\n", cross_check_vert.edges);

    if (down_new == 1) {
        if (cross_check_vert.edges & (1u << 3)) printf("Trigger 1\n");  // NE flow: N vertex cannot have a SE edge
        else printf("No trigger\n");
    } else if (down_new == 7) {
        if (cross_check_vert.edges & (1u << 5)) printf("Trigger 2\n");  // NW flow: N vertex cannot have a SW edge
        else printf("No trigger\n");
    } else if (down_new == 3) {
        if (cross_check_vert.edges & (1u << 1)) printf("Trigger 3\n");  // SE flow: S vertex cannot have a NE edge
        else printf("No trigger\n");
    } else if (down_new == 5) {
        if (cross_check_vert.edges & (1u << 7)) printf("Trigger 4\n");  // SW flow: S vertex cannot have a NW edge
        else printf("No trigger\n");
    } else {
        printf("No trigger\n");
    }

    printf("DONE\n");

    return 0;
}