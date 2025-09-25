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

    //////////////////////////////////////////////////////////////////////////////////////
    // X detection test
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
        {.drained_area = 2, .adown = 8, .edges = 17, .downstream = 4, .visited = 0},
        {.drained_area = 2, .adown = 10, .edges = 9, .downstream = 3, .visited = 0},
        {.drained_area = 2, .adown = 9, .edges = 33, .downstream = 5, .visited = 0},
        {.drained_area = 2, .adown = 11, .edges = 17, .downstream = 4, .visited = 0},
        {.drained_area = 3, .adown = 12, .edges = 17, .downstream = 4, .visited = 0},
        {.drained_area = 3, .adown = 13, .edges = 18, .downstream = 4, .visited = 0},
        {.drained_area = 3, .adown = 14, .edges = 144, .downstream = 4, .visited = 0},
        {.drained_area = 3, .adown = 15, .edges = 17, .downstream = 4, .visited = 0},
        {.drained_area = 4, .adown = 13, .edges = 5, .downstream = 2, .visited = 0},
        {.drained_area = 8, .adown = 14, .edges = 69, .downstream = 2, .visited = 0},
        {.drained_area = 12, .adown = 15, .edges = 69, .downstream = 2, .visited = 0},
        {.drained_area = 16, .adown = 0, .edges = 65, .downstream = 255, .visited = 0}
    };
    memcpy(sg.vertices, vertices, sizeof(vertices));
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

    bool triggered = false;
    if (down_new == 1) {
        if (cross_check_vert.edges & (1u << 3)) triggered = true;  // NE flow: N vertex cannot have a SE edge
    } else if (down_new == 7) {
        if (cross_check_vert.edges & (1u << 5)) triggered = true;  // NW flow: N vertex cannot have a SW edge
    } else if (down_new == 3) {
        if (cross_check_vert.edges & (1u << 1)) triggered = true;  // SE flow: S vertex cannot have a NE edge
    } else if (down_new == 5) {
        if (cross_check_vert.edges & (1u << 7)) triggered = true;  // SW flow: S vertex cannot have a NW edge
    }
    assert(triggered);
    
    a = 6;
    vert = sg.vertices[a];
    down_new = vert.downstream;

    check_row_col = sg_lin_to_cart(a, dims);
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

    triggered = false;
    if (down_new == 1) {
        if (cross_check_vert.edges & (1u << 3)) triggered = true;  // NE flow: N vertex cannot have a SE edge
    } else if (down_new == 7) {
        if (cross_check_vert.edges & (1u << 5)) triggered = true;  // NW flow: N vertex cannot have a SW edge
    } else if (down_new == 3) {
        if (cross_check_vert.edges & (1u << 1)) triggered = true;  // SE flow: S vertex cannot have a NE edge
    } else if (down_new == 5) {
        if (cross_check_vert.edges & (1u << 7)) triggered = true;  // SW flow: S vertex cannot have a NW edge
    }
    assert(triggered);

    ////////////////////////////////////////////////////////////////////////////////////////
    // Cycle dectection test
    // contains a cycle
    Vertex vertices_new[16] = {
        {.drained_area = 1, .adown = 4, .edges = 16, .downstream = 4, .visited = 0},
        {.drained_area = 1, .adown = 5, .edges = 16, .downstream = 4, .visited = 0},
        {.drained_area = 1, .adown = 6, .edges = 16, .downstream = 4, .visited = 0},
        {.drained_area = 1, .adown = 7, .edges = 16, .downstream = 4, .visited = 0},
        {.drained_area = 2, .adown = 8, .edges = 17, .downstream = 4, .visited = 0},
        {.drained_area = 2, .adown = 6, .edges = 21, .downstream = 4, .visited = 0},  // 5 goes to 6
        {.drained_area = 2, .adown = 10, .edges = 81, .downstream = 4, .visited = 0},  // 6 goes to 10
        {.drained_area = 2, .adown = 11, .edges = 17, .downstream = 4, .visited = 0},
        {.drained_area = 3, .adown = 12, .edges = 17, .downstream = 4, .visited = 0},
        {.drained_area = 3, .adown = 5, .edges = 5, .downstream = 4, .visited = 0},  // 9 goes to 5
        {.drained_area = 3, .adown = 9, .edges = 65, .downstream = 4, .visited = 0},  // 10 goes to 9
        {.drained_area = 3, .adown = 15, .edges = 17, .downstream = 4, .visited = 0},
        {.drained_area = 4, .adown = 13, .edges = 5, .downstream = 2, .visited = 0},
        {.drained_area = 8, .adown = 14, .edges = 69, .downstream = 2, .visited = 0},
        {.drained_area = 12, .adown = 15, .edges = 69, .downstream = 2, .visited = 0},
        {.drained_area = 16, .adown = 0, .edges = 65, .downstream = 255, .visited = 0}
    };
    memcpy(sg.vertices, vertices_new, sizeof(vertices_new));
    sg.energy = 64.0;

    linidx_t a_to_test[] = {1, 2, 5, 6, 9, 10};
    for (int i = 0; i < sizeof(a_to_test) / sizeof(linidx_t); i++){
        for (linidx_t j = 0; j < (dims.row * dims.col); j++) sg.vertices[j].visited = 0;
        Status code = sg_flow_downstream_safe(&sg, a_to_test[i], 1);
        assert(code == MALFORMED_GRAPH_WARNING);
    }

    sg_destroy_safe(&sg);

    return 0;
}