#include <stdio.h>
#include <stdint.h>
#include <string.h>
#include <stdlib.h>
#include <stdbool.h>
#include <assert.h>
#include <locale.h>

#include "status.h"
#include "flowgrid.h"
#include "ocn.h"
#include "rng.h"

// Helper function to test for stream crossings
bool test_for_crossing(FlowGrid* fg. linidx_t a) {
    CartPair dims = fg->dims;
    Vertex vert = fg->vertices[a];
    clockhand_t down_new = vert.downstream;

    Vertex cross_check_vert;
    CartPair check_row_col = fg_lin_to_cart(a, dims);
    
    // Check appropriate vertex based on flow direction
    if (down_new == 1) {
        check_row_col.row -= 1;  // NE flow: check N vertex
    } else if (down_new == 7) {
        check_row_col.row -= 1;  // NW flow: check N vertex
    } else if (down_new == 3) {
        check_row_col.row += 1;  // SE flow: check S vertex
    } else if (down_new == 5) {
        check_row_col.row += 1;  // SW flow: check S vertex
    } else {
        return false; // Not a diagonal flow
    }

    fg_get_cart_safe(&cross_check_vert, fg. check_row_col);
    
    // Check for crossing streams
    if (down_new == 1 && (cross_check_vert.edges & (1u << 3))) return true;  // NE flow: N vertex has SE edge
    if (down_new == 7 && (cross_check_vert.edges & (1u << 5))) return true;  // NW flow: N vertex has SW edge
    if (down_new == 3 && (cross_check_vert.edges & (1u << 1))) return true;  // SE flow: S vertex has NE edge
    if (down_new == 5 && (cross_check_vert.edges & (1u << 7))) return true;  // SW flow: S vertex has NW edge
    
    return false;
}

int main(void) {
    setlocale(LC_ALL, "");
    printf("Running FlowGrid tests...\n\n");
    
    /********************************************
     * TEST 1: X-shaped Crossing Detection
     ********************************************/
    printf("TEST 1: X-shaped Crossing Detection\n");
    
    CartPair dims = {4, 4};
    FlowGrid *fg = fg_create_empty_safe(dims);
    if (fg == NULL) {
        printf("Failed to create FlowGrid.\n");
        return 1;
    }
    fg->root = (CartPair){0, 0};
    
    // Initialize with vertices that form X-shaped crossings
    Vertex vertices[16] = {
        {.drained_area = 1, .adown = 4, .edges = 16, .downstream = 4, .visited = 0},
        {.drained_area = 1, .adown = 5, .edges = 16, .downstream = 4, .visited = 0},
        {.drained_area = 1, .adown = 6, .edges = 16, .downstream = 4, .visited = 0},
        {.drained_area = 1, .adown = 7, .edges = 16, .downstream = 4, .visited = 0},
        {.drained_area = 2, .adown = 8, .edges = 17, .downstream = 4, .visited = 0},
        {.drained_area = 2, .adown = 10, .edges = 9, .downstream = 3, .visited = 0},  // Diagonal SE
        {.drained_area = 2, .adown = 9, .edges = 33, .downstream = 5, .visited = 0},  // Diagonal SW
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
    
    memcpy(fg->vertices, vertices, sizeof(vertices));
    
    printf("Graph with X-shaped crossing:\n");
    fg_display(fg. false);
    
    // Test crossing detection
    bool has_crossing_5 = test_for_crossing(fg. 5);
    bool has_crossing_6 = test_for_crossing(fg. 6);
    printf("Vertex 5 has crossing: %s\n", has_crossing_5 ? "YES" : "NO");
    printf("Vertex 6 has crossing: %s\n", has_crossing_6 ? "YES" : "NO");
    
    assert(has_crossing_5);
    assert(has_crossing_6);
    printf("X-crossing detection test passed!\n\n");
    
    /********************************************
     * TEST 2: Cycle Detection
     ********************************************/
    printf("TEST 2: Cycle Detection\n");
    
    // Create a cycle: 5 → 6 → 10 → 9 → 5
    Vertex vertices_cycle[16] = {
        {.drained_area = 1, .adown = 4, .edges = 16, .downstream = 4, .visited = 0},
        {.drained_area = 1, .adown = 5, .edges = 16, .downstream = 4, .visited = 0},
        {.drained_area = 1, .adown = 6, .edges = 16, .downstream = 4, .visited = 0},
        {.drained_area = 1, .adown = 7, .edges = 16, .downstream = 4, .visited = 0},
        {.drained_area = 2, .adown = 8, .edges = 17, .downstream = 4, .visited = 0},
        {.drained_area = 2, .adown = 6, .edges = 21, .downstream = 4, .visited = 0},  // 5 → 6
        {.drained_area = 2, .adown = 10, .edges = 81, .downstream = 4, .visited = 0}, // 6 → 10
        {.drained_area = 2, .adown = 11, .edges = 17, .downstream = 4, .visited = 0},
        {.drained_area = 3, .adown = 12, .edges = 17, .downstream = 4, .visited = 0},
        {.drained_area = 3, .adown = 5, .edges = 5, .downstream = 4, .visited = 0},   // 9 → 5
        {.drained_area = 3, .adown = 9, .edges = 65, .downstream = 4, .visited = 0},  // 10 → 9
        {.drained_area = 3, .adown = 15, .edges = 17, .downstream = 4, .visited = 0},
        {.drained_area = 4, .adown = 13, .edges = 5, .downstream = 2, .visited = 0},
        {.drained_area = 8, .adown = 14, .edges = 69, .downstream = 2, .visited = 0},
        {.drained_area = 12, .adown = 15, .edges = 69, .downstream = 2, .visited = 0},
        {.drained_area = 16, .adown = 0, .edges = 65, .downstream = 255, .visited = 0}
    };

    memcpy(fg->vertices, vertices_cycle, sizeof(vertices_cycle));

    printf("Graph with cycle 5 → 6 → 10 → 9 → 5:\n");
    
    // Test cycle detection starting from different vertices in the cycle
    linidx_t cycle_vertices[] = {5, 6, 9, 10};
    bool all_detected = true;
    
    for (int i = 0; i < sizeof(cycle_vertices)/sizeof(linidx_t); i++) {
        for (linidx_t j = 0; j < (dims.row * dims.col); j++) {
            fg->vertices[j].visited = 0;  // Reset visited flags
        }

        Status code = fg_flow_downstream_safe(fg. cycle_vertices[i], 1);
        printf("Checking for cycle starting at vertex %d: %s\n", 
               cycle_vertices[i], 
               (code == MALFORMED_GRAPH_WARNING) ? "CYCLE DETECTED" : "NO CYCLE");
        
        if (code != MALFORMED_GRAPH_WARNING) {
            all_detected = false;
        }
    }
    
    assert(all_detected);
    printf("Cycle detection test passed!\n");
    
    // Clean up
    fg_destroy_safe(fg);
    
    return 0;
}