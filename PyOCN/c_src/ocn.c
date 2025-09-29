#ifndef OCN_C
#define OCN_C

#include <stdint.h>
#include <stdlib.h>
#include <stdbool.h>
#include <math.h>
#include <string.h>
#include <stdio.h>
#include <wchar.h>
#include <locale.h>

#include "status.h"
#include "flowgrid.h"
#include "ocn.h"

/**
 * @brief Simulated annealing acceptance criterion.
 * @param energy_new The energy of the new state.
 * @param energy_old The energy of the old state.
 * @param temperature The current temperature.
 * @return true if the new state is accepted, false otherwise.
 */
bool simulate_annealing(double energy_new, double energy_old, double temperature){
    double u = (double)rand() / (double)RAND_MAX;
    double delta_energy = energy_new - energy_old;
    double p = exp(-delta_energy / temperature);
    return u < p;
}

/**
 * @brief Update the drained area along the downstream path from a given vertex. Unsafe.
 * @param G Pointer to the FlowGrid.
 * @param da_inc The increment to add to the drained area.
 * @param a The linear index of the starting vertex.
 * @return Status code indicating success or failure
 */
Status update_drained_area(FlowGrid *G, drainedarea_t da_inc, linidx_t a){
    Vertex vert;
    do {
        vert = fg_get_lin(G, a);
        vert.drained_area += da_inc;
        fg_set_lin(G, vert, a);
        a = vert.adown;
    } while (vert.downstream != IS_ROOT);

    return SUCCESS;  
}

double ocn_compute_energy(FlowGrid *G, double gamma){
    double energy = 0.0;
    for (linidx_t i = 0; i < (linidx_t)G->dims.row * (linidx_t)G->dims.col; i++){
        energy += pow(G->vertices[i].drained_area, gamma);
    }
    return energy;
}

/**
 * @brief Update the energy of the flowgrid along a single downstream path from a given vertex. Unsafe.
 * This function only works correctly if there is a single root in the flowgrid.
 * Having one root allows for this more efficient computation.
 * @param G Pointer to the FlowGrid. Modified in-place.
 * @param da_inc The increment to add to the drained area.
 * @param a The linear index of the starting vertex.
 * @param gamma The exponent used in the energy calculation.
 * @return Status code indicating success or failure
 */
Status update_energy_single_root(FlowGrid *G, drainedarea_t da_inc, linidx_t a, double gamma){
    Vertex vert;
    double energy_old = 0.0;
    double energy_new = 0.0;
    do {
        vert = fg_get_lin(G, a);

        energy_old += pow(vert.drained_area, gamma);
        vert.drained_area += da_inc;
        energy_new += pow(vert.drained_area, gamma);
        fg_set_lin(G, vert, a);
        a = vert.adown;
    } while (vert.downstream != IS_ROOT);
    G->energy += energy_new - energy_old;
    return SUCCESS;
}

Status ocn_single_erosion_event(FlowGrid *G, double gamma, double temperature){
    Status code;

    Vertex vert;
    clockhand_t down_old, down_new;
    linidx_t a, a_down_old, a_down_new;
    int32_t a_step_dir;
    drainedarea_t da_inc;
    CartPair dims = G->dims;
    linidx_t nverts = (linidx_t)dims.row * (linidx_t)dims.col;

    double energy_old, energy_new;
    energy_old = G->energy;
    
    a = rand() % nverts;  // pick a random vertex
    a_step_dir = (rand() % 2)*2 - 1;  // pick a random direction to step in if we need a new vertex
    down_new = rand() % 8;  // pick a random new downstream direction


    for (linidx_t nverts_tried = 0; nverts_tried < nverts; nverts_tried++){  // try a new vertex each time, up to the number of vertices in the graph
        // clunky way to wrap around, since apparently % on negative numbers is confusing as hell in C
        if (a == 0 && a_step_dir == -1) a = nverts - 1;
        else if (a == nverts - 1 && a_step_dir == 1) a = 0;
        else a = (linidx_t)((int32_t)a + a_step_dir);
        // a = (linidx_t)((int32_t)a + a_step_dir) % ((int32_t)dims.row * (int32_t)dims.col);

        vert = fg_get_lin(G, a);  // unsafe is ok here because a is guaranteed to be in bounds
    
        down_old = vert.downstream;
        a_down_old = vert.adown;
        da_inc = vert.drained_area;
        
        for (uint8_t ntries = 0; ntries < 8; ntries++){  // try a new direction each time, up to 8 times. Count these as separate tries.
            down_new  = (down_new + 1) % 8;

            code = fg_change_vertex_outflow(G, a, down_new);
            if (code != SUCCESS) continue;
            
            // retrieve the downstream vertices
            vert = fg_get_lin(G, a);
            a_down_new = vert.adown;

            // confirm that the new graph is well-formed (no cycles, still reaches root)
            for (linidx_t i = 0; i < nverts; i++) G->vertices[i].visited = 0;
            code = fg_flow_downstream_safe(G, a_down_old, 1);
            if (code != SUCCESS){
                fg_change_vertex_outflow(G, a, down_old);  // undo the swap, try again
                continue;
            }
            // for (linidx_t i = 0; i < nverts; i++) G->vertices[i].visited = 0;
            code = fg_flow_downstream_safe(G, a, 2);
            if (code != SUCCESS){
                fg_change_vertex_outflow(G, a, down_old);  // undo the swap, try again
                continue;
            }

            if (code == SUCCESS) goto mh_eval;  // if we reached here, the swap resulted in a well-formed graph, so we can move on the acceptance step
        }
    }
    return MALFORMED_GRAPH_WARNING; // we tried every vertex and every direction and couldn't find a valid swap.

    
    mh_eval:
    /*
    TODO:CRITICAL PERFORMANCE ISSUE:
    This function is supposed to update the energy of the flowgrid G after a 
    change in drained area along the path starting at vertex a. The previous method
    did not work: it would update the drained area along a path, then try to recompute
    the energy by adding pow((da + da_inc), gamma) - pow(da, gamma) for each vertex
    along the path. This is incorrect because the energy of each vertex depends on the
    drained area of all upstream vertices, not just its own drained area. This is fine 
    when gamma = 1, gamma = 0, or the number of roots = 1, but otherwise it is wrong.

    Simple but inefficient fix (current): recompute the *entire* energy of the flowgrid from scratch
    each time this function is called.

    More complex fix: find the set of all upstream vertices that flow into a and compute
    their summed contribution to the energy. Pass this value (sum of (da^gamma) for all
    upstream vertices) into this function, instead of just passing da_inc.
    */
    if ((G->nroots > 1) && (gamma < 1.0)){
        // energy_old = ocn_compute_energy(G, gamma);  // recompute energy from scratch
        update_drained_area(G, -da_inc, a_down_old);  // remove drainage from old path
        update_drained_area(G, da_inc, a_down_new);  // add drainage to new path
        energy_new = ocn_compute_energy(G, gamma);  // recompute energy from scratch
        // simulated annealing: accept with prob = exp(-delta_energy / temperature). note that p > 1 if energy decreases.
        // printf("Old energy: %f, New energy: %f, Delta E: %f\n", energy_old, energy_new, energy_new - energy_old);
        if (simulate_annealing(energy_new, energy_old, temperature)){
            G->energy = energy_new;
            return SUCCESS;
        }
        // reject swap: undo everything and try again
        update_drained_area(G, da_inc, a_down_old);  // add removed drainage back to old path
        update_drained_area(G, -da_inc, a_down_new);  // remove added drainage from new path
        fg_change_vertex_outflow(G, a, down_old);  // undo the outflow change
    } else {  // if there's only one root, we can use a more efficient method
        update_energy_single_root(G, -da_inc, a_down_old, gamma);  // remove drainage from old path and update energy
        update_energy_single_root(G, da_inc, a_down_new, gamma);  // add drainage to new path and update energy
        energy_new = G->energy;
        if (simulate_annealing(energy_new, energy_old, temperature)){
            return SUCCESS;
        }
        // reject swap: undo everything and try again
        update_energy_single_root(G, da_inc, a_down_old, gamma);  // add removed drainage back to old path and update energy
        update_energy_single_root(G, -da_inc, a_down_new, gamma);  // remove added drainage from new path and update energy
        fg_change_vertex_outflow(G, a, down_old);  // undo the outflow change
    }
    
    
    return EROSION_FAILURE;  // if we reach here, we failed to find a valid swap in many, many tries
}

Status ocn_outer_ocn_loop(FlowGrid *G, uint32_t niterations, double gamma, double *annealing_schedule){
    Status code;
    for (uint32_t i = 0; i < niterations; i++){
        code = ocn_single_erosion_event(G, gamma, annealing_schedule[i]);
        if ((code != SUCCESS) && (code != EROSION_FAILURE)) return code;
    }
    return SUCCESS;
}

#endif // OCN_C
