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
#include "rng.h"

/**
 * @brief Simulated annealing acceptance criterion.
 * @param energy_new The energy of the new state.
 * @param energy_old The energy of the old state.
 * @param temperature The current temperature.
 * @return true if the new state is accepted, false otherwise.
 */
static inline bool simulate_annealing(double energy_new, double energy_old, double temperature, rng_state_t *rng){
    const double accpt_prob = rng_uniformdouble(rng);
    const double delta_energy = energy_new - energy_old;
    if (delta_energy <= 0.0) return true;  // Always accept improvements
    const double p = exp(-delta_energy / temperature);
    return accpt_prob < p;
}

/**
 * @brief Update the drained area along the downstream path from a given vertex. Unsafe.
 * @param G Pointer to the FlowGrid.
 * @param da_inc The increment to add to the drained area.
 * @param a The linear index of the starting vertex.
 * @return Status code indicating success or failure
 */
static Status update_drained_area(FlowGrid *G, drainedarea_t da_inc, linidx_t a){
    Vertex vert;
    do {
        Status code = fg_get_lin(&vert, G, a);
        if (code == OOB_ERROR) return OOB_ERROR;
        vert.drained_area += da_inc;
        code = fg_set_lin(G, vert, a);
        if (code == OOB_ERROR) return OOB_ERROR;
        a = vert.adown;
    } while (vert.downstream != IS_ROOT);

    return SUCCESS;  
}

static inline double compute_energy(FlowGrid *G, double gamma){
    double energy = 0.0;
    for (linidx_t i = 0; i < (linidx_t)G->dims.row * (linidx_t)G->dims.col; i++){
        energy += pow(G->vertices[i].drained_area, gamma);
    }
    return energy;
}

double ocn_compute_energy(FlowGrid *G, double gamma){
    return compute_energy(G, gamma);
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
static Status update_energy_single_root(FlowGrid *G, drainedarea_t da_inc, linidx_t a, double gamma){
    Vertex vert;
    double energy_old = 0.0;
    double energy_new = 0.0;
    do {
        if (fg_get_lin(&vert, G, a) == OOB_ERROR) return OOB_ERROR;

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
    int8_t down_step_dir;
    linidx_t a, a_down_old, a_down_new;
    int32_t a_step_dir;
    drainedarea_t da_inc;
    CartPair dims = G->dims;
    linidx_t nverts = (linidx_t)dims.row * (linidx_t)dims.col;

    double energy_old, energy_new;
    energy_old = G->energy;
    
    a = rng_randint32(&G->rng) % nverts;  // pick a random vertex
    a_step_dir = (rng_randint32(&G->rng) % 2)*2 - 1;  // pick a random direction to step in if we need a new vertex
    down_new = rng_randint32(&G->rng) >> 29;  // bit shift of 29 gives a number from 0-7
    down_step_dir = (rng_randint32(&G->rng) % 2)*2 - 1;  // pick a random direction to step in if we need a new direction

    for (linidx_t nverts_tried = 0; nverts_tried < nverts; nverts_tried++){  // try a new vertex each time, up to the number of vertices in the graph
        // clunky way to wrap around, since apparently % on negative numbers is confusing as hell in C
        if (a == 0 && a_step_dir == -1) a = nverts - 1;
        else if (a == nverts - 1 && a_step_dir == 1) a = 0;
        else a = (linidx_t)((int32_t)a + a_step_dir);
        if (down_new == 0 && down_step_dir == -1) down_new = 7;
        else if (down_new == 7 && down_step_dir == 1) down_new = 0;
        else down_new = (clockhand_t)((int8_t)down_new + down_step_dir);

        code = fg_get_lin(&vert, G, a);
        if (code == OOB_ERROR) return OOB_ERROR;
    
        down_old = vert.downstream;
        a_down_old = vert.adown;
        da_inc = vert.drained_area;
        
        for (uint8_t ndirections_tried = 0; ndirections_tried < 8; ndirections_tried++){  // try a new direction each time.
            down_new  = (down_new + 1) % 8;

            // propose to rerout the outflow of vertex a to direction down_new
            code = fg_change_vertex_outflow(G, a, down_new);
            if (code != SUCCESS) continue;

            // retrieve the new downstream vertex index
            code = fg_get_lin(&vert, G, a);
            if (code == OOB_ERROR) return OOB_ERROR;
            a_down_new = vert.adown;

            // confirm that the new graph is well-formed (no cycles, still reaches root)
            for (linidx_t i = 0; i < nverts; i++) G->vertices[i].visited = 0;  // reset visitation flags
            code = fg_check_for_cycles(G, a_down_old, 1);
            if (code != SUCCESS){
                fg_change_vertex_outflow(G, a, down_old);  // undo the swap, try again
                continue;
            }
            code = fg_check_for_cycles(G, a, 2);  // no need to reset visitation flags, instead use a different check number
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
    TODO: PERFORMANCE ISSUE:
    This function is supposed to update the energy of the flowgrid G after a 
    change in drained area along the path starting at vertex a.

    Simple but inefficient fix (current): recompute the *entire* energy of the flowgrid from scratch
    each time this function is called.

    More complex fix: find the set of all upstream vertices that flow into a and compute
    their summed contribution to the energy. Pass this value (sum of (da^gamma) for all
    upstream vertices) into this function, instead of just passing da_inc.
    */
    if ((G->nroots > 1) && (gamma < 1.0)){
        // energy_old = compute_energy(G, gamma);  // recompute energy from scratch
        update_drained_area(G, -da_inc, a_down_old);  // remove drainage from old path
        update_drained_area(G, da_inc, a_down_new);  // add drainage to new path
        energy_new = compute_energy(G, gamma);  // recompute energy from scratch
        if (simulate_annealing(energy_new, energy_old, temperature, &G->rng)){
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
        if (simulate_annealing(energy_new, energy_old, temperature, &G->rng)){
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
