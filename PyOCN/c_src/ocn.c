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
#include "streamgraph.h"
#include "ocn.h"

Status ocn_update_energy(StreamGraph *G, drainedarea_t da_inc, linidx_t a, double gamma){
    Vertex vert;
    drainedarea_t da;
    vert = sg_get_lin(G, a);
    int sanity_counter = 0;
    while (vert.downstream != IS_ROOT){
        // update drained area of the vertex and energy of the graph.
        da = vert.drained_area;
        G->energy += pow((double)(da + da_inc), gamma) - pow((double)da, gamma);  // update energy
        vert.drained_area += da_inc;  // update drained area of vertex
        sg_set_lin(G, vert, a);

        // get next vertex
        a = vert.adown;
        vert = sg_get_lin(G, a);

        sanity_counter++;
        if (sanity_counter > (G->dims.row * G->dims.col)){
            printf("Sanity check failed in energy update loop.\n");
        }
    }
    // update root vertex
    da = vert.drained_area;
    G->energy += pow((double)(da + da_inc), gamma) - pow((double)da, gamma);  // update energy
    vert.drained_area += da_inc;  // update drained area of vertex
    sg_set_lin(G, vert, a);

    return SUCCESS;
}

Status ocn_single_erosion_event(StreamGraph *G, double gamma, double temperature){
    Status code;

    Vertex vert;//, vert_down_old, vert_down_new;  unused?
    clockhand_t down_old, down_new;
    linidx_t a, a_down_old, a_down_new;
    drainedarea_t da_inc;
    CartPair dims = G->dims;

    double energy_old, energy_new;
    energy_old = G->energy;
    
    a = rand() % (dims.row * dims.col);  // pick a random vertex
    down_new = rand() % 8;  // pick a random new downstream direction

    
    for (linidx_t nverts_tried = 0; nverts_tried < (dims.row * dims.col); nverts_tried++){  // try a new vertex each time, up to the number of vertices in the graph
        a = (a + 1) % (dims.row * dims.col);
        vert = sg_get_lin(G, a);  // unsafe is ok here because a is guaranteed to be in bounds
    
        down_old = vert.downstream;
        a_down_old = vert.adown;
        da_inc = vert.drained_area;
        
        for (uint8_t ntries = 0; ntries < 8; ntries++){  // try a new direction each time, up to 8 times. Count these as separate tries.
            down_new  = (down_new + 1) % 8;

            code = sg_change_vertex_outflow(G, a, down_new);
            if (code != SUCCESS) continue;
            
            // retrieve the downstream vertices
            vert = sg_get_lin(G, a);
            a_down_new = vert.adown;

            // confirm that the new graph is well-formed (no cycles, still reaches root)
            for (linidx_t i = 0; i < (dims.row * dims.col); i++) G->vertices[i].visited = 0;
            code = sg_flow_downstream_safe(G, a_down_old, 1);
            if (code != SUCCESS){
                sg_change_vertex_outflow(G, a, down_old);  // undo the swap, try again
                continue;
            }
            // for (linidx_t i = 0; i < (dims.row * dims.col); i++) G->vertices[i].visited = 0;
            code = sg_flow_downstream_safe(G, a, 2);  // can use n_calls = 1 again because we reset visited flags
            if (code != SUCCESS){
                sg_change_vertex_outflow(G, a, down_old);  // undo the swap, try again
                continue;
            }

            if (code == SUCCESS) goto mh_eval;  // if we reached here, the swap resulted in a well-formed graph, so we can move on the acceptance step
        }
    }
    return MALFORMED_GRAPH_WARNING; // we tried every vertex and every direction and couldn't find a valid swap.

    // update drained area and energy along both paths
    mh_eval:
    ocn_update_energy(G, -da_inc, a_down_old, gamma);  // decrement drained area along old path
    ocn_update_energy(G, da_inc, a_down_new, gamma);  // increment drained area along new path
    energy_new = G->energy;

    // simulated annealing: accept with prob = exp(-delta_energy / temperature). note that p > 1 if energy decreases.
    double u = (double)rand() / (double)RAND_MAX;
    double delta_energy = energy_new - energy_old;
    double p = exp(-delta_energy / temperature);
    if (u < p) return SUCCESS;

    // reject swap: undo everything and try again
    ocn_update_energy(G, da_inc, a_down_old, gamma);  // undo the decrement
    ocn_update_energy(G, -da_inc, a_down_new, gamma);  // undo the increment
    sg_change_vertex_outflow(G, a, down_old);  // undo the outflow change
    
    return EROSION_FAILURE;  // if we reach here, we failed to find a valid swap in many, many tries
}

Status ocn_outer_ocn_loop(StreamGraph *G, double *energy_history, uint32_t niterations, double gamma, double *annealing_schedule){
    Status code;

    if (energy_history == NULL) return NULL_POINTER_ERROR;

    for (uint32_t i = 0; i < niterations; i++){
        code = ocn_single_erosion_event(G, gamma, annealing_schedule[i]);
        if ((code != SUCCESS) && (code != EROSION_FAILURE)) return code;
        energy_history[i] = G->energy;
    }
    return SUCCESS;
}

#endif // OCN_C