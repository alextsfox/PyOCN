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

Status ocn_single_erosion_event(StreamGraph *G, uint32_t *total_tries, double gamma, double temperature){
    Status code;

    Vertex vert;//, vert_down_old, vert_down_new;  unused?
    clockhand_t down_old, down_new;
    linidx_t a, a_down_old, a_down_new;
    
    CartPair dims = G->dims;

    double energy_old, energy_new;
    energy_old = G->energy;
    bool accept_bad_value = ((double)rand() / (double)RAND_MAX) > (1 - temperature);  // Metro-Hastings criterion. High temperature = more likely to accept bad values
    
    a = rand() % (dims.row * dims.col);  // pick a random vertex
    
    for ((*total_tries) = 0; (*total_tries) < (dims.row) * (dims.col) * 8; (*total_tries) += 8){
        a = (a + 1) % (dims.row * dims.col);  // try each vertex in turn, wrapping around to 0 after reaching the end
        down_new = rand() % 8;  // pick a random new downstream direction


        vert = sg_get_lin(G, a);  // unsafe is ok here because a is guaranteed to be in bounds

        down_old = vert.downstream;
        a_down_old = vert.adown;
        drainedarea_t da_inc = vert.drained_area;
        
        for (uint8_t ntries = 0; ntries < 8; ntries++){  // try a new direction each time, up to 8 times before giving up and picking a new vertex

            down_new  = (down_new + 1) % 8;

            code = sg_change_vertex_outflow(G, a, down_new);
            if (code != SUCCESS) continue;
            
            // retrieve the downstream vertices
            vert = sg_get_lin(G, a);
            a_down_new = vert.adown;
            // unused?
            // vert_down_new = sg_get_lin(G, a_down_new);  // bounds check was performed by sg_change_vertex_outflow
            // vert_down_old = sg_get_lin(G, a_down_old);  // bounds check was performed by sg_change_vertex_outflow

            // TODO: do we actually need to check for cycles from a_down_old?
            // confirm that the new graph is well-formed (no cycles, still reaches root)
            // for (linidx_t i = 0; i < (dims.row * dims.col); i++) G->vertices[i].visited = 0;
            // code = sg_flow_downstream_safe(G, a_down_old, 1);
            // if (code != SUCCESS){
            //     sg_change_vertex_outflow(G, a, down_old);  // undo the swap, try again
            //     continue;
            // }
            // for (linidx_t i = 0; i < (dims.row * dims.col); i++) G->vertices[i].visited = 0;
            // code = sg_flow_downstream_safe(G, a, 1);  // use ncalls = 2. This way, if the two paths intersect, it won't trigger a cycle warning.
            // if (code != SUCCESS){
            //     sg_change_vertex_outflow(G, a, down_old);  // undo the swap, try again
            //     continue;
            // }

            // more comprehensive cycle check: check *all* vertices
            for (linidx_t i = 0; i < (dims.row * dims.col); i++){
                for (linidx_t j = 0; j < (dims.row * dims.col); j++) G->vertices[j].visited = 0;
                code = sg_flow_downstream_safe(G, i, 1);
                if (code != SUCCESS) break;  // found a cycle, break out to try a new direction
            }
            if (code != SUCCESS){
                sg_change_vertex_outflow(G, a, down_old);  // undo the swap, try again
                continue;
            }

            break;  // if we reached here, the swap resulted in a well-formed graph, so we can move on the acceptance step
        }

        if (code != SUCCESS) continue;  // if we exhausted all 8 directions without a successful swap, try again. do not go to acceptance step

        // update drained area and energy along both paths
        ocn_update_energy(G, -da_inc, a_down_old, gamma);  // decrement drained area along old path
        ocn_update_energy(G, da_inc, a_down_new, gamma);  // increment drained area along new path
        energy_new = G->energy;

        if (energy_new < energy_old || accept_bad_value) return SUCCESS;

        // reject swap: undo everything and try again
        ocn_update_energy(G, da_inc, a_down_old, gamma);  // undo the decrement
        ocn_update_energy(G, -da_inc, a_down_new, gamma);  // undo the increment
        sg_change_vertex_outflow(G, a, a_down_old);  // undo the outflow change
    }
    return EROSION_FAILURE;  // if we reach here, we failed to find a valid swap in many, many tries
}

Status ocn_outer_ocn_loop(StreamGraph *G, uint32_t niterations, double gamma, double *annealing_schedule){
    Status code;
    uint32_t i = 0;
    uint32_t total_tries = 0;
    while (i < niterations){
        code = ocn_single_erosion_event(G, &total_tries, gamma, annealing_schedule[i]);
        if ((code != SUCCESS) && (code != EROSION_FAILURE)) return code;  // malformed graph?
        i += total_tries;
    }
    return SUCCESS;
}

#endif // OCN_C