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
    while (vert.downstream != IS_ROOT){
        // update drained area of the vertex and energy of the graph.
        da = vert.drained_area;
        G->energy += pow((double)(da + da_inc), gamma) - pow((double)da, gamma);  // update energy
        vert.drained_area += da_inc;  // update drained area of vertex
        sg_set_lin(G, vert, a);

        // get next vertex
        a = vert.adown;
        vert = sg_get_lin(G, a);
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

    Vertex vert, vert_down_old, vert_down_new;
    clockhand_t down_old, down_new;
    linidx_t a, a_down_old, a_down_new;
    
    CartPair dims = G->dims;

    drainedarea_t da_inc;
    double energy_old, energy_new;
    energy_old = G->energy;
    bool accept_bad_value = ((double)rand() / (double)RAND_MAX) > (1 - temperature);  // Metro-Hastings criterion. High temperature = more likely to accept bad values

    linidx_t i;
    uint8_t ntries;
    for (uint16_t total_tries = 0; total_tries < 1000; total_tries++){
        // pick a random vertex and a random new downstream direction,
        down_new = rand() % 8;
        a = rand() % (dims.row * dims.col);
        vert = sg_get_lin(G, a);  // unsafe is ok here because a is guaranteed to be in bounds
        
        down_old = vert.downstream;
        a_down_old = vert.adown;
        da_inc = vert.drained_area;

        for (ntries = 0; ntries < 8; ntries++){  // try a new direction each time, up to 8 times before giving up and picking a new vertex
            down_new  = (down_new + 1) % 8;

            code = sg_change_vertex_outflow(G, a, down_new);
            if (code != SUCCESS) continue;
            
            // retrieve the downstream vertices
            vert = sg_get_lin(G, a);
            a_down_new = vert.adown;
            vert_down_new = sg_get_lin(G, a_down_new);  // bounds check was performed by sg_change_vertex_outflow
            vert_down_old = sg_get_lin(G, a_down_old);  // bounds check was performed by sg_change_vertex_outflow

            // zero the cycle tracker before updating drained area
            for (i = 0; i < (dims.row * dims.col); i++) G->vertices[i].visited = 0;

            // confirm that both paths lead to the root without cycles
            sg_flow_downstream_safe(G, a_down_old, 1);
            if (code != SUCCESS){
                sg_change_vertex_outflow(G, a, down_old);  // undo the swap, try again
                continue;
            }
            code = sg_flow_downstream_safe(G, a_down_new, 2);  // use ncalls = 2. This way, if the two paths intersect, it won't trigger a cycle warning.
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

    return MALFORMED_GRAPH_WARNING;  // if we reach here, we failed to find a valid swap in many, many tries
}

Status ocn_outer_ocn_loop(StreamGraph *G, uint32_t niterations, double gamma, double *annealing_schedule){
    Status code;
    uint32_t i = 0;
    for (i = 0; i < niterations; i++){
        code = ocn_single_erosion_event(G, gamma, annealing_schedule[i]);
        if (code != SUCCESS) return code;  // malformed graph?
    }
    return SUCCESS;
}

#endif // OCN_C