#ifndef OCN_H
#define OCN_H

#include "status.h"
#include "streamgraph.h"

/**
 * @brief Update the energy of the streamgraph after a change in drained area at a specific vertex. Unsafe! Assumed a well-formed graph and valid index/pointer.
 * @param G Pointer to the StreamGraph.
 * @param da_inc The change in drained area (can be positive or negative).
 * @param a The linear index of the vertex where the drained area change occurs.
 * @param gamma The exponent used in the energy calculation.
 * @return Status code indicating success or failure
 */
Status ocn_update_energy(StreamGraph *G, drainedarea_t da_inc, linidx_t a, double gamma);

/**
 * @brief Perform a single erosion event on the streamgraph, attempting to change the outflow of a random vertex.
 * Selects a random vertex and a random new direction and attempts to modify the graph accordingly. 
 * If the modification results in a malformed graph, it is undone and the process is retried up to a set number of times.
 * In case of a failure, each vertex is tried 8 different times, once for each possible new direction.
 * If no valid modification is found for a given vertex, a new vertex is selected and the process repeats.
 * If no valid modification is found after trying a set number of vertices, the function exits with a warning.
 * Unsafe! Assumed a well-formed graph.
 * @param G Pointer to the StreamGraph.
 * @param gamma The exponent used in the energy calculation.
 * @param temperature The temperature parameter for the Metropolis-Hastings acceptance criterion.
 * @return Status code indicating success or failure.
 */
Status ocn_single_erosion_event(StreamGraph *G, double gamma, double temperature);

/**
 * @brief Perform multiple erosion events on the streamgraph.
 * @param G Pointer to the StreamGraph.
 * @param niterations The number of erosion events to perform.
 * @param gamma The exponent used in the energy calculation.
 * @param annealing_schedule An array of temperatures (ranging from 0-1) to use for each iteration. Length must be at least niterations.
 * @return Status code indicating success or failure
 */
Status ocn_outer_ocn_loop(StreamGraph *G, uint32_t niterations, double gamma, double *annealing_schedule);

#endif // OCN_H