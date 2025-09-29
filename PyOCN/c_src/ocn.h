/**
 * @file ocn.h
 * @author Alexander S Fox
 * @brief Header file for OCN optimization.
 */

//#TODO I'm worried that our method of choosing a new vertex if the last one is invalid introduces bias. Either choose a random direction to walk, or just try every vertex in random order.

#ifndef OCN_H
#define OCN_H

#include "status.h"
#include "flowgrid.h"

/**
 * @brief Compute the total energy of the flowgrid. Unsafe.
 * @param G Pointer to the FlowGrid. Will not be modified.
 * @param gamma The exponent used in the energy calculation.
 * @return The total energy of the flowgrid.
 */
double ocn_compute_energy(FlowGrid *G, double gamma);

/**
 * @brief Perform a single erosion event on the flowgrid, attempting to change the outflow of a random vertex.
 * Selects a random vertex and a random new direction and attempts to modify the graph accordingly. 
 * If the modification results in a malformed graph, it is undone and the process is retried up to a set number of times.
 * In case of a failure, each vertex is tried 8 different times, once for each possible new direction.
 * If no valid modification is found for a given vertex, a new vertex is selected and the process repeats.
 * If no valid modification is found after trying a set number of vertices, the function exits with a warning.
 * Unsafe! Assumed a well-formed graph.
 * @param G Pointer to the FlowGrid.
 * @param gamma The exponent used in the energy calculation.
 * @param temperature The temperature parameter for the Metropolis-Hastings acceptance criterion.
 * @return Status code indicating success or failure.
 */
Status ocn_single_erosion_event(FlowGrid *G, double gamma, double temperature);

/**
 * @brief Perform multiple erosion events on the flowgrid.
 * @param G Pointer to the FlowGrid.
 * @param energy_history An array to store the energy of the graph after each iteration. Length must be at least niterations.
 * @param niterations The number of erosion events to perform.
 * @param gamma The exponent used in the energy calculation.
 * @param annealing_schedule An array of temperatures (ranging from 0-1) to use for each iteration. Length must be at least niterations.
 * @return Status code indicating success or failure
 */
Status ocn_outer_ocn_loop(FlowGrid *G, uint32_t niterations, double gamma, double *annealing_schedule);

#endif // OCN_H
