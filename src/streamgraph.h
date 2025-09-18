#ifndef STREAMGRAPH_H
#define STREAMGRAPH_H

/**
 * @file streamgraph.h
 * @brief Data structures and functions for managing stream graphs.
 * A stream graph represents a tree of vertices with directed edges,
 * where each vertex has a drained area and flows downstream to another vertex.
 * The graph is stored in a 2D grid layout for spatial representation, and laid out in memory using block tiling.
 * The edges connected to each vertex are represented using an 8-bit integer, where each bit corresponds to one of the 8 possible neighboring directions.
 * The downstream direction is represented as a clock hand (0-7), with 255 indicating no downstream (root).
 */

#include <stdint.h>
#include <stdbool.h>
#include "returncodes.h"

// Type definitions
typedef int32_t drainedarea_t;
typedef uint16_t cartidx_t;
typedef uint32_t linidx_t;
typedef uint8_t localedges_t;
typedef uint8_t clockhand_t;

#define IS_ROOT (clockhand_t)255
#define TILE_SIZE 2

typedef struct {
    drainedarea_t drained_area;
    linidx_t adown;
    localedges_t edges;
    clockhand_t downstream;
    uint8_t visited;
} Vertex;

typedef struct {
    cartidx_t m, n;
    cartidx_t i_root, j_root;
    double energy;
    Vertex *vertices;
} StreamGraph;

Status sg_create(StreamGraph *G, cartidx_t m, cartidx_t n);
Status sg_destroy(StreamGraph *G);

linidx_t sg_cart_to_lin(cartidx_t row, cartidx_t col, cartidx_t m, cartidx_t n);

Status sg_get_cart_safe(StreamGraph G, Vertex *out, cartidx_t row, cartidx_t col);
Status sg_get_cart(StreamGraph G, Vertex *out, cartidx_t row, cartidx_t col);
Status sg_set_cart_safe(Vertex vert, StreamGraph *G, cartidx_t row, cartidx_t col);
Status sg_set_cart(Vertex vert, StreamGraph *G, cartidx_t row, cartidx_t col);

Status sg_get_lin_safe(StreamGraph G, Vertex *out, linidx_t a);
Status sg_get_lin(StreamGraph G, Vertex *out, linidx_t a);
Status sg_set_lin_safe(Vertex vert, StreamGraph *G, linidx_t a);
Status sg_set_lin(Vertex vert, StreamGraph *G, linidx_t a);

Status sg_change_vertex_outflow(linidx_t a, StreamGraph *G, clockhand_t down_new);

Status sg_increment_downstream_safe(drainedarea_t inc, StreamGraph *G, linidx_t a, uint8_t ncalls);
Status sg_increment_downstream(drainedarea_t inc, StreamGraph *G, linidx_t a, uint8_t ncalls);

Status sg_single_erosion_event(StreamGraph *G);
Status sg_outer_ocn_loop(uint32_t niterations, StreamGraph *G, drainedarea_t tol);

#endif // STREAMGRAPH_H