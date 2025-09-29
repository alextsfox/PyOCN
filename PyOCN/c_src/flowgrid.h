/**
 * @file flowgrid.h
 * @author Alexander S Fox
 * @brief Header file for FlowGrid structure and related functions.
 */

#ifndef FLOWGRID_H
#define FLOWGRID_H

#include <stdint.h>

// basic types
typedef double drainedarea_t;
typedef uint16_t cartidx_t;
typedef struct {cartidx_t row, col;} CartPair;
typedef uint32_t linidx_t;
typedef uint8_t localedges_t;
typedef uint8_t clockhand_t;
extern clockhand_t IS_ROOT;

/**
 * @brief Vertex structure representing a node in the flow grid.
 * 
 * - drained_area: The area drained by this vertex.
 * 
 * - adown: The linear index of the downstream vertex.
 * 
 * - edges: A bitmask representing the presence of edges in the 8 possible directions.
 * 
 * - downstream: The clockhand direction of the downstream flow (0-7) or IS_ROOT (255) if it is a root node.
 * 
 * - visited: A flag used for traversal algorithms
 */
typedef struct {
    drainedarea_t drained_area;  // 8B
    linidx_t adown;  // 4B
    localedges_t edges;  // 1B
    clockhand_t downstream;  // 1B
    uint8_t visited;  // 1B
} Vertex;

/**
 * @brief FlowGrid structure representing the entire flow grid.
 * 
 * - dims: The dimensions of the grid (rows, cols).
 * 
 * 
 * - energy: The energy of the flow grid.
 * 
 * - resolution: the side length of each pixel in meters.
 * 
 * - nroots: The number of root nodes in the grid. Used for more efficient energy calculations when there is only one root.
 * 
 * - vertices: A pointer to an array of Vertex structures representing the nodes in the grid.
 */
typedef struct {
    CartPair dims;
    double energy;
    double resolution;
    uint16_t nroots;
    Vertex *vertices;
} FlowGrid;

/** @brief Coordinate Transformation */
linidx_t fg_cart_to_lin(CartPair coords, CartPair dims);

/** @brief Coordinate Transformation */
CartPair fg_lin_to_cart(linidx_t a, CartPair dims);

/**
 * @brief Given a linear index and a clockhand direction, find the linear index of the vertex in that direction.
 * @param a_down Pointer to store the resulting linear index.
 * @param a The starting linear index.
 * @param down The clockhand direction to move in.
 * @param dims The dimensions of the graph.
 * @return Status code indicating success or failure
 */
Status fg_clockhand_to_lin_safe(linidx_t *a_down, linidx_t a, clockhand_t down, CartPair dims);

/**
 * @brief Get the vertex at the given Cartesian coordinates safely.
 * @param out Pointer to store the resulting vertex.
 * @param G Pointer to the FlowGrid.
 * @param coords The Cartesian coordinates of the vertex to retrieve.
 * @return Status code indicating success or failure
 */
Status fg_get_cart_safe(Vertex *out, FlowGrid *G, CartPair coords);

/**
 * @brief Get the vertex at the given Cartesian coordinates unsafely.
 * @param G Pointer to the FlowGrid.
 * @param coords The Cartesian coordinates of the vertex to retrieve.
 * @return The vertex at the specified coordinates.
 */
Vertex fg_get_cart(FlowGrid *G, CartPair coords);

/**
 * @brief Set the vertex at the given Cartesian coordinates safely.
 * @param G Pointer to the FlowGrid.
 * @param vert The vertex to use to update the graph with.
 * @param coords The Cartesian coordinates where the vertex should be set.
 * @return Status code indicating success or failure
 */
Status fg_set_cart_safe(FlowGrid *G, Vertex vert, CartPair coords);

/**
 * @brief Set the vertex at the given Cartesian coordinates unsafely.
 * @param G Pointer to the FlowGrid.
 * @param vert The vertex to use to update the graph with.
 * @param coords The Cartesian coordinates where the vertex should be set.
 */
void fg_set_cart(FlowGrid *G, Vertex vert, CartPair coords);

/**
 * @brief Get the vertex at the given linear index safely.
 * @param out Pointer to store the resulting vertex.
 * @param G Pointer to the FlowGrid.
 * @param a The linear index of the vertex to retrieve.
 * @return Status code indicating success or failure
 */
Status fg_get_lin_safe(Vertex *out, FlowGrid *G, linidx_t a);

/**
 * @brief Get the vertex at the given linear index unsafely.
 * @param G Pointer to the FlowGrid.
 * @param a The linear index of the vertex to retrieve.
 * @return The vertex at the specified linear index.
 */
Vertex fg_get_lin(FlowGrid *G, linidx_t a);

/**
 * @brief Set the vertex at the given linear index safely.
 * @param G Pointer to the FlowGrid.
 * @param vert The vertex to use to update the graph with.
 * @param a The linear index where the vertex should be set.
 * @return Status code indicating success or failure
 */
Status fg_set_lin_safe(FlowGrid *G, Vertex vert, linidx_t a);

/**
 * @brief Set the vertex at the given linear index unsafely.
 * @param G Pointer to the FlowGrid.
 * @param vert The vertex to use to update the graph with.
 * @param a The linear index where the vertex should be set.
 */
void fg_set_lin(FlowGrid *G, Vertex vert, linidx_t a);

/**
 * @brief Create an empty flowgrid with given dimensions safely.
 * @param dims The dimensions of the graph (rows, cols).
 * @return The created FlowGrid. Returns NULL if dimensions are invalid or memory allocation fails.
 */
FlowGrid *fg_create_empty_safe(CartPair dims);

/**
 * @brief Create a deep copy of a flowgrid safely.
 * @param G Pointer to the FlowGrid to copy.
 * @return A deep copy of the FlowGrid. Returns NULL if G is NULL or memory allocation fails.
 */
FlowGrid *fg_copy_safe(FlowGrid *G);

/**
 * @brief Safely destroy a flowgrid, freeing its resources.
 * @param G Pointer to the FlowGrid to destroy.
 * @return Status code indicating success or failure
 */
Status fg_destroy_safe(FlowGrid *G);

/**
 * @brief Change the outflow direction of a vertex safely.
 * @param G Pointer to the FlowGrid.
 * @param a The linear index of the vertex to modify.
 * @param down_new The new clockhand direction for the vertex's outflow.
 * @return Status code indicating success or failure. Returns MALFORMED_GRAPH_WARNING if the change would malform the graph in an immediately obvious way. Does not check for large cycles or root access.
 */
Status fg_change_vertex_outflow(FlowGrid *G, linidx_t a, clockhand_t down_new);

/** 
 * @brief Follow the downstream path from a given vertex, marking each vertex as visited.
 * @param G Pointer to the FlowGrid.
 * @param a The linear index of the starting vertex.
 * @param ncalls A unique identifier for this traversal to mark visited vertices.
 * @return Status code indicating success or failure. Returns MALFORMED_GRAPH_WARNING if a cycle is detected.
 */
Status fg_flow_downstream_safe(FlowGrid *G, linidx_t a, uint8_t ncalls);

/**
 * @brief Display the flowgrid in the terminal using ASCII or UTF-8 characters.
 * @param G Pointer to the FlowGrid to display.
 * @param use_utf8 If true, use UTF-8 characters for better visuals; otherwise, use ASCII.
 */
void fg_display(FlowGrid *G, bool use_utf8);

#endif // FLOWGRID_H
