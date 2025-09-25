#ifndef STREAMGRAPH_H
#define STREAMGRAPH_H

#include <stdint.h>

// basic types
typedef uint32_t drainedarea_t;
typedef uint16_t cartidx_t;
typedef struct {cartidx_t row, col;} CartPair;
typedef uint32_t linidx_t;
typedef uint8_t localedges_t;
typedef uint8_t clockhand_t;
extern clockhand_t IS_ROOT;

typedef struct {
    drainedarea_t drained_area;  // 4B
    linidx_t adown;  // 4B
    localedges_t edges;  // 1B
    clockhand_t downstream;  // 1B
    uint8_t visited;  // 1B
} Vertex;

typedef struct {
    CartPair dims;
    CartPair root;
    double energy;
    Vertex *vertices;
} StreamGraph;

/** @brief Coordinate Transformation */
linidx_t sg_cart_to_lin(CartPair coords, CartPair dims);

/** @brief Coordinate Transformation */
CartPair sg_lin_to_cart(linidx_t a, CartPair dims);

/**
 * @brief Given a linear index and a clockhand direction, find the linear index of the vertex in that direction.
 * @param a_down Pointer to store the resulting linear index.
 * @param a The starting linear index.
 * @param down The clockhand direction to move in.
 * @param dims The dimensions of the graph.
 * @return Status code indicating success or failure
 */
Status sg_clockhand_to_lin_safe(linidx_t *a_down, linidx_t a, clockhand_t down, CartPair dims);

/**
 * @brief Get the vertex at the given Cartesian coordinates safely.
 * @param out Pointer to store the resulting vertex.
 * @param G Pointer to the StreamGraph.
 * @param coords The Cartesian coordinates of the vertex to retrieve.
 * @return Status code indicating success or failure
 */
Status sg_get_cart_safe(Vertex *out, StreamGraph *G, CartPair coords);

/**
 * @brief Get the vertex at the given Cartesian coordinates unsafely.
 * @param G Pointer to the StreamGraph.
 * @param coords The Cartesian coordinates of the vertex to retrieve.
 * @return The vertex at the specified coordinates.
 */
Vertex sg_get_cart(StreamGraph *G, CartPair coords);

/**
 * @brief Set the vertex at the given Cartesian coordinates safely.
 * @param G Pointer to the StreamGraph.
 * @param vert The vertex to use to update the graph with.
 * @param coords The Cartesian coordinates where the vertex should be set.
 * @return Status code indicating success or failure
 */
Status sg_set_cart_safe(StreamGraph *G, Vertex vert, CartPair coords);

/**
 * @brief Set the vertex at the given Cartesian coordinates unsafely.
 * @param G Pointer to the StreamGraph.
 * @param vert The vertex to use to update the graph with.
 * @param coords The Cartesian coordinates where the vertex should be set.
 */
void sg_set_cart(StreamGraph *G, Vertex vert, CartPair coords);

/**
 * @brief Get the vertex at the given linear index safely.
 * @param out Pointer to store the resulting vertex.
 * @param G Pointer to the StreamGraph.
 * @param a The linear index of the vertex to retrieve.
 * @return Status code indicating success or failure
 */
Status sg_get_lin_safe(Vertex *out, StreamGraph *G, linidx_t a);

/**
 * @brief Get the vertex at the given linear index unsafely.
 * @param G Pointer to the StreamGraph.
 * @param a The linear index of the vertex to retrieve.
 * @return The vertex at the specified linear index.
 */
Vertex sg_get_lin(StreamGraph *G, linidx_t a);

/**
 * @brief Set the vertex at the given linear index safely.
 * @param G Pointer to the StreamGraph.
 * @param vert The vertex to use to update the graph with.
 * @param a The linear index where the vertex should be set.
 * @return Status code indicating success or failure
 */
Status sg_set_lin_safe(StreamGraph *G, Vertex vert, linidx_t a);

/**
 * @brief Set the vertex at the given linear index unsafely.
 * @param G Pointer to the StreamGraph.
 * @param vert The vertex to use to update the graph with.
 * @param a The linear index where the vertex should be set.
 */
void sg_set_lin(StreamGraph *G, Vertex vert, linidx_t a);

/**
 * @brief Create an empty streamgraph with given dimensions safely.
 * @param dims The dimensions of the graph (rows, cols).
 * @return The created StreamGraph. Returns NULL if dimensions are invalid or memory allocation fails.
 */
StreamGraph *sg_create_empty_safe(CartPair dims);

/**
 * @brief Safely destroy a streamgraph, freeing its resources.
 * @param G Pointer to the StreamGraph to destroy.
 * @return Status code indicating success or failure
 */
Status sg_destroy_safe(StreamGraph *G);

/**
 * @brief Change the outflow direction of a vertex safely.
 * @param G Pointer to the StreamGraph.
 * @param a The linear index of the vertex to modify.
 * @param down_new The new clockhand direction for the vertex's outflow.
 * @return Status code indicating success or failure. Returns MALFORMED_GRAPH_WARNING if the change would malform the graph in an immediately obvious way. Does not check for large cycles or root access.
 */
Status sg_change_vertex_outflow(StreamGraph *G, linidx_t a, clockhand_t down_new);

/** 
 * @brief Follow the downstream path from a given vertex, marking each vertex as visited.
 * @param G Pointer to the StreamGraph.
 * @param a The linear index of the starting vertex.
 * @param ncalls A unique identifier for this traversal to mark visited vertices.
 * @return Status code indicating success or failure. Returns MALFORMED_GRAPH_WARNING if a cycle is detected.
 */
Status sg_flow_downstream_safe(StreamGraph *G, linidx_t a, uint8_t ncalls);

/**
 * @brief Display the streamgraph in the terminal using ASCII or UTF-8 characters.
 * @param G Pointer to the StreamGraph to display.
 * @param use_utf8 If true, use UTF-8 characters for better visuals; otherwise, use ASCII.
 */
void sg_display(StreamGraph *G, bool use_utf8);

/**
 * @brief Create a simple test graph for demonstration purposes.
 * O  O  O  O
 * |  |  |  |
 * O  O  O  O
 * |  |  |  |
 * O  O  O  O
 * |  |  |  |
 * O--O--O--X
 * @return Pointer to the newly created StreamGraph.
 */
StreamGraph *sg_make_test_graph();

#endif // STREAMGRAPH_H