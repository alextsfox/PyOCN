#ifndef STREAMGRAPH_H
#define STREAMGRAPH_H

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