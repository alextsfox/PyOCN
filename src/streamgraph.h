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
    drainedarea_t drained_area;
    linidx_t adown;
    localedges_t edges;
    clockhand_t downstream;
    uint8_t visited;
} Vertex;

typedef struct {
    CartPair dims;
    CartPair root;
    double energy;
    Vertex *vertices;
} StreamGraph;


#endif // STREAMGRAPH_H