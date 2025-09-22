#ifndef RETURNCODES_H
#define RETURNCODES_H

#include <stdint.h>

typedef uint8_t Status;

#define SUCCESS (Status)0
// vertex return codes
// valid values for .downstream are 0-7
// use remaining values for error codes
#define OOB_ERROR (Status)1
#define NO_EDGE_WARNING (Status)2
#define NULL_POINTER_ERROR (Status)3
#define SWAP_WARNING (Status)4
#define MALFORMED_GRAPH_WARNING (Status)5
#define CYCLE_WARNING (Status)6

#endif // RETURNCODES_H