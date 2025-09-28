#ifndef STATUS_C
#define STATUS_C

#include <stdint.h>

#include "status.h"

Status SUCCESS = 0;
Status EROSION_FAILURE = 1;
Status OOB_ERROR = 2;
Status NULL_POINTER_ERROR = 3;
Status SWAP_WARNING = 4;
Status MALFORMED_GRAPH_WARNING = 5;

#endif // STATUS_C
