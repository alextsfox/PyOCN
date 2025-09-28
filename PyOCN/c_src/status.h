/**
 * @file status.h
 * @author Alexander S Fox
 * @brief Header file for status/return codes used in the library.
 */

#ifndef STATUS_H
#define STATUS_H

#include <stdint.h>

typedef uint8_t Status;
extern Status SUCCESS;
extern Status EROSION_FAILURE;
extern Status OOB_ERROR;
extern Status NULL_POINTER_ERROR;
extern Status SWAP_WARNING;
extern Status MALFORMED_GRAPH_WARNING;

#endif // STATUS_H
