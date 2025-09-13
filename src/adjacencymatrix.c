#include <stdint.h>
#include <stdlib.h>
#include <assert.h>
#include <stdbool.h>
#include "adjacencymatrix.h"

static void assert_in_bounds(AdjMat A, int32_t i, int32_t j) {
    assert(i >= 0 && i < A.rows && j >= 0 && j < A.cols);
}

uint8_t get_AdjMat(AdjMat A, int32_t i, int32_t j) {
    assert_inbounds(A, i, j);
    return A.data[i * A.cols + j];
}

void set_AdjMat(AdjMat *A, int32_t i, int32_t j, uint8_t neighbors) {
    assert_inbounds(*A, i, j);
    A->data[i * A->cols + j] = neighbors;
}

NeighborVertices get_neighbor_vertices(AdjMat A, int32_t i, int32_t j) {
    assert_in_bounds(A, i, j);
    uint8_t neighbors = get_AdjMat(A, i, j);
    
    uint8_t bitmask = 1; 
    uint8_t nneighbors = 0;  // number of neighbors found
    int32_t row[8], col[8];  // arrays to hold row and column indices of neighbors
    for (int _ = 0; _ < 8; _++) {
        switch(bitmask & neighbors){
            case 0x01:  // top neighbor
                row[nneighbors] = i - 1; col[nneighbors] = j;
                nneighbors++;
                break;
            case 0x02:  // top-right neighbor
                row[nneighbors] = i - 1; col[nneighbors] = j + 1;
                nneighbors++;
                break;
            case 0x04:  // right neighbor
                row[nneighbors] = i; col[nneighbors] = j + 1;
                nneighbors++;
                break;
            case 0x08:  // bottom-right neighbor
                row[nneighbors] = i + 1; col[nneighbors] = j + 1;
                nneighbors++;
                break;
            case 0x10:  // bottom neighbor
                row[nneighbors] = i + 1; col[nneighbors] = j;
                nneighbors++;
                break;
            case 0x20:  // bottom-left neighbor
                row[nneighbors] = i + 1; col[nneighbors] = j - 1;
                nneighbors++;
                break;
            case 0x40:  // left neighbor
                row[nneighbors] = i; col[nneighbors] = j - 1;
                nneighbors++;
                break;
            case 0x80:  // top-left neighbor
                row[nneighbors] = i - 1; col[nneighbors] = j - 1;
                nneighbors++;
                break;
            default:
                break;
        }
        bitmask <<= 1;
    }
    NeighborVertices neighbor_verts = {
        .n = nneighbors,
        .row = row,
        .col = col
    };

    return neighbor_verts;
}

int32_t *get_clockhand_neighbor(AdjMat A, int32_t i, int32_t j, uint8_t clockhand) {
    // returns the index of the neighbor pointed to by the clockhand (0-7) of the vertex at (i,j).
    uint8_t neighbors = get_AdjMat(A, i, j);
    int32_t neighbor_index[2];
    if (neighbors & (1 << clockhand)) {
        switch(clockhand) {
            case 0: neighbor_index[0] = i - 1; neighbor_index[1] = j; break;       // top
            case 1: neighbor_index[0] = i - 1; neighbor_index[1] = j + 1; break;   // top-right
            case 2: neighbor_index[0] = i; neighbor_index[1] = j + 1; break;       // right
            case 3: neighbor_index[0] = i + 1; neighbor_index[1] = j + 1; break;   // bottom-right
            case 4: neighbor_index[0] = i + 1; neighbor_index[1] = j; break;       // bottom
            case 5: neighbor_index[0] = i + 1; neighbor_index[1] = j - 1; break;   // bottom-left
            case 6: neighbor_index[0] = i; neighbor_index[1] = j - 1; break;       // left
            case 7: neighbor_index[0] = i - 1; neighbor_index[1] = j - 1; break;  // top-left
            default: return NULL; // invalid clockhand
        }
        return neighbor_index;
    }
    return NULL; // no neighbor in that direction
}

int permute_vertices(AdjMat *A, Permutation perm) {
    // permutes a single edge in AdjMat according to the permutation mapping perm.

    uint8_t source = get_AdjMat(*A, perm.i, perm.j);
    
    int validity = 0;
    if (!(source & (1 << perm.oldclockhand))){  // old edge does not exist
        validity += 1;
    }
    if (source & (1 << perm.newclockhand)){ // new edge already exists
        validity += 2;
    }
    if (validity == 0){  // only proceed if permutation is valid
        // 1. The source vertex (remove old edge to target, add new edge to new target)
        // 2. The old target vertex (remove old edge to source)
        // 3. The new target vertex (add new edge to source)
        
        uint8_t sourcemask = (1 << perm.oldclockhand) | (1 << perm.newclockhand);  
        uint8_t oldtargetmask = (1 << perm.oldclockhand) + 4;
        uint8_t newtargetmask = (1 << perm.newclockhand) + 4;

        int32_t *oldtargetidx = get_clockhand_neighbor(*A, perm.i, perm.j, perm.oldclockhand);
        int32_t *newtargetidx = get_clockhand_neighbor(*A, perm.i, perm.j, perm.newclockhand);
        
        int32_t ioldtarget = oldtargetidx[0];
        int32_t inewtarget = newtargetidx[0];
        
        int32_t joldtarget = oldtargetidx[1];
        int32_t jnewtarget = newtargetidx[1];

        uint8_t oldtarget = get_AdjMat(*A, ioldtarget, joldtarget);
        uint8_t newtarget = get_AdjMat(*A, inewtarget, jnewtarget);

        set_AdjMat(A, perm.i, perm.j, source ^ sourcemask);
        set_AdjMat(A, ioldtarget, joldtarget, oldtarget ^ oldtargetmask);
        set_AdjMat(A, inewtarget, jnewtarget, newtarget ^ newtargetmask);
    }
    return validity;
}

void free_AdjMat(AdjMat *A) {
    if (A->data != NULL) {
        free(A->data);
        A->data = NULL;
    }
}




