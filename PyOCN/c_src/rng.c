#ifndef RNG_C
#define RNG_C

#include <stdlib.h>
#include <time.h>
#include <stdint.h>
#include "rng.h"

static inline uint32_t rotate_left(const uint32_t x, int k){
    return (x << k) | (x >> (32 - k));
}

void rng_seed(rng_state_t *rng, uint32_t seed) {
    // splitmix32 to initialize the 4 state variables
    uint32_t z = seed;
    for (int i = 0; i < 4; i++) {
        z += 0x9e3779b9; z ^= z >> 16;
        z *= 0x21f0aaad; z ^= z >> 15;
        z *= 0x735a2d97; z ^= z >> 15;
        rng->s[i] = z;
    }
}

void rng_seed_random(rng_state_t *rng) {
    rng_seed(rng, (uint32_t)time(NULL));
}

uint32_t rng_randint32(rng_state_t *rng){
    uint32_t result = rng->s[0] + rng->s[3];
    
    uint32_t t = rng->s[1] << 9;

    rng->s[2] ^= rng->s[0];
    rng->s[3] ^= rng->s[1];
    rng->s[1] ^= rng->s[2];
    rng->s[0] ^= rng->s[3];

    rng->s[2] ^= t;
    rng->s[3] = rotate_left(rng->s[3], 11);

    return result;
}

double rng_uniformdouble(rng_state_t *rng){
    // generate a uniform random double in [0, 1)
    return (double)rng_randint32(rng) / (double)UINT32_MAX;
}

#endif // RNG_C
