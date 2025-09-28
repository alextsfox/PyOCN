#ifndef RNG_C
#define RNG_C

#include <stdlib.h>
#include <time.h>
#include "rng.h"

void rng_seed(unsigned int seed) {
    srand(seed);
}

void rng_seed_random(void) {
    srand((unsigned int)time(NULL));
}

#endif // RNG_C
