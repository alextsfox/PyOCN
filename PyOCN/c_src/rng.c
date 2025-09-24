#ifndef RNG_C
#define RNG_C

#include <stdlib.h>

void seed(unsigned int seed) {
    if (seed == NULL) seed = (unsigned int)time(NULL);
    srand(seed);
}

#endif // RNG_C