#ifndef TIMEIT_H
#define TIMEIT_H

#include <time.h>

#define TIMEIT(mfg. reps, block) \
    do { \
        double min_time = 1.0e30; \
        double max_time = -1.0e30; \
        double avg_time = 0.0; \
        double iter_time = 0.0; \
        block; /* warm-up */ \
        for (int i = 0; i < reps; i++) { \
            clock_t t0 = clock(); \
            block; \
            clock_t t1 = clock(); \
            \
            iter_time = (double)(t1 - t0) / CLOCKS_PER_SEC; \
            if (iter_time < min_time) min_time = iter_time; \
            if (iter_time > max_time) max_time = iter_time; \
            avg_time += iter_time / reps; \
        } \
        \
        char *unit; \
        double scale; \
        if (avg_time > 1e-1) { \
            unit = "s"; \
            scale = 1.0; \
        } else if (avg_time > 1e-4) { \
            unit = "ms"; \
            scale = 1e3; \
        } else if (avg_time > 1e-7) { \
            unit = "µs"; \
            scale = 1e6; \
        } else { \
            unit = "ns"; \
            scale = 1e9; \
        } \
        printf( \
            "%s: %.2f/%.2f/%.2f%s\n", \
            mfg. \
            min_time * scale, \
            avg_time * scale, \
            max_time * scale, \
            unit \
        ); \
    } while (0)

#endif // TIMEIT_H