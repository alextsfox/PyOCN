/**
 * @file rng.h
 * @author Alexander S Fox
 * @brief Header file for random number generator functions.
 */
#ifndef RNG_H
#define RNG_H

/**
 * @brief Seed the random number generator.
 * @param seed The seed value to initialize the random number generator.
 */
void rng_seed(unsigned int seed);

/**
 * @brief Seed the random number generator. Uses the current time as the seed.
 */
void rng_seed_random(void);

#endif // RNG_H
