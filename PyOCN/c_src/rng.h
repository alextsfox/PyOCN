/**
 * @file rng.h
 * @author Alexander S Fox
 * @brief Header file for random number generator functions.
 */
#ifndef RNG_H
#define RNG_H

#include <stdint.h>
#include "status.h"

typedef struct{
    uint32_t s[4];  // State array for xorshift128+
} rng_state_t;

/**
 * @brief Seed the random number generator.
 * Uses splitmix32 with constants from https://stackoverflow.com/a/52056161
 * @param state Pointer to the rng_state_t structure to be initialized.
 * @param seed The seed value to initialize the state.
 */
void rng_seed(rng_state_t *state, uint32_t seed);

/**
 * @brief Seed the random number generator. Uses the current time as the seed.
 * @param state Pointer to the rng_state_t structure to be initialized.
 */
void rng_seed_random(rng_state_t *state);

/**
 * @brief Generate the next random number using xoshiro128+ algorithm. 
 * From https://prng.di.unimi.it/xoshiro128plus.c, David Blackman, Sebastiano Vigna.
 * @param state Pointer to the rng_state_t structure holding the current state.
 * @return A pseudo-random 32-bit unsigned integer.
 */
uint32_t rng_randint32(rng_state_t *state);

/**
 * @brief Generate a uniform random number in the range [0, 1].
 * @param state Pointer to the rng_state_t structure holding the current state.
 * @return A pseudo-random double in the range [0, 1).
 */
double rng_uniformdouble(rng_state_t *state);

#endif // RNG_H
