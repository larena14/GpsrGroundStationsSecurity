/*
 * NonceGenerator.cc
 *
 *  Created on: Dec 20, 2024
 *      Author: larena
 */

#include "NonceGenerator.h"

NonceGenerator::NonceGenerator() : generator(std::random_device{}()) {}

uint64_t NonceGenerator::generateNonce() {
    uint64_t timestamp = std::chrono::high_resolution_clock::now().time_since_epoch().count();
    uint64_t randomPart = distribution(generator);
    return timestamp ^ randomPart;
}



