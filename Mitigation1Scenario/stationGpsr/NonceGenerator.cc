/*
 * NonceGenerator.cc
 *
 *  Created on: Dec 20, 2024
 *      Author: larena
 */

#include "NonceGenerator.h"

// Costruttore: inizializza il generatore di numeri casuali
NonceGenerator::NonceGenerator() : generator(std::random_device{}()) {}

// Metodo per generare il nonce
uint64_t NonceGenerator::generateNonce() {
    uint64_t timestamp = std::chrono::high_resolution_clock::now().time_since_epoch().count();
    uint64_t randomPart = distribution(generator); // Genera la parte casuale
    return timestamp ^ randomPart; // Combina timestamp e parte casuale per formare il nonce
}



