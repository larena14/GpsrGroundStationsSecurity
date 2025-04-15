/*
 * NonceGenerator.h
 *
 *  Created on: Dec 20, 2024
 *      Author: larena
 */

#ifndef STATIONGPSR_NONCEGENERATOR_H_
#define STATIONGPSR_NONCEGENERATOR_H_


#include <random>
#include <chrono>

class NonceGenerator {
public:
    NonceGenerator();

    // Metodo per generare il nonce
    uint64_t generateNonce();

private:
    std::mt19937_64 generator;                  // Generatore di numeri casuali
    std::uniform_int_distribution<uint64_t> distribution; // Distribuzione uniforme per generare la parte casuale
};



#endif /* STATIONGPSR_NONCEGENERATOR_H_ */
