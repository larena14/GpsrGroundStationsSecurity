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

    uint64_t generateNonce();

private:
    std::mt19937_64 generator;
    std::uniform_int_distribution<uint64_t> distribution;
};



#endif
