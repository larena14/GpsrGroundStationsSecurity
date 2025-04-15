/*
 * NonceGenerator.cc
 *
 *  Created on: Dec 20, 2024
 *      Author: larena
 */

#include "DistanceBuffer.h"
#include <vector>
#include <limits>
#include <numeric>
#include <cmath>

using namespace std;

DistanceBuffer::DistanceBuffer(size_t windowSize) : buffer(windowSize) {}

void DistanceBuffer::update(double distance) {
    if (distance > MAX_DISTANCE)
        MAX_DISTANCE = distance;
    buffer[index] = distance;
    index = (index + 1) % buffer.size();
    if(!filled && index == 0) filled = true;
}

double DistanceBuffer::avgDistance(){
    size_t count = filled ? buffer.size() : index;
    if (count == 0) return std::numeric_limits<double>::max();
    return accumulate(buffer.begin(), buffer.begin() + count, 0.0) / count;
}

double DistanceBuffer::stdDevDistance() {
    size_t count = filled ? buffer.size() : index;
    if (count == 0) return 0.0;

    double mean = avgDistance();
    double varianceSum = 0.0;

    for (size_t i = 0; i < count; ++i) {
        double diff = buffer[i] - mean;
        varianceSum += diff * diff;
    }

    return std::sqrt(varianceSum / count);
}

double DistanceBuffer::getMaxDistance(){
    return MAX_DISTANCE;
}




