/*
 * NonceGenerator.cc
 *
 *  Created on: Dec 20, 2024
 *      Author: larena
 */

#include "MovingAverageFilter.h"
#include <vector>
#include "inet/common/geometry/common/Coord.h"


MovingAverageFilter::MovingAverageFilter(size_t windowSize) : buffer(windowSize) {}

void MovingAverageFilter::update(const inet::Coord& newPos, double timestamp) {
    buffer[index] = {newPos, timestamp};
    index = (index + 1) % buffer.size();
    if(!filled && index == 0) filled = true;
}

inet::Coord MovingAverageFilter::filter(const inet::Coord& newPos, double timestamp){

    if (index == 0 && !filled)
        return newPos; // Se il buffer è vuoto, restituisco direttamente la nuova posizione

    inet::Coord weightedAvg = newPos;
    double weightSum = 1.0; // Inizialmente la nuova posizione ha peso 1
    double latestTimestamp = timestamp; // Consideriamo il timestamp della nuova posizione come riferimento

    size_t count = filled ? buffer.size() : index;
    for (size_t i = 0; i < count; ++i) {
        double age = latestTimestamp - buffer[i].second; // Differenza temporale tra nuovo dato e buffer
        double weight = 1.0 / (1.0 + age);  // Peso inversamente proporzionale al tempo

        weightedAvg.x += buffer[i].first.x * weight;
        weightedAvg.y += buffer[i].first.y * weight;
        weightedAvg.z += buffer[i].first.z * weight;
        weightSum += weight;
    }

    // Normalizzazione
    weightedAvg.x /= weightSum;
    weightedAvg.y /= weightSum;
    weightedAvg.z /= weightSum;

    return weightedAvg;
}



