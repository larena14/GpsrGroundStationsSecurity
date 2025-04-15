/*
 * NonceGenerator.h
 *
 *  Created on: Dec 20, 2024
 *      Author: larena
 */

#ifndef STATIONGPSR_DISTANCEBUFFER_H_
#define STATIONGPSR_DISTANCEBUFFER_H_


#include <vector>
#include <numeric>
#include <cmath>
#include <numeric>


using namespace std;

class DistanceBuffer {

public:
    DistanceBuffer() : DistanceBuffer(100) {}
    DistanceBuffer(size_t windowSize);

    void update(double distance);
    double avgDistance();
    double stdDevDistance();
    double getMaxDistance();

private:
    std::vector<double> buffer;
    size_t index = 0;
    bool filled = false;
    double MAX_DISTANCE = 0.1;
};

#endif /* STATIONGPSR_DISTANCEBUFFER_H_ */
