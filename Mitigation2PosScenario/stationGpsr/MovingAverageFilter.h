/*
 * NonceGenerator.h
 *
 *  Created on: Dec 20, 2024
 *      Author: larena
 */

#ifndef STATIONGPSR_MOVINGAVERAGEFILTER_H_
#define STATIONGPSR_MOVINGAVERAGEFILTER_H_


#include <vector>
#include "inet/common/geometry/common/Coord.h"


class MovingAverageFilter {

public:
    MovingAverageFilter() : MovingAverageFilter(5) {}
    MovingAverageFilter(size_t windowSize);

    void update(const inet::Coord& newPos, double timestamp);
    inet::Coord filter(const inet::Coord& newPos, double timestamp);

private:
    std::vector<std::pair<inet::Coord, double>> buffer;
    size_t index = 0;
    bool filled = false;
};

#endif /* STATIONGPSR_MOVINGAVERAGEFILTER_H_ */
