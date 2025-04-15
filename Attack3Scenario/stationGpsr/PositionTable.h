//
// Copyright (C) 2013 OpenSim Ltd.
//
// SPDX-License-Identifier: LGPL-3.0-or-later
//


#ifndef __INET_POSITIONTABLE_H
#define __INET_POSITIONTABLE_H

#include <map>
#include <vector>

#include "inet/common/geometry/common/Coord.h"
#include "inet/networklayer/common/L3Address.h"
#include <string>
#include <sstream>

namespace inet {

/**
 * This class provides a mapping between node addresses and their positions.
 */

class INET_API PositionTable
{
  private:
    typedef std::pair<simtime_t, Coord> AddressToPositionMapValue;
    typedef std::map<L3Address, AddressToPositionMapValue> AddressToPositionMap;
    AddressToPositionMap addressToPositionMap;

  public:
    struct PositionWithTimestamp {
        simtime_t timestamp;
        Coord position;

        std::string toString() const {
            std::ostringstream oss;
            oss << "Position: (" << position.x << ", " << position.y << ", " << position.z << "), "
                << "Timestamp: " << timestamp;
            return oss.str();
        }
    };

    PositionTable() {}

    std::vector<L3Address> getAddresses() const;

    bool hasPosition(const L3Address& address) const;
    Coord getPosition(const L3Address& address) const;
    PositionWithTimestamp getPositionWithTimestamp(const L3Address& address) const;
    void setPosition(const L3Address& address, const Coord& coord);

    void removePosition(const L3Address& address);
    void removeOldPositions(simtime_t timestamp);

    void clear();

    simtime_t getOldestPosition() const;

    friend std::ostream& operator<<(std::ostream& o, const PositionTable& t);
};

} // namespace inet

#endif

