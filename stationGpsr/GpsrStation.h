//
// Copyright (C) 2013 OpenSim Ltd.
//
// SPDX-License-Identifier: LGPL-3.0-or-later
//


#ifndef __INET_GPSRSTATION_H
#define __INET_GPSRSTATION_H

#include "inet/common/ModuleRefByPar.h"
#include "inet/common/geometry/common/Coord.h"
#include "inet/common/packet/Packet.h"
#include "inet/mobility/contract/IMobility.h"
#include "inet/networklayer/contract/IL3AddressType.h"
#include "inet/networklayer/contract/INetfilter.h"
#include "inet/networklayer/contract/IRoutingTable.h"
#include "inet/routing/base/RoutingProtocolBase.h"
#include "inet/routing/gpsr/Gpsr_m.h"
#include "PositionTable.h"
#include "inet/transportlayer/udp/UdpHeader_m.h"
#include "messages/Station_m.h"
#include <map>
#include <vector>
#include <queue>
#include <utility>

namespace inet {

/**
 * This class implements the Greedy Perimeter Stateless Routing for Wireless Networks.
 * The implementation supports both GG and RNG planarization algorithms.
 *
 * For more information on the routing algorithm, see the GPSR paper
 * http://www.eecs.harvard.edu/~htk/publication/2000-mobi-karp-kung.pdf
 */
// TODO optimize internal data structures for performance to use less lookups and be more prepared for routing a packet
// TODO implement position piggybacking that is all packets should carry the position of the sender, all packets act as a beacon and reset beacon timer
// TODO implement promiscuous mode, all receivers should process all packets with respect to neighbor positions
// KLUDGE implement position registry protocol instead of using a global variable
class INET_API GpsrStation : public RoutingProtocolBase, public cListener, public NetfilterBase::HookBase
{
    enum SendMode {
            FIXED_MAXIMUM = 0,
            INTERVAL = 1
        };
  private:
    // GPSR parameters
    GpsrPlanarizationMode planarizationMode = static_cast<GpsrPlanarizationMode>(-1);
    const char *interfaces = nullptr;
    simtime_t beaconInterval;
    simtime_t maxJitter;
    simtime_t neighborValidityInterval;
    bool displayBubbles;

    // context
    cModule *host = nullptr;
    opp_component_ptr<IMobility> mobility;
    IL3AddressType *addressType = nullptr;
    ModuleRefByPar<IInterfaceTable> interfaceTable;
    const char *outputInterface = nullptr;
    ModuleRefByPar<IRoutingTable> routingTable; // TODO delete when necessary functions are moved to interface table
    ModuleRefByPar<INetfilter> networkProtocol;
    static PositionTable globalPositionTable; // KLUDGE implement position registry protocol

    // packet size
    int positionByteLength = -1;

    // internal
    cMessage *beaconTimer = nullptr;
    cMessage *purgeNeighborsTimer = nullptr;
    PositionTable neighborPositionTable;

    int n;
    int m;
    bool isMalicious;
    bool isStation;

    cMessage *sendTimerInterval = nullptr;


    static double timeInPerimeter;
    static double timeInGreedy;
    static bool printed;
    static std::map<std::string, int> mapPacketHops;
    static std::map<int, L3Address> stationMap;
    static std::map<std::string, L3Address> hostnameMap;

    cMessage *purgeStationTimer = nullptr;
    cMessage *purgeDestinationTimer = nullptr;
    PositionTable quadNodesPositionTable;
    PositionTable destinationsPositionTable;
    int currentStation = -1;
    simtime_t positionValidityInterval;
    simtime_t destinationValidityInterval;

    std::multimap<L3Address, Packet *> targetAddressToDelayedPackets;

    std::queue<std::string> receivers;

    L3Address destinationAddress;
    bool destinationSetted = false;

    simtime_t minInterval;
    simtime_t maxInterval;
    SendMode sendMode = static_cast<SendMode>(-1);

    static std::map<std::string, L3Address> possibleReceiversMap;

    static int msgNumber;
    static int totalMsg;
    static int totalArrivedMsg;
    static std::map<std::string, std::pair<simtime_t, int>> packetInfo;

  public:
    GpsrStation();
    virtual ~GpsrStation();

  protected:
    // module interface
    virtual int numInitStages() const override { return NUM_INIT_STAGES; }
    void initialize(int stage) override;
    void handleMessageWhenUp(cMessage *message) override;
    void finish() override;

  private:
    int getStationIndex(Coord pos);

    const Ptr<StationNotice> createStationNotice(bool deregisterFlag);
    void sendStationNotice(const Ptr<StationNotice>& res, L3Address dest, const char* moduleName);
    void processStationNotice(Packet *packet);

    const Ptr<PositionRequest> createPositionRequest(L3Address address);
    void sendPositionRequest(const Ptr<PositionRequest>& res, L3Address dest, const char* moduleName);
    void processPositionRequest(Packet *packet);

    const Ptr<PositionResponse> createPositionResponse(L3Address address);
    const Ptr<PositionResponse> createPositionResponseFromAnotherStation(L3Address address, Coord position, simtime_t timestamp);
    void sendPositionResponse(const Ptr<PositionResponse>& res, L3Address dest, const char* moduleName);
    void processPositionResponse(Packet *packet);

    const Ptr<S2SPositionRequest> createS2SPositionRequest(L3Address applicant, const char* applicantModuleName, L3Address address);
    void sendS2SPositionRequest(const Ptr<S2SPositionRequest>& res, L3Address dest, const char* moduleName);
    void processS2SPositionRequest(Packet *packet);

    const Ptr<S2SPositionResponse> createS2SPositionResponse(L3Address applicant, const char* applicantModuleName, L3Address address);
    void sendS2SPositionResponse(const Ptr<S2SPositionResponse>& res, L3Address dest, const char* moduleName);
    void processS2SPositionResponse(Packet *packet);

    void schedulePurgeStationTimer();
    void processPurgeStationTimer();

    void schedulePurgeDestinationTimer();
    void processPurgeDestinationTimer();

    void processSendTimer();
    void processSendTimerInterval();

    // handling messages
    void processSelfMessage(cMessage *message);
    void processMessage(cMessage *message);

    // handling beacon timers
    void scheduleBeaconTimer();
    void processBeaconTimer();

    // handling purge neighbors timers
    void schedulePurgeNeighborsTimer();
    void processPurgeNeighborsTimer();

    // handling UDP packets
    void sendUdpPacket(Packet *packet);
    void processUdpPacket(Packet *packet);

    Coord generateFalsePosition();

    // handling beacons
    const Ptr<GpsrBeacon> createBeacon();
    void sendBeacon(const Ptr<GpsrBeacon>& beacon);
    void processBeacon(Packet *packet);

    // handling packets
    GpsrOption *createGpsrOption(L3Address destination);
    int computeOptionLength(GpsrOption *gpsrOption);
    void setGpsrOptionOnNetworkDatagram(Packet *packet, const Ptr<const NetworkHeaderBase>& networkHeader, GpsrOption *gpsrOption);

    // returns nullptr if not found
    GpsrOption *findGpsrOptionInNetworkDatagramForUpdate(const Ptr<NetworkHeaderBase>& networkHeader);
    const GpsrOption *findGpsrOptionInNetworkDatagram(const Ptr<const NetworkHeaderBase>& networkHeader) const;

    // throws an error when not found
    GpsrOption *getGpsrOptionFromNetworkDatagramForUpdate(const Ptr<NetworkHeaderBase>& networkHeader);
    const GpsrOption *getGpsrOptionFromNetworkDatagram(const Ptr<const NetworkHeaderBase>& networkHeader) const;

    // configuration
    void configureInterfaces();

    // position
    Coord lookupPositionInGlobalRegistry(const L3Address& address) const;
    void storePositionInGlobalRegistry(const L3Address& address, const Coord& position) const;
    void storeSelfPositionInGlobalRegistry() const;
    Coord computeIntersectionInsideLineSegments(Coord& begin1, Coord& end1, Coord& begin2, Coord& end2) const;
    Coord getNeighborPosition(const L3Address& address) const;

    // angle
    double getVectorAngle(Coord vector) const;
    double getNeighborAngle(const L3Address& address) const;

    // address
    std::string getHostName() const;
    L3Address getSelfAddress() const;
    L3Address getSenderNeighborAddress(const Ptr<const NetworkHeaderBase>& networkHeader) const;

    // neighbor
    simtime_t getNextNeighborExpiration();
    void purgeNeighbors();
    std::vector<L3Address> getPlanarNeighbors() const;
    std::vector<L3Address> getPlanarNeighborsCounterClockwise(double startAngle) const;

    // next hop
    L3Address findNextHop(const L3Address& destination, GpsrOption *gpsrOption);
    L3Address findGreedyRoutingNextHop(const L3Address& destination, GpsrOption *gpsrOption);
    L3Address findPerimeterRoutingNextHop(const L3Address& destination, GpsrOption *gpsrOption);

    // routing
    Result routeDatagram(Packet *datagram, GpsrOption *gpsrOption);

    // netfilter
    virtual Result datagramPreRoutingHook(Packet *datagram) override;
    virtual Result datagramForwardHook(Packet *datagram) override { return ACCEPT; }
    virtual Result datagramPostRoutingHook(Packet *datagram) override { return ACCEPT; }
    virtual Result datagramLocalInHook(Packet *datagram) override { return ACCEPT; }
    virtual Result datagramLocalOutHook(Packet *datagram) override;

    void delayDatagram(Packet *datagram);
    void sendDelayedDatagram(Packet *datagram);
    bool hasDelayedDatagrams(const L3Address& target);
    void eraseDelayedDatagrams(const L3Address& target);

    // lifecycle
    virtual void handleStartOperation(LifecycleOperation *operation) override;
    virtual void handleStopOperation(LifecycleOperation *operation) override;
    virtual void handleCrashOperation(LifecycleOperation *operation) override;

    // notification
    virtual void receiveSignal(cComponent *source, simsignal_t signalID, cObject *obj, cObject *details) override;

};

} // namespace inet

#endif

