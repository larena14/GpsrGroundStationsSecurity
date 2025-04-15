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
#include "NonceGenerator.h"
#include "cpp-base64/base64.h"
#include <sodium.h>
#include <tuple>

#include <openssl/evp.h>
#include <openssl/ec.h>
#include <openssl/sha.h>
#include <openssl/err.h>
#include <openssl/rand.h>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <cstring>
#include <openssl/bio.h>
#include <openssl/buffer.h>

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

    enum SignAlgorithm {
        ECDSA = 0,
        EDDSA = 1
    };
  private:
    GpsrPlanarizationMode planarizationMode = static_cast<GpsrPlanarizationMode>(-1);
    const char *interfaces = nullptr;
    simtime_t beaconInterval;
    simtime_t maxJitter;
    simtime_t neighborValidityInterval;
    bool displayBubbles;

    cModule *host = nullptr;
    opp_component_ptr<IMobility> mobility;
    const IL3AddressType *addressType = nullptr;
    ModuleRefByPar<IInterfaceTable> interfaceTable;
    const char *outputInterface = nullptr;
    ModuleRefByPar<IRoutingTable> routingTable;
    ModuleRefByPar<INetfilter> networkProtocol;
    static PositionTable globalPositionTable;

    int positionByteLength = -1;

    cMessage *beaconTimer = nullptr;
    cMessage *purgeNeighborsTimer = nullptr;
    PositionTable neighborPositionTable;

    int n;
    int m;
    bool isMalicious;
    bool isStation;

    static double timeInPerimeter;
    static double timeInGreedy;
    static std::map<std::string, int> mapPacketHops;
    static std::map<std::string, simtime_t> mapPacketSendTime;
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

    static int msgNumber;
    static int totalMsg;
    static int totalArrivedMsg;
    static std::map<std::string, std::pair<simtime_t, int>> packetInfo;

    std::map<L3Address, std::tuple<Coord, std::string, uint64_t, simtime_t>> quadNodesSignatureMap;

    NonceGenerator nonceGen;

    EVP_PKEY* keyPairECDSA;
    static std::map<int, EVP_PKEY*> keyPairStationMapECDSA;
    static std::map<L3Address, EVP_PKEY*> keyPairDroneMapECDSA;

    std::vector<unsigned char> publicKeyEDDSA;
    std::vector<unsigned char> privateKeyEDDSA;
    static std::map<int, std::vector<unsigned char>> publicKeyStationMapEDDSA;
    static std::map<L3Address, std::vector<unsigned char>> publicKeyDroneMapEDDSA;

    static std::vector<std::string> maliciousStations;

    cMessage *sendTimerInterval = nullptr;

    SendMode sendMode = static_cast<SendMode>(-1);
    SignAlgorithm signAlgorithm = static_cast<SignAlgorithm>(-1);

    simtime_t minInterval;
    simtime_t maxInterval;

    static std::map<std::string, L3Address> possibleReceiversMap;

    static double verifyTime;
    static double signTime;

    static int verifyCount;
    static int signCount;

    static std::map<std::string, bool> finished;

    static double falsePositive;
    static double truePositive;
    static double falseNegative;
    static double trueNegative;

    cMessage *sendPositionRequestMessage = nullptr;

    static std::map<std::string, std::pair<double, double>> batteriesInfo;



  public:
    GpsrStation();
    virtual ~GpsrStation();
    static simsignal_t signSignal;
    static simsignal_t verifySignal;

  protected:
    virtual int numInitStages() const override { return NUM_INIT_STAGES; }
    void initialize(int stage) override;
    void handleMessageWhenUp(cMessage *message) override;
    void finish() override;

  private:
    double getBatteryLevel();

    int getStationIndex(Coord pos);

    void generateKeysECDSA();
    std::string signMessageECDSA(EVP_PKEY* pkey, std::string message, uint64_t nonce);
    bool verifySignatureECDSA(EVP_PKEY* pkey, std::string message, uint64_t nonce, std::string base64_signature);


    void printHex(const std::vector<unsigned char>& key);
    void generateKeysEDDSA();
    std::string signMessageEDDSA(std::string messageString, const std::vector<unsigned char> &privateKey, uint64_t nonce);
    std::string generateFakeSignatureEDDSA(std::string messageString, uint64_t nonce);
    bool verifySignatureEDDSA(std::string messageString, std::string signatureString, const std::vector<unsigned char> &publicKey, uint64_t nonce);

    const Ptr<StationNotice> createStationNotice(bool deregisterFlag);
    void sendStationNotice(const Ptr<StationNotice>& res, L3Address dest, const char* moduleName);
    void processStationNotice(Packet *packet);

    const Ptr<PositionRequest> createPositionRequest(L3Address address);
    void sendPositionRequest(const Ptr<PositionRequest>& res, L3Address dest, const char* moduleName);
    void processPositionRequest(Packet *packet);

    const Ptr<PositionResponse> createPositionResponse(L3Address address);
    const Ptr<PositionResponse> createPositionResponseFromAnotherStation(
            L3Address address, Coord position, simtime_t timestamp,
            std::string droneSignature, uint64_t droneNonce,
            std::string stationSignature, uint64_t stationNonce
    );
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

    void processSendPositionRequest();

    void processSelfMessage(cMessage *message);
    void processMessage(cMessage *message);

    void scheduleBeaconTimer();
    void processBeaconTimer();

    void schedulePurgeNeighborsTimer();
    void processPurgeNeighborsTimer();

    void sendUdpPacket(Packet *packet);
    void processUdpPacket(Packet *packet);

    const Ptr<GpsrBeacon> createBeacon();
    void sendBeacon(const Ptr<GpsrBeacon>& beacon);
    void processBeacon(Packet *packet);

    GpsrOption *createGpsrOption(L3Address destination);
    int computeOptionLength(GpsrOption *gpsrOption);
    void setGpsrOptionOnNetworkDatagram(Packet *packet, const Ptr<const NetworkHeaderBase>& networkHeader, GpsrOption *gpsrOption);

    GpsrOption *findGpsrOptionInNetworkDatagramForUpdate(const Ptr<NetworkHeaderBase>& networkHeader);
    const GpsrOption *findGpsrOptionInNetworkDatagram(const Ptr<const NetworkHeaderBase>& networkHeader) const;

    GpsrOption *getGpsrOptionFromNetworkDatagramForUpdate(const Ptr<NetworkHeaderBase>& networkHeader);
    const GpsrOption *getGpsrOptionFromNetworkDatagram(const Ptr<const NetworkHeaderBase>& networkHeader) const;

    void configureInterfaces();

    Coord lookupPositionInGlobalRegistry(const L3Address& address) const;
    void storePositionInGlobalRegistry(const L3Address& address, const Coord& position) const;
    void storeSelfPositionInGlobalRegistry() const;
    Coord computeIntersectionInsideLineSegments(Coord& begin1, Coord& end1, Coord& begin2, Coord& end2) const;
    Coord getNeighborPosition(const L3Address& address) const;
    Coord generateFalsePositionStation();
    double randomDouble(double min, double max);

    double getVectorAngle(Coord vector) const;
    double getNeighborAngle(const L3Address& address) const;

    std::string getHostName() const;
    L3Address getSelfAddress() const;
    L3Address getSenderNeighborAddress(const Ptr<const NetworkHeaderBase>& networkHeader) const;

    simtime_t getNextNeighborExpiration();
    void purgeNeighbors();
    std::vector<L3Address> getPlanarNeighbors() const;
    std::vector<L3Address> getPlanarNeighborsCounterClockwise(double startAngle) const;

    L3Address findNextHop(const L3Address& destination, GpsrOption *gpsrOption);
    L3Address findGreedyRoutingNextHop(const L3Address& destination, GpsrOption *gpsrOption);
    L3Address findPerimeterRoutingNextHop(const L3Address& destination, GpsrOption *gpsrOption);

    Result routeDatagram(Packet *datagram, GpsrOption *gpsrOption);

    virtual Result datagramPreRoutingHook(Packet *datagram) override;
    virtual Result datagramForwardHook(Packet *datagram) override { return ACCEPT; }
    virtual Result datagramPostRoutingHook(Packet *datagram) override { return ACCEPT; }
    virtual Result datagramLocalInHook(Packet *datagram) override { return ACCEPT; }
    virtual Result datagramLocalOutHook(Packet *datagram) override;

    void delayDatagram(Packet *datagram);
    void sendDelayedDatagram(Packet *datagram);
    bool hasDelayedDatagrams(const L3Address& target);
    void eraseDelayedDatagrams(const L3Address& target);

    virtual void handleStartOperation(LifecycleOperation *operation) override;
    virtual void handleStopOperation(LifecycleOperation *operation) override;
    virtual void handleCrashOperation(LifecycleOperation *operation) override;

    virtual void receiveSignal(cComponent *source, simsignal_t signalID, cObject *obj, cObject *details) override;

};

}

#endif

