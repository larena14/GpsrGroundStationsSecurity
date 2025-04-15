//
// Copyright (C) 2013 OpenSim Ltd.
//
// SPDX-License-Identifier: LGPL-3.0-or-later
//

#include "inet/physicallayer/wireless/common/contract/packetlevel/IRadio.h"
#include "inet/physicallayer/wireless/common/contract/packetlevel/RadioControlInfo_m.h"
#include "inet/common/INETDefs.h"
#include <inet/common/ModuleAccess.h>
#include <inet/linklayer/ieee80211/mac/Ieee80211Frame_m.h>
#include "inet/physicallayer/wireless/common/contract/packetlevel/SignalTag_m.h"
#include "inet/physicallayer/wireless/common/contract/packetlevel/IReception.h"
#include <numeric>
#include <cmath>
#include "stationGpsr/GpsrStation.h"
#include <algorithm>
#include <ctime>
#include <chrono>
#include <cstring> // Per strstr
#include <map>
#include <queue>
#include "nlohmann/json.hpp"
#include <iostream>
#include <fstream>
#include <iomanip>
#include <sstream>
#include <openssl/sha.h>
#include <string>
#include <random>
#include <openssl/evp.h>
#include <openssl/pem.h>
#include <openssl/ec.h>
#include <openssl/ecdsa.h>
#include <openssl/err.h>
#include <openssl/sha.h>
#include <openssl/err.h>
#include <vector>
#include <sstream>
#include <stdexcept>
#include <iostream>
#include <cstdlib>
#include <fstream>
#include <array>
#include <cctype>
#include <memory> // Per std::shared_ptr
#include <thread>
#include <random>


#include "inet/common/INETUtils.h"
#include "inet/common/IProtocolRegistrationListener.h"
#include "inet/common/ModuleAccess.h"
#include "inet/common/ProtocolTag_m.h"
#include "inet/common/TimeTag_m.h"
#include "inet/common/TimeTag.h"

#include "inet/common/lifecycle/ModuleOperations.h"
#include "inet/linklayer/common/InterfaceTag_m.h"
#include "inet/networklayer/common/HopLimitTag_m.h"
#include "inet/networklayer/common/IpProtocolId_m.h"
#include "inet/networklayer/common/L3AddressTag_m.h"
#include "inet/networklayer/common/L3Tools.h"
#include "inet/networklayer/common/NextHopAddressTag_m.h"
#include "inet/networklayer/contract/IInterfaceTable.h"
#include "messages/Station_m.h"

#include "inet/power/storage/SimpleEpEnergyStorage.h"
#include "inet/transportlayer/contract/udp/UdpSocket.h"



#ifdef INET_WITH_IPv4
#include "inet/networklayer/ipv4/Ipv4Header_m.h"
#endif

#ifdef INET_WITH_IPv6
#include "inet/networklayer/ipv6/Ipv6ExtensionHeaders_m.h"
#include "inet/networklayer/ipv6/Ipv6InterfaceData.h"
#endif

#ifdef INET_WITH_NEXTHOP
#include "inet/networklayer/nexthop/NextHopForwardingHeader_m.h"
#endif

using json = nlohmann::json;

namespace inet {

Define_Module(GpsrStation);

simsignal_t GpsrStation::signSignal = cComponent::registerSignal("signSignal");
simsignal_t GpsrStation::verifySignal = cComponent::registerSignal("verifySignal");

std::map<std::string, std::pair<double, double>> GpsrStation::batteriesInfo;
std::map<std::string, bool>  GpsrStation::finished;

double GpsrStation::timeInPerimeter;
double GpsrStation::timeInGreedy;
bool GpsrStation::printed;
std::map<std::string, int> GpsrStation::mapPacketHops;
std::map<int, L3Address> GpsrStation::stationMap;
std::map<std::string, L3Address> GpsrStation::hostnameMap;

std::vector<L3Address> GpsrStation::maliciousNode;
std::set<L3Address> GpsrStation::maliciousIndividuati;
std::set<L3Address> GpsrStation::maliciousTest;
std::set<L3Address> GpsrStation::falsiPositivi;

int GpsrStation::msgNumber;
int GpsrStation::totalMsg;
int GpsrStation::totalArrivedMsg;

double GpsrStation::verifyTime;
double GpsrStation::signTime;
int GpsrStation::verifyCount;
int GpsrStation::signCount;

double GpsrStation::falsePositive;
double GpsrStation::truePositive;
double GpsrStation::falseNegative;
double GpsrStation::trueNegative;

double GpsrStation::totalDifference;
std::string GpsrStation:: alg;


std::map<std::string, std::pair<simtime_t, int>> GpsrStation::packetInfo;
std::map<std::string,double> GpsrStation::timeMessages;

std::map<int, std::vector<unsigned char>> GpsrStation::publicKeyMapEDDSA;
std::map<int, EVP_PKEY*> GpsrStation::keyPairStationMapECDSA;

std::map<std::string, L3Address> GpsrStation::possibleReceiversMap;
static inline double determinant(double a1, double a2, double b1, double b2)
{
    return a1 * b2 - a2 * b1;
}

GpsrStation::GpsrStation()
{

}

GpsrStation::~GpsrStation()
{
    cancelAndDelete(beaconTimer);
    cancelAndDelete(purgeNeighborsTimer);
    cancelAndDelete(purgeStationTimer);
    cancelAndDelete(purgeDestinationTimer);
    cancelAndDelete(sendTimerInterval);
}

//
// module interface
//

void GpsrStation::initialize(int stage)
{
    if (stage == INITSTAGE_ROUTING_PROTOCOLS)
        addressType = getSelfAddress().getAddressType();

    RoutingProtocolBase::initialize(stage);

    if (stage == INITSTAGE_LOCAL) {
        // Gpsr parameters
        const char *planarizationModeString = par("planarizationMode");
        if (!strcmp(planarizationModeString, ""))
            planarizationMode = GPSR_NO_PLANARIZATION;
        else if (!strcmp(planarizationModeString, "GG"))
            planarizationMode = GPSR_GG_PLANARIZATION;
        else if (!strcmp(planarizationModeString, "RNG"))
            planarizationMode = GPSR_RNG_PLANARIZATION;
        else
            throw cRuntimeError("Unknown planarization mode");
        interfaces = par("interfaces");
        beaconInterval = par("beaconInterval");
        maxJitter = par("maxJitter");
        neighborValidityInterval = par("neighborValidityInterval");
        displayBubbles = par("displayBubbles");


        // context
        host = getContainingNode(this);
        interfaceTable.reference(this, "interfaceTableModule", true);
        outputInterface = par("outputInterface");
        mobility = check_and_cast<IMobility *>(host->getSubmodule("mobility"));
        routingTable.reference(this, "routingTableModule", true);
        networkProtocol.reference(this, "networkProtocolModule", true);
        // internal
        beaconTimer = new cMessage("BeaconTimer");
        purgeNeighborsTimer = new cMessage("PurgeNeighborsTimer");
        // packet size
        positionByteLength = par("positionByteLength");
        // KLUDGE implement position registry protocol
        globalPositionTable.clear();

        n = par("n");
        m = par("m");
        isMalicious = par("isMalicious");
        isStation = par("isStation");

        timeInPerimeter = 0.0;
        timeInGreedy = 0.0;
        printed = false;
        mapPacketHops.clear();
        stationMap.clear();
        hostnameMap.clear();

        purgeStationTimer = new cMessage("PurgeStationTimer");
        positionValidityInterval = par("positionValidityInterval");

        purgeDestinationTimer = new cMessage("PurgeDestinationTimer");
        destinationValidityInterval = par("destinationValidityInterval");

        const char *sendModeString = par("sendMode");
                        if (!strcmp(sendModeString, "FixedMaximum")){
                            sendMode = FIXED_MAXIMUM;
                            alg="FixedMaximum";
                        }
                        else if (!strcmp(sendModeString, "Interval")){
                            sendMode = INTERVAL;
                            alg="Interval";
                        }
                        else
                            throw cRuntimeError("Unknown send mode");

                        const char *signAlgorithmString = par("signAlgorithm");
                        if (!strcmp(signAlgorithmString, "ECDSA")){
                            signAlgorithm = ECDSA;
                        }
                        else if (!strcmp(signAlgorithmString, "EDDSA")){
                            signAlgorithm = EDDSA;
                        }
                        else
                            throw cRuntimeError("Unknown sign algorithm");

                 minInterval=par("minInterval");
                 maxInterval=par("maxInterval");

                 if (sendMode == FIXED_MAXIMUM){
                                      std::vector<std::string> receiversStr = cStringTokenizer(par("receivers")).asVector();

                                      for (auto i : receiversStr) {
                                          receivers.push(i);
                                      }

                                      std::vector<std::string> timesStr = cStringTokenizer(par("sendTimes")).asVector();
                                      for (auto i : timesStr) {

                                          cMessage* sendTimer = new cMessage("SendTimer");
                                          simtime_t t = simtime_t::parse(i.c_str());
                                          scheduleAt(t, sendTimer);

                                      }
                                  }
                                  else
                                      sendTimerInterval = new cMessage("SendTimerInterval");

                 msgNumber = 0;
                         totalMsg = 0;
                         totalArrivedMsg = 0;

                         verifyTime = 0.0;
                         signTime = 0.0;
                         verifyCount = 0;
                         signCount = 0;
                         totalDifference=0.0;

                         falsePositive = 0.0;
                         truePositive = 0.0;
                         falseNegative = 0.0;
                         trueNegative = 0.0;

                         timeMessages.clear();
                         packetInfo.clear();
                         trustness.clear();

                         publicKeyMapEDDSA.clear();
                         keyPairStationMapECDSA.clear();

                         possibleReceiversMap.clear();

                         batteriesInfo.clear();
                         finished.clear();




    }
    else if (stage == INITSTAGE_ROUTING_PROTOCOLS) {
        registerProtocol(Protocol::manet, gate("ipOut"), gate("ipIn"));
        host->subscribe(linkBrokenSignal, this);
        networkProtocol->registerHook(0, this);
        WATCH(neighborPositionTable);
        WATCH(quadNodesPositionTable);
        WATCH(destinationsPositionTable);
    }
}

void GpsrStation::finish(){
    if(!isStation) batteriesInfo[host->getFullName()].second = getBatteryLevel();
           finished[host->getFullName()] = true;
           bool allTrue = true;
           for (const auto& pair : finished) {
               if (!pair.second) {
                   allTrue = false;
                   break;
               }
           }




    if(allTrue){

        for (const auto& pair : timeMessages) {
                            totalDifference += pair.second;
                }

                json data;
                double batteryRes;
                data["Algorithm"]=alg;
                data["PerimeterTime"] = timeInPerimeter;
                data["GreedyTime"] = timeInGreedy;
                data["totalMsg"] = totalMsg;
                data["totalArrivedMsg"] = totalArrivedMsg;
                data["Malicious Found"] = maliciousIndividuati.size();
                data["False positive"] = falsiPositivi.size();
                data["meanMessagesTime"]=totalDifference/totalArrivedMsg;
                data["falsePositive"] = falsePositive;
                data["truePositive"] = truePositive;
                data["falseNegative"] = falseNegative;
                data["trueNegative"] = trueNegative;

                if (signCount > 0 && verifyCount > 0) {
                                   data["meanSignTime"] = (signTime/signCount);
                                   data["meanVerifyTime"] = (verifyTime/verifyCount);
                               }

                json messageArray;
                for(auto& message: packetInfo){
                    for(const auto pair:timeMessages){
                        if(message.first==pair.first){
                    json info = {
                        {"msgName", message.first},
                        {"time", message.second.first.str()},
                        {"hops", message.second.second},
                        {"receptionTime",pair.second}
                    };

                    messageArray.push_back(info);}
               }
                }
        data["messages"] = messageArray;
        json batteryArray;
                        for(auto& b: batteriesInfo){
                            json info = {
                                {"hostname", b.first},
                                {"J consumed", b.second.first - b.second.second}
                            };
                            batteryArray.push_back(info);
                        }
                        data["batteries"] = batteryArray;

        json droppedMessageArray;

        json notArrivedArray;
        for (int i = 0; i < totalMsg; ++i) {
             std::string msgName = "SimpleMessage" + std::to_string(i);
             if (packetInfo.find(msgName) == packetInfo.end()){
                   json info = {
                        {"msgName", msgName},
                        {"hops", mapPacketHops[msgName]}
                   };
             notArrivedArray.push_back(info);
             }
        }
        data["notArrivedMessages"] = notArrivedArray;

        std::ofstream file("results/json/output.json");
        if (file.is_open()) {
            file << data.dump(4);
            file.close();
            std::cout << "File JSON salvato con successo in 'output.json'." << std::endl;
        } else {
            std::cerr << "Errore: impossibile salvare il file JSON." << std::endl;
        }

        printed = true;
    }
}

void GpsrStation::handleMessageWhenUp(cMessage *message)
{
    if (message->isSelfMessage())
        processSelfMessage(message);
    else
        processMessage(message);
}

double GpsrStation::getBatteryLevel() {
    // Trova il modulo di energia del drone
    auto energyStorage = getModuleByPath("^.energyStorage");

    if (energyStorage) {
        // Converte in un puntatore alla classe corretta
        auto battery = dynamic_cast<inet::power::SimpleEpEnergyStorage*>(energyStorage);

        if (battery) {
            return battery->getResidualEnergyCapacity().get();
        }
    }
    return 0;
}



int GpsrStation::getStationIndex(Coord pos){
    double l = n / m;

    int qX = static_cast<int>(std::floor(pos.getX() / l));
    int qY = static_cast<int>(std::floor(pos.getY() / l));

    if (qX >= m)
        qX = m - 1;
    if (qY >= m)
        qY = m - 1;

    int q = qY * m + qX;

    return q;
}

Coord GpsrStation::generateFalsePosition(){
    double prob = uniform(0.0, 1.0);
            if (prob <= 0.33){

                /* Opposta nella mappa */

                std::vector<Coord> angles;
                Coord ul = Coord();
                Coord ur = Coord();
                Coord dl = Coord();
                Coord dr = Coord();

                ul.setX(0);
                ul.setY(0);
                ul.setZ(0);

                ur.setX(n);
                ur.setY(0);
                ur.setZ(0);

                dl.setX(0);
                dl.setY(n);
                dl.setZ(0);

                dr.setX(n);
                dr.setY(n);
                dr.setZ(0);

                angles.push_back(ul);
                angles.push_back(ur);
                angles.push_back(dl);
                angles.push_back(dr);

                Coord pos = mobility->getCurrentPosition();

                double maxDistance = -1;
                int maxAngleIndex = 0;
                for (size_t i = 0; i < angles.size(); ++i) {
                    double currDist = pos.distance(angles[i]);
                    if (currDist > maxDistance){
                        maxDistance = currDist;
                        maxAngleIndex = i;
                    }
                }

                Coord f = Coord();
                Coord maxAngle = angles[maxAngleIndex];
                double delta = n / m;
                if (maxAngle == ul){
                    /* UP LEFT */
                    f.setX(uniform(0, delta));
                    f.setY(uniform(0, delta));
                    f.setZ(uniform(100, 120));
                }
                else if (maxAngle == ur){
                    /* UP RIGHT */
                    f.setX(uniform(n - delta, n));
                    f.setY(uniform(0, delta));
                    f.setZ(uniform(100, 120));
                }
                else if (maxAngle == dl){
                    /* DOWN LEFT */
                    f.setX(uniform(0, delta));
                    f.setY(uniform(n - delta, n));
                    f.setZ(uniform(100, 120));
                }
                else{
                    /* DOWN RIGHT */
                    f.setX(uniform(n - delta, n));
                    f.setY(uniform(n - delta, n));
                    f.setZ(uniform(100, 120));
                }


                return f;
            }
            else if (prob > 0.33 && prob <= 0.66){
                /* Opposta nel quadrante */

                Coord pos = mobility->getCurrentPosition();
                double x = pos.getX();
                double y = pos.getY();

                double l = n / m;

                int qX = std::floor(x / l);
                int qY = std::floor(y / l);

                if (qX >= m) qX = m - 1;
                if (qY >= m) qY = m - 1;

                double TL_x = qX * l;
                double TL_y = qY * l;
                double TR_x = (qX + 1) * l;
                double TR_y = qY * l;
                double BL_x = qX * l;
                double BL_y = (qY + 1) * l;
                double BR_x = (qX + 1) * l;
                double BR_y = (qY + 1) * l;

                double center_x = (TL_x + BR_x) / 2;
                double center_y = (TL_y + BR_y) / 2;

                double opposite_x = 2 * center_x - x;
                double opposite_y = 2 * center_y - y;

                Coord falsePosition = Coord();

                falsePosition.setX(opposite_x);
                falsePosition.setY(opposite_y);
                falsePosition.setZ(uniform(100, 120));


                return falsePosition;
            }
            else{
                /* Random */

                Coord f = Coord();

                f.setX(uniform(0, n));
                f.setY(uniform(0, n));
                f.setZ(uniform(100, 120));



                return f;
            }
}


void GpsrStation::generateKeysEDDSA() {
    publicKeyEDDSA.resize(crypto_sign_PUBLICKEYBYTES);
    privateKeyEDDSA.resize(crypto_sign_SECRETKEYBYTES);

    if (crypto_sign_keypair(publicKeyEDDSA.data(), privateKeyEDDSA.data()) != 0) {
        throw std::runtime_error("Errore nella generazione delle chiavi");
    }
}

std::string GpsrStation::signMessageEDDSA(std::string messageString, const std::vector<unsigned char> &privateKey, uint64_t nonce) {
    auto start = std::chrono::high_resolution_clock::now();
    std::vector<unsigned char> message(messageString.begin(), messageString.end());
    std::vector<unsigned char> messageWithNonce(sizeof(nonce) + message.size());
    std::memcpy(messageWithNonce.data(), &nonce, sizeof(nonce));
    std::memcpy(messageWithNonce.data() + sizeof(nonce), message.data(), message.size());

    std::vector<unsigned char> signature(crypto_sign_BYTES + messageWithNonce.size());
    unsigned long long signatureLen;

    if (crypto_sign(signature.data(), &signatureLen, messageWithNonce.data(), messageWithNonce.size(), privateKey.data()) != 0) {
        throw std::runtime_error("Errore nella firma del messaggio");
    }
    auto end = std::chrono::high_resolution_clock::now();
         std::chrono::duration<double> duration = end - start;
         signTime += duration.count();
         signCount += 1;


    signature.resize(signatureLen);
    emit(signSignal, (double) message.size());
    return base64_encode(signature.data(), signature.size());
}

bool GpsrStation::verifySignatureEDDSA(std::string messageString, std::string signatureString, const std::vector<unsigned char> &publicKey, uint64_t nonce){
    auto start = std::chrono::high_resolution_clock::now();
    std::vector<unsigned char> message(messageString.begin(), messageString.end());
    std::string decodedSignature = base64_decode(signatureString);
    std::vector<unsigned char> signature(decodedSignature.begin(), decodedSignature.end());

    std::vector<unsigned char> messageWithNonce(sizeof(nonce) + message.size());
    std::memcpy(messageWithNonce.data(), &nonce, sizeof(nonce));
    std::memcpy(messageWithNonce.data() + sizeof(nonce), message.data(), message.size());

    std::vector<unsigned char> verifiedMessage(signature.size() + crypto_sign_BYTES);
    unsigned long long verifiedMessageLen;

    if (crypto_sign_open(verifiedMessage.data(), &verifiedMessageLen, signature.data(), signature.size(), publicKey.data()) != 0) {
        return false;
    }

    verifiedMessage.resize(verifiedMessageLen);
    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> duration = end - start;
    verifyTime += duration.count();
    verifyCount += 1;

    emit(verifySignal, (double) message.size());
    return std::equal(verifiedMessage.begin(), verifiedMessage.end(), messageWithNonce.begin());
}




float GpsrStation::distanza(Coord punto1, Coord punto2) {
    // Calcola la distanza Euclidea tra i due punti
    float distanza = sqrt(pow(punto2.x - punto1.x, 2) + pow(punto2.y - punto1.y, 2) + pow(punto2.z - punto1.z, 2));
    return distanza;
}

//ECDSA
void GpsrStation::generateKeysECDSA() {
    EVP_PKEY_CTX* ctx = EVP_PKEY_CTX_new_id(EVP_PKEY_EC, nullptr);
    if (!ctx) throw std::runtime_error("Errore nella generazione delle chiavi ECDSA");

    if (EVP_PKEY_keygen_init(ctx) <= 0) throw std::runtime_error("Errore nella generazione delle chiavi ECDSA");
    if (EVP_PKEY_CTX_set_ec_paramgen_curve_nid(ctx, NID_secp256k1) <= 0) throw std::runtime_error("Errore nella generazione delle chiavi ECDSA");

    keyPairECDSA = nullptr;
    if (EVP_PKEY_keygen(ctx, &keyPairECDSA) <= 0) throw std::runtime_error("Errore nella generazione delle chiavi ECDSA");

    EVP_PKEY_CTX_free(ctx);
}

std::string GpsrStation::signMessageECDSA(EVP_PKEY* pkey, std::string message, uint64_t nonce) {
    auto start = std::chrono::high_resolution_clock::now();

    EVP_MD_CTX* mdctx = EVP_MD_CTX_new();
    if (!mdctx) throw std::runtime_error("Errore nella firma ECDSA");

    if (EVP_DigestSignInit(mdctx, nullptr, EVP_sha256(), nullptr, pkey) <= 0) throw std::runtime_error("Errore nella firma ECDSA");

    std::ostringstream full_message;
    full_message << message << nonce;
    std::string final_message = full_message.str();

    if (EVP_DigestSignUpdate(mdctx, final_message.c_str(), final_message.size()) <= 0) throw std::runtime_error("Errore nella firma ECDSA");

    size_t siglen = 0;
    if (EVP_DigestSignFinal(mdctx, nullptr, &siglen) <= 0) throw std::runtime_error("Errore nella firma ECDSA");

    std::vector<unsigned char> signature(siglen);
    if (EVP_DigestSignFinal(mdctx, signature.data(), &siglen) <= 0) throw std::runtime_error("Errore nella firma ECDSA");

    EVP_MD_CTX_free(mdctx);

    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> duration = end - start;
    signTime += duration.count();
    signCount += 1;
    signature.resize(siglen);
    emit(signSignal, (double) message.size());
    return base64_encode(signature.data(), signature.size());
}

bool GpsrStation::verifySignatureECDSA(EVP_PKEY* pkey, std::string message, uint64_t nonce, std::string base64_signature) {
    auto start = std::chrono::high_resolution_clock::now();

    std::string decodedSignature = base64_decode(base64_signature);
    std::vector<unsigned char> signature(decodedSignature.begin(), decodedSignature.end());

    EVP_MD_CTX* mdctx = EVP_MD_CTX_new();
    if (!mdctx) {
        return false; // throw std::runtime_error("Errore nella verifica ECDSA - 1");
    }

    if (EVP_DigestVerifyInit(mdctx, nullptr, EVP_sha256(), nullptr, pkey) <= 0){
        return false; // throw std::runtime_error("Errore nella verifica ECDSA - 2");
    }

    std::ostringstream full_message;
    full_message << message << nonce;
    std::string final_message = full_message.str();

    if (EVP_DigestVerifyUpdate(mdctx, final_message.c_str(), final_message.size()) <= 0){
        return false; // throw std::runtime_error("Errore nella verifica ECDSA - 3");
    }

    int result = EVP_DigestVerifyFinal(mdctx, signature.data(), signature.size());

    EVP_MD_CTX_free(mdctx);
    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> duration = end - start;
    verifyTime += duration.count();
    verifyCount += 1;
    emit(verifySignal, (double) message.size());
    return (result == 1);
}






//FINE FUNZIONI CRITTOGRAFICHE

const Ptr<StationNotice> GpsrStation::createStationNotice(bool deregisterFlag)
{
    const auto& stationNotice = makeShared<StationNotice>();
    stationNotice->setSource(getSelfAddress());
    stationNotice->setSourceModuleName(host->getName());
    std::vector<L3Address> v= neighborPositionTable.getAddresses();
    size_t size=v.size();
    stationNotice->setAddressListArraySize(size);
    stationNotice->setDistanceListArraySize(size);
    for(int i=0;i<size;i++){
        /*
         * inviamo alla stazione la lista dei vicini con le
         * relative distanze che abbiamo stimato
         *
         * */
        L3Address element=v[i];
        stationNotice->setAddressList(i,element);
        stationNotice->setDistanceList(i,mappaDistanzaIndirizzo[element]);
    }
    if(isMalicious){ // && destinatioSetted
        stationNotice->setPosition(generateFalsePosition());
    }else{
        stationNotice->setPosition(mobility->getCurrentPosition());
    }
    stationNotice->setDeregister(deregisterFlag);
    //uint64_t nonce = nonceGen.generateNonce();
    //std::string nonce = generateRandomString(16);
    stationNotice->setChunkLength(B(getSelfAddress().getAddressType()->getAddressByteLength() + positionByteLength + 64 + 1));
    return stationNotice;
}

void GpsrStation::sendStationNotice(const Ptr<StationNotice>& req, L3Address dest, const char* moduleName)
{
    Packet *udpPacket = new Packet("StationNotice");
    udpPacket->insertAtBack(req);
    auto udpHeader = makeShared<UdpHeader>();
          udpHeader->setSourcePort(GPSR_UDP_PORT);
          udpHeader->setDestinationPort(GPSR_UDP_PORT);
          udpHeader->setCrcMode(CRC_DISABLED);
          udpPacket->insertAtFront(udpHeader);
          auto addresses = udpPacket->addTag<L3AddressReq>();
          addresses->setSrcAddress(getSelfAddress());
          addresses->setDestAddress(dest);
          udpPacket->addTag<HopLimitReq>()->setHopLimit(255);
          udpPacket->addTag<PacketProtocolTag>()->setProtocol(&Protocol::manet);
          udpPacket->addTag<DispatchProtocolReq>()->setProtocol(addressType->getNetworkProtocol());
         // std::cout << host->getFullName() << " send StationNotice to " << moduleName << " with address: " << dest << endl;
          sendUdpPacket(udpPacket);


}

void GpsrStation::checkDistanceConsistency(L3Address nodeAddress, Coord reportedPosition, double PENALTY) {
    // Definiamo il peso per il nuovo valore (ad esempio 0.1)
    double alpha = 0.3;

    std::map<L3Address,double> mapA;
    for (L3Address indirizzo : neighboursMap[nodeAddress]) { // i vicini del nodo che fa la richiesta
        if (neighboursMapDistance.find(indirizzo) != neighboursMapDistance.end() && trustness[indirizzo] > notTrusted) {
            mapA = neighboursMapDistance[indirizzo]; // mappa degli indirizzi
            double distance = mapA[nodeAddress]; // distanza tra 'indirizzo' e 'nodeAddress'
            if (distance == -1) {
                // Misurazione negativa: penalizziamo nodeAddress
                double measurement = trustness[nodeAddress] - PENALTY;
                trustness[nodeAddress] = (1 - alpha) * trustness[nodeAddress] + alpha * measurement;
            } else {
                // Misurazione positiva: premiamo 'indirizzo'
                double measurement = trustness[indirizzo] + PENALTY;
                trustness[indirizzo] = (1 - alpha) * trustness[indirizzo] + alpha * measurement;
            }
        }
    }
}


void GpsrStation::processStationNotice(Packet *packet){
    packet->popAtFront<UdpHeader>();
    bool trovato=false;
    const auto& req = packet->peekAtFront<StationNotice>();
    L3Address address = req->getSource();
    Coord position = req->getPosition();
    if(req->getDeregister()){
        if(quadNodesPositionTable.hasPosition(address)){
            quadNodesPositionTable.removePosition(address);
        }
    }
    else{

        std::list<L3Address>list;
        std::map<L3Address,double> mappingIndirizzoDistanza;

        size_t size= req->getAddressListArraySize();
        for(int i=0;i<size;i++){
                     list.push_back(req->getAddressList(i)); //unelemento della lista
                     mappingIndirizzoDistanza.insert({req->getAddressList(i),req->getDistanceList(i)});
                   }
        neighboursMap[address]=list;
        neighboursMapDistance[address]=mappingIndirizzoDistanza;

        checkDistanceConsistency(address,position,1);

        if(trustness[address]>notTrusted){
            quadNodesPositionTable.setPosition(address, position);
            sendStationNoticeResponse(createStationNoticeResponse(position.getX(), position.getY(), position.getZ()),req->getSource(), req->getSourceModuleName());
            }
            else{
                if(quadNodesPositionTable.hasPosition(address)){
                    quadNodesPositionTable.removePosition(address);}

                int c=count(maliciousNode.begin(),maliciousNode.end(),address);
                if(c>0){
                    maliciousIndividuati.insert(address);
                    trueNegative += 1.0;
                }
                else{
                    falsiPositivi.insert(address);
                    falseNegative += 1.0; //false negative
                }
                neighboursMap.erase(address);
            }
        if(trustness[address]>trusted){
            int c=count(maliciousNode.begin(),maliciousNode.end(),address);
                            if(c>0){
                                falsePositive += 1.0;
                            }
                            else{
                                truePositive += 1.0;
                            }
        }


    }

    delete packet;
}



const Ptr<StationNoticeResponse> GpsrStation::createStationNoticeResponse(double x, double y, double z){

    const auto& stationResponse = makeShared<StationNoticeResponse>();
         std::string clear = std::to_string(x) + "-" + std::to_string(y) + "-" + std::to_string(z);
         uint64_t nonce = nonceGen.generateNonce();
         std::string signature;
         if (signAlgorithm == ECDSA){
             signature=signMessageECDSA(keyPairECDSA, clear, nonce);
         }
         else{
             signature = signMessageEDDSA(clear, privateKeyEDDSA, nonce);
         }
         stationResponse->setC(signature.c_str());
         stationResponse->setX(x);
         stationResponse->setY(y);
         stationResponse->setZ(z);
         stationResponse->setNonce(nonce);
         stationResponse->setChunkLength(B(signature.length() + 16*3 + sizeof(nonce) + positionByteLength));
         return stationResponse;

}
void GpsrStation::sendStationNoticeResponse(const Ptr<StationNoticeResponse>& res, L3Address dest, const char* moduleName){

    Packet *udpPacket = new Packet("StationNoticeResponse");
    udpPacket->insertAtBack(res);
    auto udpHeader = makeShared<UdpHeader>();
          udpHeader->setSourcePort(GPSR_UDP_PORT);
          udpHeader->setDestinationPort(GPSR_UDP_PORT);
          udpHeader->setCrcMode(CRC_DISABLED);
          udpPacket->insertAtFront(udpHeader);
          auto addresses = udpPacket->addTag<L3AddressReq>();
          addresses->setSrcAddress(getSelfAddress());
          addresses->setDestAddress(dest);
          udpPacket->addTag<HopLimitReq>()->setHopLimit(255);
          udpPacket->addTag<PacketProtocolTag>()->setProtocol(&Protocol::manet);
          udpPacket->addTag<DispatchProtocolReq>()->setProtocol(addressType->getNetworkProtocol());
          sendUdpPacket(udpPacket);

}

void GpsrStation::processStationNoticeResponse(Packet *packet){

    packet->popAtFront<UdpHeader>();
    const auto& res = packet->peekAtFront<StationNoticeResponse>();
    std::string clear = std::to_string(res->getX()) + "-" + std::to_string(res->getY()) + "-" + std::to_string(res->getZ());
    Coord pos = Coord();
    pos.setX(res->getX());
    pos.setY(res->getY());
    pos.setZ(res->getZ());
    sendBeacon(createBeacon(pos, res->getC(), res->getNonce()));
    storeSelfPositionInGlobalRegistry();
    delete packet;

}

const Ptr<PositionRequest> GpsrStation::createPositionRequest(L3Address address)
{
    const auto& positionRequest = makeShared<PositionRequest>();
    positionRequest->setSource(getSelfAddress());
    const char* sourceModuleName = host->getFullName();
    positionRequest->setSourceModuleName(sourceModuleName);
    positionRequest->setAddress(address);
    positionRequest->setChunkLength(B(getSelfAddress().getAddressType()->getAddressByteLength() + std::strlen(sourceModuleName) + address.getAddressType()->getAddressByteLength()));
    return positionRequest;
}

void GpsrStation::sendPositionRequest(const Ptr<PositionRequest>& req, L3Address dest, const char* moduleName)
{

    Packet *udpPacket = new Packet("PositionRequest");
    udpPacket->insertAtBack(req);
    auto udpHeader = makeShared<UdpHeader>();
          udpHeader->setSourcePort(GPSR_UDP_PORT);
          udpHeader->setDestinationPort(GPSR_UDP_PORT);
          udpHeader->setCrcMode(CRC_DISABLED);
          udpPacket->insertAtFront(udpHeader);
          auto addresses = udpPacket->addTag<L3AddressReq>();
          addresses->setSrcAddress(getSelfAddress());
          addresses->setDestAddress(dest);
          udpPacket->addTag<HopLimitReq>()->setHopLimit(255);
          udpPacket->addTag<PacketProtocolTag>()->setProtocol(&Protocol::manet);
          udpPacket->addTag<DispatchProtocolReq>()->setProtocol(addressType->getNetworkProtocol());
          sendUdpPacket(udpPacket);
}

void GpsrStation::processPositionRequest(Packet *packet){
    packet->popAtFront<UdpHeader>();
    const auto& req = packet->peekAtFront<PositionRequest>();
    L3Address source = req->getSource();
    const char* moduleName = req->getSourceModuleName();
    L3Address address = req->getAddress(); //indirzzo della destinazione address


    if (quadNodesPositionTable.hasPosition(address)){

        sendPositionResponse(createPositionResponse(address), source, moduleName);
    }
    else{
        for (const auto& pair : stationMap) {
            std::string stationModule = "station" + std::to_string(pair.first);
            L3Address stationAddress = pair.second;
            sendS2SPositionRequest(createS2SPositionRequest(source, moduleName, address), stationAddress, stationModule.c_str());
        }
    }
    delete packet;
}

const Ptr<PositionResponse> GpsrStation::createPositionResponse(L3Address address)
{
    const auto& positionResponse = makeShared<PositionResponse>();
    // doppio controllo, la posizione potrebbe essere stata eiminata perchè troppo vecchia
    if (quadNodesPositionTable.hasPosition(address)){
        positionResponse->setSetted(true);
        positionResponse->setAddress(address);
        PositionTable::PositionWithTimestamp p = quadNodesPositionTable.getPositionWithTimestamp(address);
        positionResponse->setPosition(p.position);
        positionResponse->setTime(p.timestamp);
    }
    else{
        positionResponse->setSetted(false);
    }
    positionResponse->setChunkLength(B(1 + address.getAddressType()->getAddressByteLength() + positionByteLength + sizeof(positionResponse->getTime())));
    return positionResponse;
}

void GpsrStation::sendPositionResponse(const Ptr<PositionResponse>& req, L3Address dest, const char* moduleName)
{
    Packet *udpPacket = new Packet("PositionResponse");
    udpPacket->insertAtBack(req);
    auto udpHeader = makeShared<UdpHeader>();
          udpHeader->setSourcePort(GPSR_UDP_PORT);
          udpHeader->setDestinationPort(GPSR_UDP_PORT);
          udpHeader->setCrcMode(CRC_DISABLED);
          udpPacket->insertAtFront(udpHeader);
          auto addresses = udpPacket->addTag<L3AddressReq>();
          addresses->setSrcAddress(getSelfAddress());
          addresses->setDestAddress(dest);
          udpPacket->addTag<HopLimitReq>()->setHopLimit(255);
          udpPacket->addTag<PacketProtocolTag>()->setProtocol(&Protocol::manet);
          udpPacket->addTag<DispatchProtocolReq>()->setProtocol(addressType->getNetworkProtocol());
          sendUdpPacket(udpPacket);
}

void GpsrStation::processPositionResponse(Packet *packet){
    packet->popAtFront<UdpHeader>();
    const auto& req = packet->peekAtFront<PositionResponse>();
    bool setted = req->getSetted();
    if(setted){
        L3Address address = req->getAddress();
        Coord position = req->getPosition();
        simtime_t time = req->getTime();
        destinationsPositionTable.setPosition(address, position);
        if (hasDelayedDatagrams(address)){
            auto lt = targetAddressToDelayedPackets.lower_bound(address);
            auto ut = targetAddressToDelayedPackets.upper_bound(address);
            for (auto it = lt; it != ut; it++){
                sendDelayedDatagram(it->second);
            }
            eraseDelayedDatagrams(address);
        }

    }
    delete packet;
}


const Ptr<S2SPositionRequest> GpsrStation::createS2SPositionRequest(L3Address applicant, const char* applicantModuleName, L3Address address)
{
    const auto& positionRequest = makeShared<S2SPositionRequest>();
    positionRequest->setSource(getSelfAddress());
    const char* sourceModuleName = host->getFullName();
    positionRequest->setSourceModuleName(sourceModuleName);
    positionRequest->setApplicant(applicant);
    positionRequest->setApplicantModuleName(applicantModuleName);
    positionRequest->setAddress(address);
    positionRequest->setChunkLength(B(
            getSelfAddress().getAddressType()->getAddressByteLength()
            + std::strlen(sourceModuleName)
            + applicant.getAddressType()->getAddressByteLength()
            + std::strlen(applicantModuleName)
            + address.getAddressType()->getAddressByteLength()));
    return positionRequest;
}

void GpsrStation::sendS2SPositionRequest(const Ptr<S2SPositionRequest>& req, L3Address dest, const char* moduleName)
{
    Packet *udpPacket = new Packet("S2SPositionRequest");
    udpPacket->insertAtBack(req);
    auto udpHeader = makeShared<UdpHeader>();
    udpHeader->setSourcePort(GPSR_UDP_PORT);
    udpHeader->setDestinationPort(GPSR_UDP_PORT);
    udpHeader->setCrcMode(CRC_DISABLED);
    udpPacket->insertAtFront(udpHeader);
    auto addresses = udpPacket->addTag<L3AddressReq>();
    addresses->setSrcAddress(getSelfAddress());
    addresses->setDestAddress(dest);
    udpPacket->addTag<HopLimitReq>()->setHopLimit(255);
    udpPacket->addTag<PacketProtocolTag>()->setProtocol(&Protocol::udp);
    udpPacket->addTag<DispatchProtocolReq>()->setProtocol(&Protocol::ipv4);
    udpPacket->addTag<InterfaceReq>()->setInterfaceId(interfaceTable->findInterfaceByName("wlan0")->getInterfaceId());
    cModule *destination = getParentModule()->getParentModule()->getSubmodule(moduleName);
    cModule *stationGpsr = destination->getSubmodule("gpsr");
    sendDirect(udpPacket, stationGpsr, "directIn");
}

void GpsrStation::processS2SPositionRequest(Packet *packet){
    packet->popAtFront<UdpHeader>();
    const auto& req = packet->peekAtFront<S2SPositionRequest>();
    L3Address applicant = req->getApplicant();
    const char* applicantModuleName = req->getApplicantModuleName();
    L3Address source = req->getSource();
    const char* moduleName = req->getSourceModuleName();
    L3Address address = req->getAddress();
    if (quadNodesPositionTable.hasPosition(address)){
        sendS2SPositionResponse(createS2SPositionResponse(applicant, applicantModuleName, address), source, moduleName);
    }
    delete packet;
}

const Ptr<S2SPositionResponse> GpsrStation::createS2SPositionResponse(L3Address applicant, const char* applicantModuleName, L3Address address)
{
    const auto& positionResponse = makeShared<S2SPositionResponse>();
    if (quadNodesPositionTable.hasPosition(address)){
        positionResponse->setSetted(true);
        positionResponse->setApplicant(applicant);
        positionResponse->setApplicantModuleName(applicantModuleName);
        positionResponse->setAddress(address);
        PositionTable::PositionWithTimestamp p = quadNodesPositionTable.getPositionWithTimestamp(address);
        positionResponse->setPosition(p.position);
        positionResponse->setTime(p.timestamp);
    }
    else{
        positionResponse->setSetted(false);
    }
    positionResponse->setChunkLength(B(
            1 + address.getAddressType()->getAddressByteLength()
            + applicant.getAddressType()->getAddressByteLength()
            + std::strlen(applicantModuleName)
            + positionByteLength
            + sizeof(positionResponse->getTime())));
    return positionResponse;
}


void GpsrStation::sendS2SPositionResponse(const Ptr<S2SPositionResponse>& req, L3Address dest, const char* moduleName)
{

    Packet *udpPacket = new Packet("S2SPositionResponse");
    udpPacket->insertAtBack(req);
    auto udpHeader = makeShared<UdpHeader>();
    udpHeader->setSourcePort(GPSR_UDP_PORT);
    udpHeader->setDestinationPort(GPSR_UDP_PORT);
    udpHeader->setCrcMode(CRC_DISABLED);
    udpPacket->insertAtFront(udpHeader);
    auto addresses = udpPacket->addTag<L3AddressReq>();
    addresses->setSrcAddress(getSelfAddress());
    addresses->setDestAddress(dest);
    udpPacket->addTag<HopLimitReq>()->setHopLimit(255);
    udpPacket->addTag<PacketProtocolTag>()->setProtocol(&Protocol::udp);
    udpPacket->addTag<DispatchProtocolReq>()->setProtocol(&Protocol::ipv4);
    udpPacket->addTag<InterfaceReq>()->setInterfaceId(interfaceTable->findInterfaceByName("wlan0")->getInterfaceId());
    cModule *destination = getParentModule()->getParentModule()->getSubmodule(moduleName);
    cModule *stationGpsr = destination->getSubmodule("gpsr");
    sendDirect(udpPacket, stationGpsr, "directIn");
}

void GpsrStation::processS2SPositionResponse(Packet *packet){
    packet->popAtFront<UdpHeader>();
    const auto& req = packet->peekAtFront<S2SPositionResponse>();
    bool setted = req->getSetted();
    if(setted){
        L3Address applicant = req->getApplicant();
        const char* applicantModuleName = req->getApplicantModuleName();
        L3Address address = req->getAddress();
        Coord position = req->getPosition();
        simtime_t time = req->getTime();
        sendPositionResponse(createPositionResponseFromAnotherStation(address, position, time), applicant, applicantModuleName);
    }
    delete packet;
}

const Ptr<PositionResponse> GpsrStation::createPositionResponseFromAnotherStation(L3Address address, Coord position, simtime_t timestamp)
{
    const auto& positionResponse = makeShared<PositionResponse>();
    positionResponse->setSetted(true);
    positionResponse->setAddress(address);
    positionResponse->setPosition(position);
    positionResponse->setTime(timestamp);
    positionResponse->setChunkLength(B(1 + address.getAddressType()->getAddressByteLength() + positionByteLength + sizeof(positionResponse->getTime())));
    return positionResponse;
}

// FINE MESSAGGI STAZIONI


void GpsrStation::schedulePurgeStationTimer()
{
    simtime_t nextExpiration;
    simtime_t oldestPosition = quadNodesPositionTable.getOldestPosition();
    if(quadNodesPositionTable.getAddresses().size() == 0){
        nextExpiration = simTime() + positionValidityInterval;
    }
    else if (oldestPosition == SimTime::getMaxTime()){
        nextExpiration = oldestPosition;
    }
    else{
        nextExpiration = oldestPosition + positionValidityInterval;
    }
    if (nextExpiration == SimTime::getMaxTime()) {
        if (purgeStationTimer->isScheduled())
            cancelEvent(purgeStationTimer);
    }
    else {
        if (!purgeStationTimer->isScheduled())
            scheduleAt(nextExpiration, purgeStationTimer);
        else {
            if (purgeStationTimer->getArrivalTime() != nextExpiration) {
                rescheduleAt(nextExpiration, purgeStationTimer);
            }
        }
    }
}

void GpsrStation::processPurgeStationTimer()
{
    quadNodesPositionTable.removeOldPositions(simTime() - positionValidityInterval);
    schedulePurgeStationTimer();
}

void GpsrStation::schedulePurgeDestinationTimer()
{
    simtime_t nextExpiration;
    simtime_t oldestPosition = destinationsPositionTable.getOldestPosition();
    if(destinationsPositionTable.getAddresses().size() == 0){
        nextExpiration = simTime() + destinationValidityInterval;
    }
    else if (oldestPosition == SimTime::getMaxTime()){
        nextExpiration = oldestPosition;
    }
    else{
        nextExpiration = oldestPosition + destinationValidityInterval;
    }
    if (nextExpiration == SimTime::getMaxTime()) {
        if (purgeDestinationTimer->isScheduled())
            cancelEvent(purgeDestinationTimer);
    }
    else {
        if (!purgeDestinationTimer->isScheduled())
            scheduleAt(nextExpiration, purgeDestinationTimer);
        else {
            if (purgeDestinationTimer->getArrivalTime() != nextExpiration) {
                rescheduleAt(nextExpiration, purgeDestinationTimer);
            }
        }
    }
}

void GpsrStation::processPurgeDestinationTimer()
{
    destinationsPositionTable.removeOldPositions(simTime() - destinationValidityInterval);
    schedulePurgeDestinationTimer();
}

void GpsrStation::processSendTimer()
{
    std::string hostname = receivers.front();
            L3Address destination = hostnameMap[hostname];
            receivers.pop();

            const auto& simpleMessage = makeShared<SimpleMessage>();
            const char* payload = "SimplePayload";
            simpleMessage->setPayload(payload);

            simpleMessage->setChunkLength(B(std::strlen(payload)));

            std::string msgName = "SimpleMessage" + std::to_string(msgNumber);
            msgNumber++;

            // test
            totalMsg++;

            Packet *udpPacket = new Packet(msgName.c_str());
            udpPacket->insertAtBack(simpleMessage);
            auto udpHeader = makeShared<UdpHeader>();
            udpHeader->setSourcePort(GPSR_UDP_PORT);
            udpHeader->setDestinationPort(GPSR_UDP_PORT);
            udpHeader->setCrcMode(CRC_DISABLED);
            udpPacket->insertAtFront(udpHeader);
            auto addresses = udpPacket->addTag<L3AddressReq>();
            addresses->setSrcAddress(getSelfAddress());
            addresses->setDestAddress(destination);
            udpPacket->addTag<HopLimitReq>()->setHopLimit(255);
            udpPacket->addTag<PacketProtocolTag>()->setProtocol(&Protocol::manet);
            udpPacket->addTag<DispatchProtocolReq>()->setProtocol(addressType->getNetworkProtocol());
            sendUdpPacket(udpPacket);
}
void GpsrStation::processSendTimerInterval()
{
    if (possibleReceiversMap.size() > 0){
        std::vector<L3Address> randomAddresses;
        for (const auto& pair : possibleReceiversMap) {
            randomAddresses.push_back(pair.second);
        }

        std::random_device rd;
        std::mt19937 gen(rd());
        std::shuffle(randomAddresses.begin(), randomAddresses.end(), gen);

        L3Address destination = randomAddresses.front();

        const auto& simpleMessage = makeShared<SimpleMessage>();
        const char* payload = "SimplePayload";
        simpleMessage->setPayload(payload);
        simpleMessage->setChunkLength(B(std::strlen(payload)));

        std::string msgName = "SimpleMessage" + std::to_string(msgNumber);
        msgNumber++;

        // test
        totalMsg++;

        Packet *udpPacket = new Packet(msgName.c_str());
        udpPacket->insertAtBack(simpleMessage);
        auto udpHeader = makeShared<UdpHeader>();
        udpHeader->setSourcePort(GPSR_UDP_PORT);
        udpHeader->setDestinationPort(GPSR_UDP_PORT);
        udpHeader->setCrcMode(CRC_DISABLED);
        udpPacket->insertAtFront(udpHeader);
        auto addresses = udpPacket->addTag<L3AddressReq>();
        addresses->setSrcAddress(getSelfAddress());
        addresses->setDestAddress(destination);
        udpPacket->addTag<HopLimitReq>()->setHopLimit(255);
        udpPacket->addTag<PacketProtocolTag>()->setProtocol(&Protocol::manet);
        udpPacket->addTag<DispatchProtocolReq>()->setProtocol(addressType->getNetworkProtocol());
        sendUdpPacket(udpPacket);
    }
    scheduleAfter(uniform(minInterval, maxInterval), sendTimerInterval);
}


//
// handling messages
//

void GpsrStation::processSelfMessage(cMessage *message)
{
    const char* sendTimer = "SendTimer";
    if (message == beaconTimer)
        processBeaconTimer();
    else if (message == purgeNeighborsTimer)
        processPurgeNeighborsTimer();
    else if (message == purgeStationTimer)
        processPurgeStationTimer();
    else if (message == purgeDestinationTimer)
        processPurgeDestinationTimer();
    else if (strcmp(message->getName(), sendTimer) == 0){
        processSendTimer();
        delete message;
    }
    else if (message == sendTimerInterval){
               processSendTimerInterval();
           }

    else
        throw cRuntimeError("Unknown self message");
}

void GpsrStation::processMessage(cMessage *message)
{
    if (auto pk = dynamic_cast<Packet *>(message)){
        std::string notice = "StationNotice";
        std::string noticeResponse = "StationNoticeResponse";
        std::string request = "PositionRequest";
        std::string response = "PositionResponse";
        std::string s2sRequest = "S2SPositionRequest";
        std::string s2sResponse = "S2SPositionResponse";


        const char* simpleMessage = "SimpleMessage";
        if(pk->getFullName() == notice){
            if(isStation)
                processStationNotice(pk);
            else
                delete pk;
        }
        else if(pk->getFullName()== noticeResponse){
            if(!isStation)
                processStationNoticeResponse(pk);
            else
                delete pk;
        }

        else if(pk->getFullName() == request){
            if(isStation)
                processPositionRequest(pk);
            else
                delete pk;
        }
        else if(pk->getFullName() == response){
            processPositionResponse(pk);
        }
        else if(pk->getFullName() == s2sRequest){
            if(isStation)
                processS2SPositionRequest(pk);
            else
                delete pk;
        }
        else if(pk->getFullName() == s2sResponse){
            if(isStation)
                processS2SPositionResponse(pk);
            else
                delete pk;
        }
        else if(strstr(message->getName(), simpleMessage) != nullptr){
            if(!isStation){
                            simtime_t receptionTime=simTime();
                            simtime_t sendTime=message->getCreationTime();
                            std::string name=message->getName();
                            simtime_t difference = (receptionTime - sendTime);

                            timeMessages[name]=difference.dbl();
                            packetInfo[name] = {receptionTime, mapPacketHops[name]};
                            totalArrivedMsg++;
                            delete pk;
            }
            else
                delete pk;
        }
        else{
            processUdpPacket(pk);
        }
    }
    else
        throw cRuntimeError("Unknown message");
}

//
// beacon timers
//
void GpsrStation::scheduleBeaconTimer()
{
    scheduleAfter(beaconInterval + uniform(-1, 1) * maxJitter, beaconTimer);
}

void GpsrStation::processBeaconTimer()
{    int stationIndex;
    const L3Address selfAddress = getSelfAddress();
    if (!selfAddress.isUnspecified()) {
        // Invia alla stazione la posizione
        stationIndex = getStationIndex(mobility->getCurrentPosition());
        L3Address stationAddress = stationMap[stationIndex];
        std::string stationModule = "station" + std::to_string(stationIndex);

        if (!stationAddress.isUnspecified()){
            sendStationNotice(createStationNotice(false), stationAddress, stationModule.c_str());
        }
        // Deregistrati dalla precedente (eventualmente)
        if (currentStation != stationIndex){
            L3Address currentStationAddress = stationMap[currentStation];
            std::string currentStationModule = "station" + std::to_string(currentStation);
            if (!currentStationAddress.isUnspecified()){
                sendStationNotice(createStationNotice(true), currentStationAddress, currentStationModule.c_str());
            }
            currentStation = stationIndex;
        }

    }
    scheduleBeaconTimer();
    schedulePurgeNeighborsTimer();
}

// handling purge neighbors timers
void GpsrStation::schedulePurgeNeighborsTimer()
{
    simtime_t nextExpiration = getNextNeighborExpiration();
    if (nextExpiration == SimTime::getMaxTime()) {
        if (purgeNeighborsTimer->isScheduled())
            cancelEvent(purgeNeighborsTimer);
    }
    else {
        if (!purgeNeighborsTimer->isScheduled())
            scheduleAt(nextExpiration, purgeNeighborsTimer);
        else {
            if (purgeNeighborsTimer->getArrivalTime() != nextExpiration) {
                rescheduleAt(nextExpiration, purgeNeighborsTimer);
            }
        }
    }
}

void GpsrStation::processPurgeNeighborsTimer()
{
    purgeNeighbors();
    schedulePurgeNeighborsTimer();
}

//
// handling UDP packets
//

void GpsrStation::sendUdpPacket(Packet *packet)
{
    send(packet, "ipOut");
}

void GpsrStation::processUdpPacket(Packet *packet)
{
    packet->popAtFront<UdpHeader>();
    processBeacon(packet);
    schedulePurgeNeighborsTimer();
}

//
// handling beacons
//

const Ptr<GpsrBeacon> GpsrStation::createBeacon(Coord position, const char* signature, uint64_t nonce)
{
    const auto& beacon = makeShared<GpsrBeacon>();
        // MODIFICA
        beacon->setAddress(getSelfAddress());
        beacon->setPosition(position);
        beacon->setSignature(signature);
        beacon->setNonce(nonce);
        beacon->setChunkLength(B(getSelfAddress().getAddressType()->getAddressByteLength() + positionByteLength));
    return beacon;
}

void GpsrStation::sendBeacon(const Ptr<GpsrBeacon>& beacon)
{
    Packet *udpPacket = new Packet("GPSRBeacon");
    udpPacket->insertAtBack(beacon);
    auto udpHeader = makeShared<UdpHeader>();
    udpHeader->setSourcePort(GPSR_UDP_PORT);
    udpHeader->setDestinationPort(GPSR_UDP_PORT);
    udpHeader->setCrcMode(CRC_DISABLED);
    udpPacket->insertAtFront(udpHeader);
    auto addresses = udpPacket->addTag<L3AddressReq>();
    addresses->setSrcAddress(getSelfAddress());
    addresses->setDestAddress(addressType->getLinkLocalManetRoutersMulticastAddress());
    udpPacket->addTag<HopLimitReq>()->setHopLimit(255);
    udpPacket->addTag<PacketProtocolTag>()->setProtocol(&Protocol::manet);
    udpPacket->addTag<DispatchProtocolReq>()->setProtocol(addressType->getNetworkProtocol());

    udpPacket->addTag<CreationTimeTag>()->setCreationTime(simTime());

    sendUdpPacket(udpPacket);
}


/*
 * Funzione per stimare la distanza
 * a partire dalla potenza ricevuta
 *
 * */


double GpsrStation::fromRSStoDistance(double Pr){
    double Pt= 3.01;
    double frequency = 2.4;

    double c = 3e8;
    double C = 20*log10(frequency*1e9)+20*log10(4*M_PI/c);

    double noiseFactor = 1.0 + (static_cast<double>(rand()) / RAND_MAX) * 0.1; /*dallo 0% al 10% in modo randomico*/

    return pow(10,(Pt-Pr-C)/20.0)* noiseFactor;

}



void GpsrStation::processBeacon(Packet *packet)
{
    const auto& beacon = packet->peekAtFront<GpsrBeacon>();
    EV_INFO << "Processing beacon: address = " << beacon->getAddress() << ", position = " << beacon->getPosition() << endl;
    L3Address beaconAddress = beacon->getAddress();
    Coord position= beacon->getPosition();
    uint64_t nonce = beacon->getNonce();
    std::string signature = beacon->getSignature();
    double stima=0;
    double rss;

    std::string clear = std::to_string(position.getX()) + "-" + std::to_string(position.getY()) + "-" + std::to_string(position.getZ());
    int stationIndex = getStationIndex(position);
    bool verified=false;
                if (signAlgorithm == ECDSA){

                            verified = verifySignatureECDSA(keyPairStationMapECDSA[stationIndex],clear,nonce,signature);
                        }
                        else{
                            verified = verifySignatureEDDSA(clear, signature, publicKeyMapEDDSA[stationIndex], nonce);
                        }


    if (verified){
            neighborPositionTable.setPosition(beaconAddress, position);
     }

    auto signal = packet -> getTag<inet::SignalPowerInd>();
    if(signal)
        rss = math::mW2dBmW(mW(signal->getPower()).get());
    double distance = fromRSStoDistance(rss);

    double realDistance = distanza(mobility->getCurrentPosition(),beacon->getPosition());
    double error= std::abs(realDistance - distance);
    /*
     * se l'errore tra la distanza stimata e la distanza
     * euclidea è maggiore di una soglia allora il nodo
     * sta falsando la propria posizione, di conseguenza inseriamo un valore fittizio
     * come flag in modo tale da passarlo alla stazione che andrà a modificare
     * i livelli di trustness
     *
     * */
    if(isMalicious){
        mappaDistanzaIndirizzo[beaconAddress]=-1;
    }
    else if(error > 14){
        mappaDistanzaIndirizzo[beaconAddress]=-1; //valore fittio della distanza

    }else{
    mappaDistanzaIndirizzo[beaconAddress]=distance;
    }

        if(trustness[beaconAddress] <= notTrusted1){
            delete packet;
        }
        else if(trustness[beaconAddress] >= trusted1){
            neighborPositionTable.setPosition(beacon->getAddress(), beacon->getPosition());
            delete packet;
        }
        else{
            delete packet;
        }
}



//
// handling packets
//

// TODO: handling answer position and receive response

GpsrOption *GpsrStation::createGpsrOption(L3Address destination)
{

    GpsrOption *gpsrOption = new GpsrOption();
    gpsrOption->setRoutingMode(GPSR_GREEDY_ROUTING);
    gpsrOption->setDestinationPosition(destinationsPositionTable.getPosition(destination));
    gpsrOption->setLength(computeOptionLength(gpsrOption));
    return gpsrOption;
}

int GpsrStation::computeOptionLength(GpsrOption *option)
{
    // routingMode
    int routingModeBytes = 1;
    // destinationPosition, perimeterRoutingStartPosition, perimeterRoutingForwardPosition
    int positionsBytes = 3 * positionByteLength;
    // currentFaceFirstSenderAddress, currentFaceFirstReceiverAddress, senderAddress
    int addressesBytes = 3 * getSelfAddress().getAddressType()->getAddressByteLength();
    // type and length
    int tlBytes = 1 + 1;

    return tlBytes + routingModeBytes + positionsBytes + addressesBytes;
}

//
// configuration
//

void GpsrStation::configureInterfaces()
{
    // join multicast groups
    cPatternMatcher interfaceMatcher(interfaces, false, true, false);
    for (int i = 0; i < interfaceTable->getNumInterfaces(); i++) {
        NetworkInterface *networkInterface = interfaceTable->getInterface(i);
        if (networkInterface->isMulticast() && interfaceMatcher.matches(networkInterface->getInterfaceName()))
            networkInterface->joinMulticastGroup(addressType->getLinkLocalManetRoutersMulticastAddress());
    }
}

//
// position
//

// KLUDGE implement position registry protocol
PositionTable GpsrStation::globalPositionTable;

Coord GpsrStation::lookupPositionInGlobalRegistry(const L3Address& address) const
{
    // KLUDGE implement position registry protocol
    return globalPositionTable.getPosition(address);
}

void GpsrStation::storePositionInGlobalRegistry(const L3Address& address, const Coord& position) const
{
    // KLUDGE implement position registry protocol
    globalPositionTable.setPosition(address, position);
}

void GpsrStation::storeSelfPositionInGlobalRegistry() const
{
    auto selfAddress = getSelfAddress();
    if (!selfAddress.isUnspecified())
        storePositionInGlobalRegistry(selfAddress, mobility->getCurrentPosition());
}

Coord GpsrStation::computeIntersectionInsideLineSegments(Coord& begin1, Coord& end1, Coord& begin2, Coord& end2) const
{
    // NOTE: we must explicitly avoid computing the intersection points inside due to double instability
    if (begin1 == begin2 || begin1 == end2 || end1 == begin2 || end1 == end2)
        return Coord::NIL;
    else {
        double x1 = begin1.x;
        double y1 = begin1.y;
        double x2 = end1.x;
        double y2 = end1.y;
        double x3 = begin2.x;
        double y3 = begin2.y;
        double x4 = end2.x;
        double y4 = end2.y;
        double a = determinant(x1, y1, x2, y2);
        double b = determinant(x3, y3, x4, y4);
        double c = determinant(x1 - x2, y1 - y2, x3 - x4, y3 - y4);
        double x = determinant(a, x1 - x2, b, x3 - x4) / c;
        double y = determinant(a, y1 - y2, b, y3 - y4) / c;
        if ((x <= x1 && x <= x2) || (x >= x1 && x >= x2) || (x <= x3 && x <= x4) || (x >= x3 && x >= x4) ||
            (y <= y1 && y <= y2) || (y >= y1 && y >= y2) || (y <= y3 && y <= y4) || (y >= y3 && y >= y4))
            return Coord::NIL;
        else
            return Coord(x, y, 0);
    }
}

Coord GpsrStation::getNeighborPosition(const L3Address& address) const
{
    return neighborPositionTable.getPosition(address);
}

//
// angle
//

double GpsrStation::getVectorAngle(Coord vector) const
{
    ASSERT(vector != Coord::ZERO);
    double angle = atan2(-vector.y, vector.x);
    if (angle < 0)
        angle += 2 * M_PI;
    return angle;
}

double GpsrStation::getNeighborAngle(const L3Address& address) const
{
    return getVectorAngle(getNeighborPosition(address) - mobility->getCurrentPosition());
}

//
// address
//

std::string GpsrStation::getHostName() const
{
    return host->getFullName();
}

L3Address GpsrStation::getSelfAddress() const
{
    // TODO choose self address based on a new 'interfaces' parameter
    L3Address ret = routingTable->getRouterIdAsGeneric();
#ifdef INET_WITH_IPv6
    if (ret.getType() == L3Address::IPv6) {
        for (int i = 0; i < interfaceTable->getNumInterfaces(); i++) {
            NetworkInterface *ie = interfaceTable->getInterface(i);
            if ((!ie->isLoopback())) {
                if (auto ipv6Data = ie->findProtocolData<Ipv6InterfaceData>()) {
                    ret = ipv6Data->getPreferredAddress();
                    break;
                }
            }
        }
    }
#endif
    return ret;
}

L3Address GpsrStation::getSenderNeighborAddress(const Ptr<const NetworkHeaderBase>& networkHeader) const
{
    const GpsrOption *gpsrOption = getGpsrOptionFromNetworkDatagram(networkHeader);
    return gpsrOption->getSenderAddress();
}

//
// neighbor
//

simtime_t GpsrStation::getNextNeighborExpiration()
{
    simtime_t oldestPosition = neighborPositionTable.getOldestPosition();
    if (oldestPosition == SimTime::getMaxTime())
        return oldestPosition;
    else
        return oldestPosition + neighborValidityInterval;
}

void GpsrStation::purgeNeighbors()
{
    neighborPositionTable.removeOldPositions(simTime() - neighborValidityInterval);
}

std::vector<L3Address> GpsrStation::getPlanarNeighbors() const
{
    std::vector<L3Address> planarNeighbors;
    std::vector<L3Address> neighborAddresses = neighborPositionTable.getAddresses();
    Coord selfPosition = mobility->getCurrentPosition();
    for (auto it = neighborAddresses.begin(); it != neighborAddresses.end(); it++) {
        auto neighborAddress = *it;
        Coord neighborPosition = neighborPositionTable.getPosition(neighborAddress);
        if (planarizationMode == GPSR_NO_PLANARIZATION)
            return neighborAddresses;
        else if (planarizationMode == GPSR_RNG_PLANARIZATION) {
            double neighborDistance = (neighborPosition - selfPosition).length();
            for (auto& witnessAddress : neighborAddresses) {
                Coord witnessPosition = neighborPositionTable.getPosition(witnessAddress);
                double witnessDistance = (witnessPosition - selfPosition).length();
                double neighborWitnessDistance = (witnessPosition - neighborPosition).length();
                if (neighborAddress == witnessAddress)
                    continue;
                else if (neighborDistance > std::max(witnessDistance, neighborWitnessDistance))
                    goto eliminate;
            }
        }
        else if (planarizationMode == GPSR_GG_PLANARIZATION) {
            Coord middlePosition = (selfPosition + neighborPosition) / 2;
            double neighborDistance = (neighborPosition - middlePosition).length();
            for (auto& witnessAddress : neighborAddresses) {
                Coord witnessPosition = neighborPositionTable.getPosition(witnessAddress);
                double witnessDistance = (witnessPosition - middlePosition).length();
                if (neighborAddress == witnessAddress)
                    continue;
                else if (witnessDistance < neighborDistance)
                    goto eliminate;
            }
        }
        else
            throw cRuntimeError("Unknown planarization mode");
        planarNeighbors.push_back(*it);
      eliminate:;
    }
    return planarNeighbors;
}

std::vector<L3Address> GpsrStation::getPlanarNeighborsCounterClockwise(double startAngle) const
{
    std::vector<L3Address> neighborAddresses = getPlanarNeighbors();
    std::sort(neighborAddresses.begin(), neighborAddresses.end(), [&] (const L3Address& address1, const L3Address& address2) {
        // NOTE: make sure the neighbor at startAngle goes to the end
        auto angle1 = getNeighborAngle(address1) - startAngle;
        auto angle2 = getNeighborAngle(address2) - startAngle;
        if (angle1 <= 0)
            angle1 += 2 * M_PI;
        if (angle2 <= 0)
            angle2 += 2 * M_PI;
        return angle1 < angle2;
    });
    return neighborAddresses;
}

//
// next hop
//

L3Address GpsrStation::findNextHop(const L3Address& destination, GpsrOption *gpsrOption)
{
    switch (gpsrOption->getRoutingMode()) {
        case GPSR_GREEDY_ROUTING: return findGreedyRoutingNextHop(destination, gpsrOption);
        case GPSR_PERIMETER_ROUTING: return findPerimeterRoutingNextHop(destination, gpsrOption);
        default: throw cRuntimeError("Unknown routing mode");
    }
}

L3Address GpsrStation::findGreedyRoutingNextHop(const L3Address& destination, GpsrOption *gpsrOption)
{
    auto start = std::chrono::high_resolution_clock::now();
    L3Address selfAddress = getSelfAddress();
    Coord selfPosition = mobility->getCurrentPosition();
    Coord destinationPosition = gpsrOption->getDestinationPosition();
    double bestDistance = (destinationPosition - selfPosition).length();
    L3Address bestNeighbor;
    std::vector<L3Address> neighborAddresses = neighborPositionTable.getAddresses();
    for (auto& neighborAddress : neighborAddresses) {
        Coord neighborPosition = neighborPositionTable.getPosition(neighborAddress);
        double neighborDistance = (destinationPosition - neighborPosition).length();
        if (neighborDistance < bestDistance) {
            bestDistance = neighborDistance;
            bestNeighbor = neighborAddress;
        }
    }
    if (isMalicious || bestNeighbor.isUnspecified()) {

        gpsrOption->setRoutingMode(GPSR_PERIMETER_ROUTING);
        gpsrOption->setPerimeterRoutingStartPosition(selfPosition);
        gpsrOption->setPerimeterRoutingForwardPosition(selfPosition);
        gpsrOption->setCurrentFaceFirstSenderAddress(selfAddress);
        gpsrOption->setCurrentFaceFirstReceiverAddress(L3Address());
        auto end = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> duration = end - start;
        timeInGreedy += duration.count();
        return findPerimeterRoutingNextHop(destination, gpsrOption);
    }
    else{
        auto end = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> duration = end - start;
        timeInGreedy += duration.count();
        return bestNeighbor;
    }
}
std::string verifyResponse = "VerifyPositionResponse";
L3Address GpsrStation::findPerimeterRoutingNextHop(const L3Address& destination, GpsrOption *gpsrOption)
{
    auto start = std::chrono::high_resolution_clock::now();
    L3Address selfAddress = getSelfAddress();
    Coord selfPosition = mobility->getCurrentPosition();
    Coord perimeterRoutingStartPosition = gpsrOption->getPerimeterRoutingStartPosition();
    Coord destinationPosition = gpsrOption->getDestinationPosition();
    double selfDistance = (destinationPosition - selfPosition).length();
    double perimeterRoutingStartDistance = (destinationPosition - perimeterRoutingStartPosition).length();
    if (selfDistance < perimeterRoutingStartDistance && !isMalicious) {

        gpsrOption->setRoutingMode(GPSR_GREEDY_ROUTING);
        gpsrOption->setPerimeterRoutingStartPosition(Coord());
        gpsrOption->setPerimeterRoutingForwardPosition(Coord());
        gpsrOption->setCurrentFaceFirstSenderAddress(L3Address());
        gpsrOption->setCurrentFaceFirstReceiverAddress(L3Address());
        auto end = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> duration = end - start;
        timeInPerimeter += duration.count();
        return findGreedyRoutingNextHop(destination, gpsrOption);
    }
    else {
        const L3Address& firstSenderAddress = gpsrOption->getCurrentFaceFirstSenderAddress();
        const L3Address& firstReceiverAddress = gpsrOption->getCurrentFaceFirstReceiverAddress();
        auto senderNeighborAddress = gpsrOption->getSenderAddress();
        auto neighborAngle = senderNeighborAddress.isUnspecified() ? getVectorAngle(destinationPosition - mobility->getCurrentPosition()) : getNeighborAngle(senderNeighborAddress);
        L3Address selectedNeighborAddress;
        std::vector<L3Address> neighborAddresses = getPlanarNeighborsCounterClockwise(neighborAngle);
        for (auto& neighborAddress : neighborAddresses) {
            Coord neighborPosition = getNeighborPosition(neighborAddress);
            Coord intersection = computeIntersectionInsideLineSegments(perimeterRoutingStartPosition, destinationPosition, selfPosition, neighborPosition);
            if (std::isnan(intersection.x)) {
                selectedNeighborAddress = neighborAddress;
                break;
            }
            else {
                gpsrOption->setCurrentFaceFirstSenderAddress(selfAddress);
                gpsrOption->setCurrentFaceFirstReceiverAddress(L3Address());
                gpsrOption->setPerimeterRoutingForwardPosition(intersection);
            }
        }
        if (selectedNeighborAddress.isUnspecified()) {
            auto end = std::chrono::high_resolution_clock::now();
            std::chrono::duration<double> duration = end - start;
            timeInPerimeter += duration.count();
            return L3Address();
        }
        else if (firstSenderAddress == selfAddress && firstReceiverAddress == selectedNeighborAddress) {
            auto end = std::chrono::high_resolution_clock::now();
            std::chrono::duration<double> duration = end - start;
            timeInPerimeter += duration.count();
            return L3Address();
        }
        else {
            if (gpsrOption->getCurrentFaceFirstReceiverAddress().isUnspecified())
                gpsrOption->setCurrentFaceFirstReceiverAddress(selectedNeighborAddress);
            auto end = std::chrono::high_resolution_clock::now();
            std::chrono::duration<double> duration = end - start;
            timeInPerimeter += duration.count();
            return selectedNeighborAddress;
        }
    }
}

//
// routing
//

INetfilter::IHook::Result GpsrStation::routeDatagram(Packet *datagram, GpsrOption *gpsrOption)
{
    const auto& networkHeader = getNetworkProtocolHeader(datagram);
    const L3Address& source = networkHeader->getSourceAddress();
    const L3Address& destination = networkHeader->getDestinationAddress();

    if (isMalicious){
        if (!destinationSetted){
            destinationAddress = destination;
            destinationSetted = true;
        }
        else if (destinationSetted && destinationAddress != destination){
            destinationAddress = destination;
        }
    }
    auto nextHop = findNextHop(destination, gpsrOption);
    datagram->addTagIfAbsent<NextHopAddressReq>()->setNextHopAddress(nextHop);
    if (nextHop.isUnspecified())
        return DROP;
    else {
        gpsrOption->setSenderAddress(getSelfAddress());
        auto networkInterface = CHK(interfaceTable->findInterfaceByName(outputInterface));
        datagram->addTagIfAbsent<InterfaceReq>()->setInterfaceId(networkInterface->getInterfaceId());

        mapPacketHops[datagram->getFullName()]++;

        return ACCEPT;
    }
}

void GpsrStation::setGpsrOptionOnNetworkDatagram(Packet *packet, const Ptr<const NetworkHeaderBase>& networkHeader, GpsrOption *gpsrOption)
{
    packet->trimFront();
#ifdef INET_WITH_IPv4
    if (dynamicPtrCast<const Ipv4Header>(networkHeader)) {
        auto ipv4Header = removeNetworkProtocolHeader<Ipv4Header>(packet);
        gpsrOption->setType(IPOPTION_TLV_GPSR);
        B oldHlen = ipv4Header->calculateHeaderByteLength();
        ASSERT(ipv4Header->getHeaderLength() == oldHlen);
        ipv4Header->addOption(gpsrOption);
        B newHlen = ipv4Header->calculateHeaderByteLength();
        ipv4Header->setHeaderLength(newHlen);
        ipv4Header->addChunkLength(newHlen - oldHlen);
        ipv4Header->setTotalLengthField(ipv4Header->getTotalLengthField() + newHlen - oldHlen);
        insertNetworkProtocolHeader(packet, Protocol::ipv4, ipv4Header);
    }
    else
#endif
#ifdef INET_WITH_IPv6
    if (dynamicPtrCast<const Ipv6Header>(networkHeader)) {
        auto ipv6Header = removeNetworkProtocolHeader<Ipv6Header>(packet);
        gpsrOption->setType(IPv6TLVOPTION_TLV_GPSR);
        B oldHlen = ipv6Header->calculateHeaderByteLength();
        Ipv6HopByHopOptionsHeader *hdr = check_and_cast_nullable<Ipv6HopByHopOptionsHeader *>(ipv6Header->findExtensionHeaderByTypeForUpdate(IP_PROT_IPv6EXT_HOP));
        if (hdr == nullptr) {
            hdr = new Ipv6HopByHopOptionsHeader();
            hdr->setByteLength(B(8));
            ipv6Header->addExtensionHeader(hdr);
        }
        hdr->getTlvOptionsForUpdate().appendTlvOption(gpsrOption);
        hdr->setByteLength(B(utils::roundUp(2 + B(hdr->getTlvOptions().getLength()).get(), 8)));
        B newHlen = ipv6Header->calculateHeaderByteLength();
        ipv6Header->addChunkLength(newHlen - oldHlen);
        insertNetworkProtocolHeader(packet, Protocol::ipv6, ipv6Header);
    }
    else
#endif
#ifdef INET_WITH_NEXTHOP
    if (dynamicPtrCast<const NextHopForwardingHeader>(networkHeader)) {
        auto nextHopHeader = removeNetworkProtocolHeader<NextHopForwardingHeader>(packet);
        gpsrOption->setType(NEXTHOP_TLVOPTION_TLV_GPSR);
        int oldHlen = nextHopHeader->getTlvOptions().getLength();
        nextHopHeader->getTlvOptionsForUpdate().appendTlvOption(gpsrOption);
        int newHlen = nextHopHeader->getTlvOptions().getLength();
        nextHopHeader->addChunkLength(B(newHlen - oldHlen));
        insertNetworkProtocolHeader(packet, Protocol::nextHopForwarding, nextHopHeader);
    }
    else
#endif
    {
    }
}

const GpsrOption *GpsrStation::findGpsrOptionInNetworkDatagram(const Ptr<const NetworkHeaderBase>& networkHeader) const
{
    const GpsrOption *gpsrOption = nullptr;

#ifdef INET_WITH_IPv4
    if (auto ipv4Header = dynamicPtrCast<const Ipv4Header>(networkHeader)) {
        gpsrOption = check_and_cast_nullable<const GpsrOption *>(ipv4Header->findOptionByType(IPOPTION_TLV_GPSR));
    }
    else
#endif
#ifdef INET_WITH_IPv6
    if (auto ipv6Header = dynamicPtrCast<const Ipv6Header>(networkHeader)) {
        const Ipv6HopByHopOptionsHeader *hdr = check_and_cast_nullable<const Ipv6HopByHopOptionsHeader *>(ipv6Header->findExtensionHeaderByType(IP_PROT_IPv6EXT_HOP));
        if (hdr != nullptr) {
            int i = (hdr->getTlvOptions().findByType(IPv6TLVOPTION_TLV_GPSR));
            if (i >= 0)
                gpsrOption = check_and_cast<const GpsrOption *>(hdr->getTlvOptions().getTlvOption(i));
        }
    }
    else
#endif
#ifdef INET_WITH_NEXTHOP
    if (auto nextHopHeader = dynamicPtrCast<const NextHopForwardingHeader>(networkHeader)) {
        int i = (nextHopHeader->getTlvOptions().findByType(NEXTHOP_TLVOPTION_TLV_GPSR));
        if (i >= 0)
            gpsrOption = check_and_cast<const GpsrOption *>(nextHopHeader->getTlvOptions().getTlvOption(i));
    }
    else
#endif
    {
    }
    return gpsrOption;
}

GpsrOption *GpsrStation::findGpsrOptionInNetworkDatagramForUpdate(const Ptr<NetworkHeaderBase>& networkHeader)
{
    GpsrOption *gpsrOption = nullptr;

#ifdef INET_WITH_IPv4
    if (auto ipv4Header = dynamicPtrCast<Ipv4Header>(networkHeader)) {
        gpsrOption = check_and_cast_nullable<GpsrOption *>(ipv4Header->findMutableOptionByType(IPOPTION_TLV_GPSR));
    }
    else
#endif
#ifdef INET_WITH_IPv6
    if (auto ipv6Header = dynamicPtrCast<Ipv6Header>(networkHeader)) {
        Ipv6HopByHopOptionsHeader *hdr = check_and_cast_nullable<Ipv6HopByHopOptionsHeader *>(ipv6Header->findExtensionHeaderByTypeForUpdate(IP_PROT_IPv6EXT_HOP));
        if (hdr != nullptr) {
            int i = (hdr->getTlvOptions().findByType(IPv6TLVOPTION_TLV_GPSR));
            if (i >= 0)
                gpsrOption = check_and_cast<GpsrOption *>(hdr->getTlvOptionsForUpdate().getTlvOptionForUpdate(i));
        }
    }
    else
#endif
#ifdef INET_WITH_NEXTHOP
    if (auto nextHopHeader = dynamicPtrCast<NextHopForwardingHeader>(networkHeader)) {
        int i = (nextHopHeader->getTlvOptions().findByType(NEXTHOP_TLVOPTION_TLV_GPSR));
        if (i >= 0)
            gpsrOption = check_and_cast<GpsrOption *>(nextHopHeader->getTlvOptionsForUpdate().getTlvOptionForUpdate(i));
    }
    else
#endif
    {
    }
    return gpsrOption;
}

const GpsrOption *GpsrStation::getGpsrOptionFromNetworkDatagram(const Ptr<const NetworkHeaderBase>& networkHeader) const
{
    const GpsrOption *gpsrOption = findGpsrOptionInNetworkDatagram(networkHeader);
    if (gpsrOption == nullptr)
        throw cRuntimeError("Gpsr option not found in datagram!");
    return gpsrOption;
}

GpsrOption *GpsrStation::getGpsrOptionFromNetworkDatagramForUpdate(const Ptr<NetworkHeaderBase>& networkHeader)
{
    GpsrOption *gpsrOption = findGpsrOptionInNetworkDatagramForUpdate(networkHeader);
    if (gpsrOption == nullptr)
        throw cRuntimeError("Gpsr option not found in datagram!");
    return gpsrOption;
}

//
// netfilter
//

INetfilter::IHook::Result GpsrStation::datagramPreRoutingHook(Packet *datagram)
{
    Enter_Method("datagramPreRoutingHook");
    const auto& networkHeader = getNetworkProtocolHeader(datagram);
    const L3Address& destination = networkHeader->getDestinationAddress();
    if (destination.isMulticast() || destination.isBroadcast() || routingTable->isLocalAddress(destination))
        return ACCEPT;
    else {
        const auto& ipv4Header = datagram->peekAtFront<Ipv4Header>();
        if(ipv4Header->getTimeToLive() <= 1)
            return DROP;
        auto gpsrOption = const_cast<GpsrOption *>(getGpsrOptionFromNetworkDatagram(networkHeader));
        return routeDatagram(datagram, gpsrOption);
    }
}

INetfilter::IHook::Result GpsrStation::datagramLocalOutHook(Packet *packet)
{
    Enter_Method("datagramLocalOutHook");
    const auto& networkHeader = getNetworkProtocolHeader(packet);
    const L3Address& destination = networkHeader->getDestinationAddress();
    const char* simpleMessage = "SimpleMessage";
         if (destination.isMulticast() || destination.isBroadcast() || routingTable->isLocalAddress(destination) ||
                     strstr(packet->getName(), simpleMessage) == nullptr)
                 return ACCEPT;
    else {
        if (destinationsPositionTable.hasPosition(destination)){
            GpsrOption *gpsrOption = createGpsrOption(networkHeader->getDestinationAddress());
            setGpsrOptionOnNetworkDatagram(packet, networkHeader, gpsrOption);
            return routeDatagram(packet, gpsrOption);
        }
        else{
            delayDatagram(packet);

            const L3Address selfAddress = getSelfAddress();
            if (!selfAddress.isUnspecified()) {

                int stationIndex = getStationIndex(mobility->getCurrentPosition());
                L3Address stationAddress = stationMap[stationIndex];
                std::string stationModule = "station" + std::to_string(stationIndex);

                if (!stationAddress.isUnspecified()){
                    sendPositionRequest(createPositionRequest(destination), stationAddress, stationModule.c_str());
                }
            }

            return QUEUE;
        }
    }
}

void GpsrStation::delayDatagram(Packet *datagram)
{
    const auto& networkHeader = getNetworkProtocolHeader(datagram);
    const L3Address& target = networkHeader->getDestinationAddress();
    targetAddressToDelayedPackets.insert(std::pair<L3Address, Packet *>(target, datagram));
}

void GpsrStation::sendDelayedDatagram(Packet *datagram)
{

    const auto& networkHeader = getNetworkProtocolHeader(datagram);
    const L3Address& destination = networkHeader->getDestinationAddress();
    const char* msgName = datagram->getFullName();
    datagram->popAtFront<Ipv4Header>();
    datagram->popAtFront<UdpHeader>();
    const auto& simpleMessage = datagram->peekAtFront<SimpleMessage>();
    Packet *udpPacket = new Packet(msgName);
    udpPacket->insertAtBack(simpleMessage);
    auto udpHeader = makeShared<UdpHeader>();
    udpHeader->setSourcePort(GPSR_UDP_PORT);
    udpHeader->setDestinationPort(GPSR_UDP_PORT);
    udpHeader->setCrcMode(CRC_DISABLED);
    udpPacket->insertAtFront(udpHeader);
    auto addresses = udpPacket->addTag<L3AddressReq>();
    addresses->setSrcAddress(getSelfAddress());
    addresses->setDestAddress(destination);
    udpPacket->addTag<HopLimitReq>()->setHopLimit(255);
    udpPacket->addTag<PacketProtocolTag>()->setProtocol(&Protocol::manet);
    udpPacket->addTag<DispatchProtocolReq>()->setProtocol(addressType->getNetworkProtocol());
    sendUdpPacket(udpPacket);
}

bool GpsrStation::hasDelayedDatagrams(const L3Address& target)
{
    return containsKey(targetAddressToDelayedPackets, target);
}

void GpsrStation::eraseDelayedDatagrams(const L3Address& target)
{
    auto lt = targetAddressToDelayedPackets.lower_bound(target);
    auto ut = targetAddressToDelayedPackets.upper_bound(target);
    targetAddressToDelayedPackets.erase(lt, ut);
}

//
// lifecycle
//

void GpsrStation::handleStartOperation(LifecycleOperation *operation)
{
    configureInterfaces();
    hostnameMap[host->getFullName()] = getSelfAddress();
    finished[host->getFullName()] = false;
    storeSelfPositionInGlobalRegistry();
    if(isStation){
        int myIndex = getStationIndex(mobility->getCurrentPosition());
               stationMap[myIndex] = getSelfAddress();
               if (signAlgorithm == ECDSA){
                           generateKeysECDSA();
                           keyPairStationMapECDSA[myIndex] = keyPairECDSA;
                       }
                       else{
                           generateKeysEDDSA();
                           publicKeyMapEDDSA[myIndex] = publicKeyEDDSA;
                       }
               schedulePurgeStationTimer();
    }
    else {
        auto energyStorage = getModuleByPath("^.energyStorage");
                       if (energyStorage) {
                           // Converte in un puntatore alla classe corretta
                           auto battery = dynamic_cast<inet::power::SimpleEpEnergyStorage*>(energyStorage);
                           if (battery)
                               batteriesInfo[host->getFullName()].first = getBatteryLevel();
                       }
       if(isMalicious){
            maliciousNode.push_back(getSelfAddress());
       }
       possibleReceiversMap[host->getFullName()] = getSelfAddress();
               currentStation = getStationIndex(mobility->getCurrentPosition());
               if (sendMode == INTERVAL){
               scheduleAfter(uniform(minInterval, maxInterval), sendTimerInterval);
               }
               scheduleBeaconTimer();
               schedulePurgeDestinationTimer();
           }
    }


void GpsrStation::handleStopOperation(LifecycleOperation *operation)
{
    // TODO send a beacon to remove ourself from peers neighbor position table
    destinationsPositionTable.clear();
    quadNodesPositionTable.clear();
    neighborPositionTable.clear();
    cancelEvent(beaconTimer);
    cancelEvent(purgeNeighborsTimer);
    cancelEvent(purgeStationTimer);
    cancelEvent(purgeDestinationTimer);
    cancelEvent(sendTimerInterval);
}

void GpsrStation::handleCrashOperation(LifecycleOperation *operation)
{
    destinationsPositionTable.clear();
    quadNodesPositionTable.clear();
    neighborPositionTable.clear();
    cancelEvent(beaconTimer);
    cancelEvent(purgeNeighborsTimer);
    cancelEvent(purgeStationTimer);
    cancelEvent(purgeDestinationTimer);
    cancelEvent(sendTimerInterval);
}

//
// notification
//

void GpsrStation::receiveSignal(cComponent *source, simsignal_t signalID, cObject *obj, cObject *details)
{
    Enter_Method("%s", cComponent::getSignalName(signalID));

    if (signalID == linkBrokenSignal) {
        // TODO remove the neighbor
    }
}

} // namespace inet






