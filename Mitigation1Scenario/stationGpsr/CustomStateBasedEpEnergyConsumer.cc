//
// Copyright (C) 2013 OpenSim Ltd.
//
// SPDX-License-Identifier: LGPL-3.0-or-later
//


#include "CustomStateBasedEpEnergyConsumer.h"

#include "inet/common/ModuleAccess.h"
#include "GpsrStation.h"

namespace inet {

using namespace inet::power;

Define_Module(CustomStateBasedEpEnergyConsumer);

void CustomStateBasedEpEnergyConsumer::initialize(int stage)
{
    cSimpleModule::initialize(stage);
    if (stage == INITSTAGE_LOCAL) {
        offPowerConsumption = W(par("offPowerConsumption"));
        sleepPowerConsumption = W(par("sleepPowerConsumption"));
        switchingPowerConsumption = W(par("switchingPowerConsumption"));
        receiverIdlePowerConsumption = W(par("receiverIdlePowerConsumption"));
        receiverBusyPowerConsumption = W(par("receiverBusyPowerConsumption"));
        receiverReceivingPowerConsumption = W(par("receiverReceivingPowerConsumption"));
        receiverReceivingPreamblePowerConsumption = W(par("receiverReceivingPreamblePowerConsumption"));
        receiverReceivingHeaderPowerConsumption = W(par("receiverReceivingHeaderPowerConsumption"));
        receiverReceivingDataPowerConsumption = W(par("receiverReceivingDataPowerConsumption"));
        transmitterIdlePowerConsumption = W(par("transmitterIdlePowerConsumption"));
        transmitterTransmittingPowerConsumption = W(par("transmitterTransmittingPowerConsumption"));
        transmitterTransmittingPreamblePowerConsumption = W(par("transmitterTransmittingPreamblePowerConsumption"));
        transmitterTransmittingHeaderPowerConsumption = W(par("transmitterTransmittingHeaderPowerConsumption"));
        transmitterTransmittingDataPowerConsumption = W(par("transmitterTransmittingDataPowerConsumption"));
        hoveringPowerConsumption = W(par("hoveringPowerConsumption"));
        maxMovingPowerConsumption = W(par("maxMovingPowerConsumption"));
        signPowerConsumption = W(par("signPowerConsumption"));
        verifyPowerConsumption = W(par("verifyPowerConsumption"));
        cModule *radioModule = getParentModule();
        radioModule->subscribe(physicallayer::IRadio::radioModeChangedSignal, this);
        radioModule->subscribe(physicallayer::IRadio::receptionStateChangedSignal, this);
        radioModule->subscribe(physicallayer::IRadio::transmissionStateChangedSignal, this);
        radioModule->subscribe(physicallayer::IRadio::receivedSignalPartChangedSignal, this);
        radioModule->subscribe(physicallayer::IRadio::transmittedSignalPartChangedSignal, this);
        radio = check_and_cast<physicallayer::IRadio *>(radioModule);

        cModule* mobilityModule = getParentModule()->getParentModule()->getParentModule()->getSubmodule("mobility");
        mobilityModule->subscribe(inet::IMobility::mobilityStateChangedSignal, this);

        mobility = check_and_cast<IMobility *>(mobilityModule);

        cModule* gpsrModule = getParentModule()->getParentModule()->getParentModule()->getSubmodule("gpsr");
        gpsrModule->subscribe(inet::GpsrStation::signSignal, this);
        gpsrModule->subscribe(inet::GpsrStation::verifySignal, this);

        lastPosition = mobility->getCurrentPosition();

        powerConsumption = W(0);
        energySource.reference(this, "energySourceModule", true);
        WATCH(powerConsumption);
    }
    else if (stage == INITSTAGE_POWER)
        energySource->addEnergyConsumer(this);
}

void CustomStateBasedEpEnergyConsumer::receiveSignal(cComponent *source, simsignal_t signal, cObject *obj, cObject *details){
    Enter_Method("%s", cComponent::getSignalName(signal));
    //metodo per la gestione del consumo riguardante la mobilità
    //std::cout << "receiveSignal cObject from: " << source->getFullName() << endl;
    if(signal == inet::IMobility::mobilityStateChangedSignal){
        powerConsumption = computeMovementConsumption();
        emit(powerConsumptionChangedSignal, powerConsumption.get());
    }
    else{
        throw cRuntimeError("Unknown signal, MOV");
    }
}

W CustomStateBasedEpEnergyConsumer::computeMovementConsumption() const
{
    double speed = mobility->getCurrentVelocity().length();
    //std::cout << mobility->getCurrentVelocity().length() << endl;
    //std::cout << "Consuming for movement: " << W(10*speed) + hoveringPowerConsumption << endl;
    return W(10*speed) + hoveringPowerConsumption;
}

void CustomStateBasedEpEnergyConsumer::receiveSignal(cComponent *source, simsignal_t signal, double b, cObject *details){
    Enter_Method("%s", cComponent::getSignalName(signal));
    //metodo per la gestione del consumo riguardante crittografia
    //std::cout << "receiveSignal bool " << source->getFullName() << endl;
    //std::cout << "receiveSignal cObject from: " << source->getFullName() << endl;
    if(signal == inet::GpsrStation::signSignal){
        powerConsumption = computeSignConsumption(b);
        emit(powerConsumptionChangedSignal, powerConsumption.get());
    }
    else if(signal == inet::GpsrStation::verifySignal){
        powerConsumption = computeVerifyConsumption(b);
        emit(powerConsumptionChangedSignal, powerConsumption.get());
    }
    else{
        throw cRuntimeError("Unknown signal, CRYPT");
    }
}

W CustomStateBasedEpEnergyConsumer::computeSignConsumption(double bytes) const
{
    //std::cout << "Consuming for signing: " << W(bytes*signPowerConsumption) << endl;
    return W(bytes*signPowerConsumption);
}

W CustomStateBasedEpEnergyConsumer::computeVerifyConsumption(double bytes) const
{
    //std::cout << "Consuming for verifying: " << W(bytes*signPowerConsumption) << endl;
    return W(bytes*verifyPowerConsumption);
}



void CustomStateBasedEpEnergyConsumer::receiveSignal(cComponent *source, simsignal_t signal, intval_t value, cObject *details)
{
    Enter_Method("%s", cComponent::getSignalName(signal));

    if (signal == physicallayer::IRadio::radioModeChangedSignal ||
        signal == physicallayer::IRadio::receptionStateChangedSignal ||
        signal == physicallayer::IRadio::transmissionStateChangedSignal ||
        signal == physicallayer::IRadio::receivedSignalPartChangedSignal ||
        signal == physicallayer::IRadio::transmittedSignalPartChangedSignal)
    {
        powerConsumption = computePowerConsumption();
        emit(powerConsumptionChangedSignal, powerConsumption.get());
    }
    else{
        throw cRuntimeError("Unknown signal, RAD");
    }
}

W CustomStateBasedEpEnergyConsumer::computePowerConsumption() const
{
    physicallayer::IRadio::RadioMode radioMode = radio->getRadioMode();
    if (radioMode == physicallayer::IRadio::RADIO_MODE_OFF)
        return offPowerConsumption;
    else if (radioMode == physicallayer::IRadio::RADIO_MODE_SLEEP)
        return sleepPowerConsumption;
    else if (radioMode == physicallayer::IRadio::RADIO_MODE_SWITCHING)
        return switchingPowerConsumption;
    W powerConsumption = W(0);
    physicallayer::IRadio::ReceptionState receptionState = radio->getReceptionState();
    physicallayer::IRadio::TransmissionState transmissionState = radio->getTransmissionState();
    if (radioMode == physicallayer::IRadio::RADIO_MODE_RECEIVER || radioMode == physicallayer::IRadio::RADIO_MODE_TRANSCEIVER) {
        switch (receptionState) {
            case physicallayer::IRadio::RECEPTION_STATE_IDLE:
                powerConsumption += receiverIdlePowerConsumption;
                break;
            case physicallayer::IRadio::RECEPTION_STATE_BUSY:
                powerConsumption += receiverBusyPowerConsumption;
                break;
            case physicallayer::IRadio::RECEPTION_STATE_RECEIVING: {
                auto part = radio->getReceivedSignalPart();
                switch (part) {
                    case physicallayer::IRadioSignal::SIGNAL_PART_NONE:
                        break;
                    case physicallayer::IRadioSignal::SIGNAL_PART_WHOLE:
                        powerConsumption += receiverReceivingPowerConsumption;
                        break;
                    case physicallayer::IRadioSignal::SIGNAL_PART_PREAMBLE:
                        powerConsumption += receiverReceivingPreamblePowerConsumption;
                        break;
                    case physicallayer::IRadioSignal::SIGNAL_PART_HEADER:
                        powerConsumption += receiverReceivingHeaderPowerConsumption;
                        break;
                    case physicallayer::IRadioSignal::SIGNAL_PART_DATA:
                        powerConsumption += receiverReceivingDataPowerConsumption;
                        break;
                    default:
                        throw cRuntimeError("Unknown received signal part");
                }
                break;
            }
            case physicallayer::IRadio::RECEPTION_STATE_UNDEFINED:
                break;
            default:
                throw cRuntimeError("Unknown radio reception state");
        }
    }
    if (radioMode == physicallayer::IRadio::RADIO_MODE_TRANSMITTER || radioMode == physicallayer::IRadio::RADIO_MODE_TRANSCEIVER) {
        switch (transmissionState) {
            case physicallayer::IRadio::TRANSMISSION_STATE_IDLE:
                powerConsumption += transmitterIdlePowerConsumption;
                break;
            case physicallayer::IRadio::TRANSMISSION_STATE_TRANSMITTING: {
                auto part = radio->getTransmittedSignalPart();
                switch (part) {
                    case physicallayer::IRadioSignal::SIGNAL_PART_NONE:
                        break;
                    case physicallayer::IRadioSignal::SIGNAL_PART_WHOLE:
                        powerConsumption += transmitterTransmittingPowerConsumption;
                        break;
                    case physicallayer::IRadioSignal::SIGNAL_PART_PREAMBLE:
                        powerConsumption += transmitterTransmittingPreamblePowerConsumption;
                        break;
                    case physicallayer::IRadioSignal::SIGNAL_PART_HEADER:
                        powerConsumption += transmitterTransmittingHeaderPowerConsumption;
                        break;
                    case physicallayer::IRadioSignal::SIGNAL_PART_DATA:
                        powerConsumption += transmitterTransmittingDataPowerConsumption;
                        break;
                    default:
                        throw cRuntimeError("Unknown transmitted signal part");
                }
                break;
            }
            case physicallayer::IRadio::TRANSMISSION_STATE_UNDEFINED:
                break;
            default:
                throw cRuntimeError("Unknown radio transmission state");
        }
    }
    return powerConsumption;
}

} // namespace inet

