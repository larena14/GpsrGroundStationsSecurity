//
// Copyright (C) 2013 OpenSim Ltd.
//
// SPDX-License-Identifier: LGPL-3.0-or-later
//


#ifndef __INET_CUSTOMSTATEBASEDEPENERGYCONSUMER_H
#define __INET_CUSTOMSTATEBASEDEPENERGYCONSUMER_H

#include "inet/common/ModuleRefByPar.h"
#include "inet/physicallayer/wireless/common/contract/packetlevel/IRadio.h"
#include "inet/power/contract/IEpEnergyConsumer.h"
#include "inet/power/contract/IEpEnergySource.h"
#include "inet/power/contract/IEnergyConsumer.h"



namespace inet {


/**
 * This is a simple radio power consumer model. The power consumption is
 * determined by the radio mode, the transmitter state and the receiver state
 * using constant parameters.
 *
 */
class INET_API CustomStateBasedEpEnergyConsumer : public cSimpleModule, public power::IEpEnergyConsumer, public cListener
{
  protected:
    W offPowerConsumption = W(NaN);
    W sleepPowerConsumption = W(NaN);
    W switchingPowerConsumption = W(NaN);
    W receiverIdlePowerConsumption = W(NaN);
    W receiverBusyPowerConsumption = W(NaN);
    W receiverReceivingPowerConsumption = W(NaN);
    W receiverReceivingPreamblePowerConsumption = W(NaN);
    W receiverReceivingHeaderPowerConsumption = W(NaN);
    W receiverReceivingDataPowerConsumption = W(NaN);
    W transmitterIdlePowerConsumption = W(NaN);
    W transmitterTransmittingPowerConsumption = W(NaN);
    W transmitterTransmittingPreamblePowerConsumption = W(NaN);
    W transmitterTransmittingHeaderPowerConsumption = W(NaN);
    W transmitterTransmittingDataPowerConsumption = W(NaN);
    W hoveringPowerConsumption = W(NaN);
    W maxMovingPowerConsumption = W(NaN);
    W signPowerConsumption = W(NaN);
    W verifyPowerConsumption = W(NaN);
    opp_component_ptr<physicallayer::IRadio> radio;
    ModuleRefByPar<power::IEpEnergySource> energySource;
    W powerConsumption = W(NaN);
    opp_component_ptr<IMobility> mobility;
    Coord lastPosition;

  protected:
    virtual int numInitStages() const override { return NUM_INIT_STAGES; }
    virtual void initialize(int stage) override;
    virtual W computePowerConsumption() const;
    virtual W computeMovementConsumption() const;
    virtual W computeSignConsumption(double bytes) const;
    virtual W computeVerifyConsumption(double bytes) const;

  public:
    virtual power::IEnergySource *getEnergySource() const override { return energySource; }
    virtual W getPowerConsumption() const override { return powerConsumption; }

    virtual void receiveSignal(cComponent *source, simsignal_t signal, intval_t value, cObject *details) override;
    virtual void receiveSignal(cComponent *source, simsignal_t signal, cObject *obj, cObject *details) override;
    virtual void receiveSignal(cComponent *source, simsignal_t signal, double b, cObject *details) override;


};

}

#endif

