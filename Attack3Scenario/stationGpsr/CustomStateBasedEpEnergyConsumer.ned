//
// Copyright (C) 2013 OpenSim Ltd.
//
// SPDX-License-Identifier: LGPL-3.0-or-later
//


package stationGpsr;

import inet.power.contract.IEpEnergyConsumer;

//
// This module provides a radio power consumer model. The power consumption is
// determined by the radio mode, the transmitter state and the receiver state
// using constant parameters.
//
simple CustomStateBasedEpEnergyConsumer like IEpEnergyConsumer
{
    parameters:
        string energySourceModule;
        double offPowerConsumption @unit(W) = default(0mW);
        double sleepPowerConsumption @unit(W) = default(1mW);
        double switchingPowerConsumption @unit(W) = default(1mW);
        double receiverIdlePowerConsumption @unit(W) = default(2mW);
        double receiverBusyPowerConsumption @unit(W) = default(5mW);
        double receiverReceivingPowerConsumption @unit(W) = default(10mW);
        double receiverReceivingPreamblePowerConsumption @unit(W) = default(receiverReceivingPowerConsumption);
        double receiverReceivingHeaderPowerConsumption @unit(W) = default(receiverReceivingPowerConsumption);
        double receiverReceivingDataPowerConsumption @unit(W) = default(receiverReceivingPowerConsumption);
        double transmitterIdlePowerConsumption @unit(W) = default(2mW);
        double transmitterTransmittingPowerConsumption @unit(W) = default(100mW);
        double transmitterTransmittingPreamblePowerConsumption @unit(W) = default(transmitterTransmittingPowerConsumption);
        double transmitterTransmittingHeaderPowerConsumption @unit(W) = default(transmitterTransmittingPowerConsumption);
        double transmitterTransmittingDataPowerConsumption @unit(W) = default(transmitterTransmittingPowerConsumption);
        double maxMovingPowerConsumption @unit(W) = default(250W);
        double hoveringPowerConsumption @unit(W) = default(150W);
        double signPowerConsumption @unit(W) = default(2mW);
        double verifyPowerConsumption @unit(W) = default(1.5mW);
        @display("i=block/cogwheel");
        @class(CustomStateBasedEpEnergyConsumer);
        @signal[powerConsumptionChanged](type=double);
        @statistic[powerConsumption](title="Power consumption"; source=powerConsumptionChanged; record=vector; interpolationmode=sample-hold; unit=W);
}

