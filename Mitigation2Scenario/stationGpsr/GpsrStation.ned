//
// Copyright (C) 2013 OpenSim Ltd.
//
// SPDX-License-Identifier: LGPL-3.0-or-later
//


package stationGpsr;

import inet.routing.contract.IManetRouting;

//
// This module provides Greedy Perimeter Stateless Routing for Wireless Networks.
// The implementation supports both GG and RNG planarization algorithms.
//
// For more information on the routing algorithm, see the GPSR paper
// http://www.eecs.harvard.edu/~~htk/publication/2000-mobi-karp-kung.pdf
//
simple GpsrStation like IManetRouting
{
    parameters:
        // properties
        @display("i=block/routing");

        // context parameters
        string interfaceTableModule;   // The path to the InterfaceTable module
        string routingTableModule = default("^.ipv4.routingTable");
        string networkProtocolModule = default("^.ipv4.ip");
        string outputInterface = default("wlan0");

        // GPSR parameters
        string planarizationMode @enum("", "GG", "RNG") = default("GG");
        string interfaces = default("*");
        double beaconInterval @unit(s) = default(10s);
        double maxJitter @unit(s) = default(0.5 * beaconInterval);
        double neighborValidityInterval @unit(s) = default(4.5 * beaconInterval);
        int positionByteLength @unit(B) = default(2 * 4B);

        // visualization parameters
        bool displayBubbles = default(false);   // display bubble messages about changes in routing state for packets
        
        bool isMalicious = default(false);
        bool isStation = default(false);
        int n = default(1000);
        int m = default(3);
        double minInterval @unit(s) = default(5s); 
        double maxInterval @unit(s) = default(10s); 
        
        string sendMode @enum("FixedMaximum", "Interval") = default("Interval");
        string signAlgorithm @enum("ECDSA", "EDDSA") = default("ECDSA");
        double positionValidityInterval @unit(s) = default(4.5 * beaconInterval);
        double destinationValidityInterval @unit(s) = default(4.5 * beaconInterval);
        
        string sendTimes = default("");
        string receivers = default("");
    gates:
        input ipIn;
        input directIn @directIn;
        output ipOut;
}

