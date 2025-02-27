package scenario1000;

import inet.common.scenario.ScenarioManager;
import inet.networklayer.configurator.ipv4.Ipv4NetworkConfigurator;
import inet.networklayer.ipv4.RoutingTableRecorder;
import GpsrRouter;
import inet.physicallayer.wireless.ieee80211.packetlevel.Ieee80211ScalarRadioMedium;
import inet.visualizer.common.IntegratedVisualizer;

network DronesNetwork1000_20
{
    parameters:
        @display("bgb=1000,1000");
    submodules:
        radioMedium: Ieee80211ScalarRadioMedium {
            parameters:
                @display("p=1200,211;is=s");
        }
        visualizer: IntegratedVisualizer {
            @display("p=1200,512");
        }
        configurator: Ipv4NetworkConfigurator {
            parameters:
                config = xml("<config><interface hosts='*' address='145.236.x.x' netmask='255.255.0.0'/></config>");
                @display("p=1200,110;is=s");
        }
        routingTableRecorder: RoutingTableRecorder {
            @display("p=1200,311");
        }
        scenarioManager: ScenarioManager {
            parameters:
                script = default(xml("<scenario/>"));
                @display("p=1200,412;is=s");
        }
        node0: GpsrRouter {
            gpsr.isMalicious = false;
            gpsr.isStation = false;
            gpsr.n = 1000;
            gpsr.m = 3;
            gpsr.sendTimes = "5.0s 5.6s 5.8s 6.4s 6.9s 7.0s";
            gpsr.receivers = "node44 node45 node47 node15 node56 node67";
            @display("p=0.0,0.0;i=misc/node_vs");
        }
        node1: GpsrRouter {
            gpsr.isMalicious = true;
            gpsr.isStation = false;
            gpsr.n = 1000;
            gpsr.m = 3;
            gpsr.sendTimes = "7.5s 7.9s";
            gpsr.receivers = "node80 node11";
            @display("p=0.0,111.11;i=misc/node_vs");
        }
        node2: GpsrRouter {
            gpsr.isMalicious = true;
            gpsr.isStation = false;
            gpsr.n = 1000;
            gpsr.m = 3;
            gpsr.sendTimes = "8.4s 9.1s";
            gpsr.receivers = "node13 node94";
            @display("p=0.0,222.22;i=misc/node_vs");
        }
        node3: GpsrRouter {
            gpsr.isMalicious = false;
            gpsr.isStation = false;
            gpsr.n = 1000;
            gpsr.m = 3;
            gpsr.sendTimes = "9.6s";
            gpsr.receivers = "node79";
            @display("p=0.0,333.33;i=misc/node_vs");
        }
        node4: GpsrRouter {
            gpsr.isMalicious = true;
            gpsr.isStation = false;
            gpsr.n = 1000;
            gpsr.m = 3;
            gpsr.sendTimes = "9.8s 10.6s 10.8s 10.9s 11.2s 11.3s";
            gpsr.receivers = "node87 node98 node27 node47 node22 node5";
            @display("p=0.0,444.44;i=misc/node_vs");
        }
        node5: GpsrRouter {
            gpsr.isMalicious = false;
            gpsr.isStation = false;
            gpsr.n = 1000;
            gpsr.m = 3;
            gpsr.sendTimes = "12.7s";
            gpsr.receivers = "node10";
            @display("p=0.0,555.56;i=misc/node_vs");
        }
        node6: GpsrRouter {
            gpsr.isMalicious = false;
            gpsr.isStation = false;
            gpsr.n = 1000;
            gpsr.m = 3;
            gpsr.sendTimes = "13.2s 13.4s";
            gpsr.receivers = "node38 node32";
            @display("p=0.0,666.67;i=misc/node_vs");
        }
        node7: GpsrRouter {
            gpsr.isMalicious = false;
            gpsr.isStation = false;
            gpsr.n = 1000;
            gpsr.m = 3;
            gpsr.sendTimes = "";
            gpsr.receivers = "";
            @display("p=0.0,777.78;i=misc/node_vs");
        }
        node8: GpsrRouter {
            gpsr.isMalicious = true;
            gpsr.isStation = false;
            gpsr.n = 1000;
            gpsr.m = 3;
            gpsr.sendTimes = "13.9s";
            gpsr.receivers = "node43";
            @display("p=0.0,888.89;i=misc/node_vs");
        }
        node9: GpsrRouter {
            gpsr.isMalicious = false;
            gpsr.isStation = false;
            gpsr.n = 1000;
            gpsr.m = 3;
            gpsr.sendTimes = "14.3s 14.5s";
            gpsr.receivers = "node90 node51";
            @display("p=0.0,1000.0;i=misc/node_vs");
        }
        node10: GpsrRouter {
            gpsr.isMalicious = false;
            gpsr.isStation = false;
            gpsr.n = 1000;
            gpsr.m = 3;
            gpsr.sendTimes = "14.9s";
            gpsr.receivers = "node36";
            @display("p=111.11,0.0;i=misc/node_vs");
        }
        node11: GpsrRouter {
            gpsr.isMalicious = true;
            gpsr.isStation = false;
            gpsr.n = 1000;
            gpsr.m = 3;
            gpsr.sendTimes = "15.3s";
            gpsr.receivers = "node93";
            @display("p=111.11,111.11;i=misc/node_vs");
        }
        node12: GpsrRouter {
            gpsr.isMalicious = false;
            gpsr.isStation = false;
            gpsr.n = 1000;
            gpsr.m = 3;
            gpsr.sendTimes = "15.4s 15.6s";
            gpsr.receivers = "node43 node44";
            @display("p=111.11,222.22;i=misc/node_vs");
        }
        node13: GpsrRouter {
            gpsr.isMalicious = false;
            gpsr.isStation = false;
            gpsr.n = 1000;
            gpsr.m = 3;
            gpsr.sendTimes = "16.1s 16.7s 17.0s";
            gpsr.receivers = "node82 node20 node47";
            @display("p=111.11,333.33;i=misc/node_vs");
        }
        node14: GpsrRouter {
            gpsr.isMalicious = false;
            gpsr.isStation = false;
            gpsr.n = 1000;
            gpsr.m = 3;
            gpsr.sendTimes = "17.2s";
            gpsr.receivers = "node81";
            @display("p=111.11,444.44;i=misc/node_vs");
        }
        node15: GpsrRouter {
            gpsr.isMalicious = false;
            gpsr.isStation = false;
            gpsr.n = 1000;
            gpsr.m = 3;
            gpsr.sendTimes = "18.0s 18.3s 18.4s";
            gpsr.receivers = "node32 node45 node19";
            @display("p=111.11,555.56;i=misc/node_vs");
        }
        node16: GpsrRouter {
            gpsr.isMalicious = false;
            gpsr.isStation = false;
            gpsr.n = 1000;
            gpsr.m = 3;
            gpsr.sendTimes = "18.7s 19.0s 19.1s 19.3s";
            gpsr.receivers = "node93 node72 node67 node48";
            @display("p=111.11,666.67;i=misc/node_vs");
        }
        node17: GpsrRouter {
            gpsr.isMalicious = false;
            gpsr.isStation = false;
            gpsr.n = 1000;
            gpsr.m = 3;
            gpsr.sendTimes = "";
            gpsr.receivers = "";
            @display("p=111.11,777.78;i=misc/node_vs");
        }
        node18: GpsrRouter {
            gpsr.isMalicious = false;
            gpsr.isStation = false;
            gpsr.n = 1000;
            gpsr.m = 3;
            gpsr.sendTimes = "19.7s";
            gpsr.receivers = "node47";
            @display("p=111.11,888.89;i=misc/node_vs");
        }
        node19: GpsrRouter {
            gpsr.isMalicious = false;
            gpsr.isStation = false;
            gpsr.n = 1000;
            gpsr.m = 3;
            gpsr.sendTimes = "22.8s 22.9s";
            gpsr.receivers = "node66 node34";
            @display("p=111.11,1000.0;i=misc/node_vs");
        }
        node20: GpsrRouter {
            gpsr.isMalicious = true;
            gpsr.isStation = false;
            gpsr.n = 1000;
            gpsr.m = 3;
            gpsr.sendTimes = "23.4s";
            gpsr.receivers = "node25";
            @display("p=222.22,0.0;i=misc/node_vs");
        }
        node21: GpsrRouter {
            gpsr.isMalicious = true;
            gpsr.isStation = false;
            gpsr.n = 1000;
            gpsr.m = 3;
            gpsr.sendTimes = "";
            gpsr.receivers = "";
            @display("p=222.22,111.11;i=misc/node_vs");
        }
        node22: GpsrRouter {
            gpsr.isMalicious = true;
            gpsr.isStation = false;
            gpsr.n = 1000;
            gpsr.m = 3;
            gpsr.sendTimes = "23.6s 23.8s 23.9s";
            gpsr.receivers = "node71 node86 node40";
            @display("p=222.22,222.22;i=misc/node_vs");
        }
        node23: GpsrRouter {
            gpsr.isMalicious = true;
            gpsr.isStation = false;
            gpsr.n = 1000;
            gpsr.m = 3;
            gpsr.sendTimes = "24.3s";
            gpsr.receivers = "node37";
            @display("p=222.22,333.33;i=misc/node_vs");
        }
        node24: GpsrRouter {
            gpsr.isMalicious = false;
            gpsr.isStation = false;
            gpsr.n = 1000;
            gpsr.m = 3;
            gpsr.sendTimes = "24.8s 24.9s 25.1s";
            gpsr.receivers = "node60 node32 node31";
            @display("p=222.22,444.44;i=misc/node_vs");
        }
        node25: GpsrRouter {
            gpsr.isMalicious = false;
            gpsr.isStation = false;
            gpsr.n = 1000;
            gpsr.m = 3;
            gpsr.sendTimes = "26.2s 27.1s";
            gpsr.receivers = "node42 node64";
            @display("p=222.22,555.56;i=misc/node_vs");
        }
        node26: GpsrRouter {
            gpsr.isMalicious = false;
            gpsr.isStation = false;
            gpsr.n = 1000;
            gpsr.m = 3;
            gpsr.sendTimes = "27.2s";
            gpsr.receivers = "node54";
            @display("p=222.22,666.67;i=misc/node_vs");
        }
        node27: GpsrRouter {
            gpsr.isMalicious = false;
            gpsr.isStation = false;
            gpsr.n = 1000;
            gpsr.m = 3;
            gpsr.sendTimes = "";
            gpsr.receivers = "";
            @display("p=222.22,777.78;i=misc/node_vs");
        }
        node28: GpsrRouter {
            gpsr.isMalicious = false;
            gpsr.isStation = false;
            gpsr.n = 1000;
            gpsr.m = 3;
            gpsr.sendTimes = "27.8s";
            gpsr.receivers = "node92";
            @display("p=222.22,888.89;i=misc/node_vs");
        }
        node29: GpsrRouter {
            gpsr.isMalicious = false;
            gpsr.isStation = false;
            gpsr.n = 1000;
            gpsr.m = 3;
            gpsr.sendTimes = "27.9s 28.3s";
            gpsr.receivers = "node43 node61";
            @display("p=222.22,1000.0;i=misc/node_vs");
        }
        node30: GpsrRouter {
            gpsr.isMalicious = false;
            gpsr.isStation = false;
            gpsr.n = 1000;
            gpsr.m = 3;
            gpsr.sendTimes = "";
            gpsr.receivers = "";
            @display("p=333.33,0.0;i=misc/node_vs");
        }
        node31: GpsrRouter {
            gpsr.isMalicious = false;
            gpsr.isStation = false;
            gpsr.n = 1000;
            gpsr.m = 3;
            gpsr.sendTimes = "28.4s";
            gpsr.receivers = "node73";
            @display("p=333.33,111.11;i=misc/node_vs");
        }
        node32: GpsrRouter {
            gpsr.isMalicious = false;
            gpsr.isStation = false;
            gpsr.n = 1000;
            gpsr.m = 3;
            gpsr.sendTimes = "29.8s";
            gpsr.receivers = "node81";
            @display("p=333.33,222.22;i=misc/node_vs");
        }
        node33: GpsrRouter {
            gpsr.isMalicious = true;
            gpsr.isStation = false;
            gpsr.n = 1000;
            gpsr.m = 3;
            gpsr.sendTimes = "30.0s 30.2s";
            gpsr.receivers = "node90 node73";
            @display("p=333.33,333.33;i=misc/node_vs");
        }
        node34: GpsrRouter {
            gpsr.isMalicious = false;
            gpsr.isStation = false;
            gpsr.n = 1000;
            gpsr.m = 3;
            gpsr.sendTimes = "30.3s";
            gpsr.receivers = "node37";
            @display("p=333.33,444.44;i=misc/node_vs");
        }
        node35: GpsrRouter {
            gpsr.isMalicious = false;
            gpsr.isStation = false;
            gpsr.n = 1000;
            gpsr.m = 3;
            gpsr.sendTimes = "30.6s 30.8s";
            gpsr.receivers = "node65 node48";
            @display("p=333.33,555.56;i=misc/node_vs");
        }
        node36: GpsrRouter {
            gpsr.isMalicious = false;
            gpsr.isStation = false;
            gpsr.n = 1000;
            gpsr.m = 3;
            gpsr.sendTimes = "";
            gpsr.receivers = "";
            @display("p=333.33,666.67;i=misc/node_vs");
        }
        node37: GpsrRouter {
            gpsr.isMalicious = false;
            gpsr.isStation = false;
            gpsr.n = 1000;
            gpsr.m = 3;
            gpsr.sendTimes = "";
            gpsr.receivers = "";
            @display("p=333.33,777.78;i=misc/node_vs");
        }
        node38: GpsrRouter {
            gpsr.isMalicious = false;
            gpsr.isStation = false;
            gpsr.n = 1000;
            gpsr.m = 3;
            gpsr.sendTimes = "";
            gpsr.receivers = "";
            @display("p=333.33,888.89;i=misc/node_vs");
        }
        node39: GpsrRouter {
            gpsr.isMalicious = false;
            gpsr.isStation = false;
            gpsr.n = 1000;
            gpsr.m = 3;
            gpsr.sendTimes = "30.9s";
            gpsr.receivers = "node81";
            @display("p=333.33,1000.0;i=misc/node_vs");
        }
        node40: GpsrRouter {
            gpsr.isMalicious = false;
            gpsr.isStation = false;
            gpsr.n = 1000;
            gpsr.m = 3;
            gpsr.sendTimes = "31.1s 31.2s 31.3s";
            gpsr.receivers = "node51 node72 node87";
            @display("p=444.44,0.0;i=misc/node_vs");
        }
        node41: GpsrRouter {
            gpsr.isMalicious = false;
            gpsr.isStation = false;
            gpsr.n = 1000;
            gpsr.m = 3;
            gpsr.sendTimes = "31.7s";
            gpsr.receivers = "node74";
            @display("p=444.44,111.11;i=misc/node_vs");
        }
        node42: GpsrRouter {
            gpsr.isMalicious = false;
            gpsr.isStation = false;
            gpsr.n = 1000;
            gpsr.m = 3;
            gpsr.sendTimes = "32.3s 32.5s";
            gpsr.receivers = "node70 node69";
            @display("p=444.44,222.22;i=misc/node_vs");
        }
        node43: GpsrRouter {
            gpsr.isMalicious = false;
            gpsr.isStation = false;
            gpsr.n = 1000;
            gpsr.m = 3;
            gpsr.sendTimes = "";
            gpsr.receivers = "";
            @display("p=444.44,333.33;i=misc/node_vs");
        }
        node44: GpsrRouter {
            gpsr.isMalicious = false;
            gpsr.isStation = false;
            gpsr.n = 1000;
            gpsr.m = 3;
            gpsr.sendTimes = "";
            gpsr.receivers = "";
            @display("p=444.44,444.44;i=misc/node_vs");
        }
        node45: GpsrRouter {
            gpsr.isMalicious = false;
            gpsr.isStation = false;
            gpsr.n = 1000;
            gpsr.m = 3;
            gpsr.sendTimes = "";
            gpsr.receivers = "";
            @display("p=444.44,555.56;i=misc/node_vs");
        }
        node46: GpsrRouter {
            gpsr.isMalicious = false;
            gpsr.isStation = false;
            gpsr.n = 1000;
            gpsr.m = 3;
            gpsr.sendTimes = "";
            gpsr.receivers = "";
            @display("p=444.44,666.67;i=misc/node_vs");
        }
        node47: GpsrRouter {
            gpsr.isMalicious = true;
            gpsr.isStation = false;
            gpsr.n = 1000;
            gpsr.m = 3;
            gpsr.sendTimes = "32.7s 34.3s";
            gpsr.receivers = "node55 node84";
            @display("p=444.44,777.78;i=misc/node_vs");
        }
        node48: GpsrRouter {
            gpsr.isMalicious = true;
            gpsr.isStation = false;
            gpsr.n = 1000;
            gpsr.m = 3;
            gpsr.sendTimes = "";
            gpsr.receivers = "";
            @display("p=444.44,888.89;i=misc/node_vs");
        }
        node49: GpsrRouter {
            gpsr.isMalicious = false;
            gpsr.isStation = false;
            gpsr.n = 1000;
            gpsr.m = 3;
            gpsr.sendTimes = "";
            gpsr.receivers = "";
            @display("p=444.44,1000.0;i=misc/node_vs");
        }
        node50: GpsrRouter {
            gpsr.isMalicious = false;
            gpsr.isStation = false;
            gpsr.n = 1000;
            gpsr.m = 3;
            gpsr.sendTimes = "35.5s";
            gpsr.receivers = "node80";
            @display("p=555.56,0.0;i=misc/node_vs");
        }
        node51: GpsrRouter {
            gpsr.isMalicious = false;
            gpsr.isStation = false;
            gpsr.n = 1000;
            gpsr.m = 3;
            gpsr.sendTimes = "36.6s 37.6s";
            gpsr.receivers = "node90 node78";
            @display("p=555.56,111.11;i=misc/node_vs");
        }
        node52: GpsrRouter {
            gpsr.isMalicious = false;
            gpsr.isStation = false;
            gpsr.n = 1000;
            gpsr.m = 3;
            gpsr.sendTimes = "37.9s 38.0s";
            gpsr.receivers = "node97 node82";
            @display("p=555.56,222.22;i=misc/node_vs");
        }
        node53: GpsrRouter {
            gpsr.isMalicious = false;
            gpsr.isStation = false;
            gpsr.n = 1000;
            gpsr.m = 3;
            gpsr.sendTimes = "38.1s";
            gpsr.receivers = "node95";
            @display("p=555.56,333.33;i=misc/node_vs");
        }
        node54: GpsrRouter {
            gpsr.isMalicious = false;
            gpsr.isStation = false;
            gpsr.n = 1000;
            gpsr.m = 3;
            gpsr.sendTimes = "38.3s 38.5s";
            gpsr.receivers = "node86 node60";
            @display("p=555.56,444.44;i=misc/node_vs");
        }
        node55: GpsrRouter {
            gpsr.isMalicious = false;
            gpsr.isStation = false;
            gpsr.n = 1000;
            gpsr.m = 3;
            gpsr.sendTimes = "";
            gpsr.receivers = "";
            @display("p=555.56,555.56;i=misc/node_vs");
        }
        node56: GpsrRouter {
            gpsr.isMalicious = false;
            gpsr.isStation = false;
            gpsr.n = 1000;
            gpsr.m = 3;
            gpsr.sendTimes = "";
            gpsr.receivers = "";
            @display("p=555.56,666.67;i=misc/node_vs");
        }
        node57: GpsrRouter {
            gpsr.isMalicious = false;
            gpsr.isStation = false;
            gpsr.n = 1000;
            gpsr.m = 3;
            gpsr.sendTimes = "";
            gpsr.receivers = "";
            @display("p=555.56,777.78;i=misc/node_vs");
        }
        node58: GpsrRouter {
            gpsr.isMalicious = false;
            gpsr.isStation = false;
            gpsr.n = 1000;
            gpsr.m = 3;
            gpsr.sendTimes = "38.6s 39.2s";
            gpsr.receivers = "node92 node71";
            @display("p=555.56,888.89;i=misc/node_vs");
        }
        node59: GpsrRouter {
            gpsr.isMalicious = true;
            gpsr.isStation = false;
            gpsr.n = 1000;
            gpsr.m = 3;
            gpsr.sendTimes = "39.7s";
            gpsr.receivers = "node88";
            @display("p=555.56,1000.0;i=misc/node_vs");
        }
        node60: GpsrRouter {
            gpsr.isMalicious = false;
            gpsr.isStation = false;
            gpsr.n = 1000;
            gpsr.m = 3;
            gpsr.sendTimes = "39.8s";
            gpsr.receivers = "node88";
            @display("p=666.67,0.0;i=misc/node_vs");
        }
        node61: GpsrRouter {
            gpsr.isMalicious = false;
            gpsr.isStation = false;
            gpsr.n = 1000;
            gpsr.m = 3;
            gpsr.sendTimes = "40.0s";
            gpsr.receivers = "node70";
            @display("p=666.67,111.11;i=misc/node_vs");
        }
        node62: GpsrRouter {
            gpsr.isMalicious = false;
            gpsr.isStation = false;
            gpsr.n = 1000;
            gpsr.m = 3;
            gpsr.sendTimes = "40.4s";
            gpsr.receivers = "node83";
            @display("p=666.67,222.22;i=misc/node_vs");
        }
        node63: GpsrRouter {
            gpsr.isMalicious = false;
            gpsr.isStation = false;
            gpsr.n = 1000;
            gpsr.m = 3;
            gpsr.sendTimes = "";
            gpsr.receivers = "";
            @display("p=666.67,333.33;i=misc/node_vs");
        }
        node64: GpsrRouter {
            gpsr.isMalicious = false;
            gpsr.isStation = false;
            gpsr.n = 1000;
            gpsr.m = 3;
            gpsr.sendTimes = "41.1s";
            gpsr.receivers = "node99";
            @display("p=666.67,444.44;i=misc/node_vs");
        }
        node65: GpsrRouter {
            gpsr.isMalicious = false;
            gpsr.isStation = false;
            gpsr.n = 1000;
            gpsr.m = 3;
            gpsr.sendTimes = "";
            gpsr.receivers = "";
            @display("p=666.67,555.56;i=misc/node_vs");
        }
        node66: GpsrRouter {
            gpsr.isMalicious = false;
            gpsr.isStation = false;
            gpsr.n = 1000;
            gpsr.m = 3;
            gpsr.sendTimes = "41.9s 42.1s";
            gpsr.receivers = "node80 node76";
            @display("p=666.67,666.67;i=misc/node_vs");
        }
        node67: GpsrRouter {
            gpsr.isMalicious = false;
            gpsr.isStation = false;
            gpsr.n = 1000;
            gpsr.m = 3;
            gpsr.sendTimes = "42.2s";
            gpsr.receivers = "node97";
            @display("p=666.67,777.78;i=misc/node_vs");
        }
        node68: GpsrRouter {
            gpsr.isMalicious = false;
            gpsr.isStation = false;
            gpsr.n = 1000;
            gpsr.m = 3;
            gpsr.sendTimes = "42.3s";
            gpsr.receivers = "node78";
            @display("p=666.67,888.89;i=misc/node_vs");
        }
        node69: GpsrRouter {
            gpsr.isMalicious = true;
            gpsr.isStation = false;
            gpsr.n = 1000;
            gpsr.m = 3;
            gpsr.sendTimes = "42.5s";
            gpsr.receivers = "node85";
            @display("p=666.67,1000.0;i=misc/node_vs");
        }
        node70: GpsrRouter {
            gpsr.isMalicious = false;
            gpsr.isStation = false;
            gpsr.n = 1000;
            gpsr.m = 3;
            gpsr.sendTimes = "";
            gpsr.receivers = "";
            @display("p=777.78,0.0;i=misc/node_vs");
        }
        node71: GpsrRouter {
            gpsr.isMalicious = false;
            gpsr.isStation = false;
            gpsr.n = 1000;
            gpsr.m = 3;
            gpsr.sendTimes = "";
            gpsr.receivers = "";
            @display("p=777.78,111.11;i=misc/node_vs");
        }
        node72: GpsrRouter {
            gpsr.isMalicious = false;
            gpsr.isStation = false;
            gpsr.n = 1000;
            gpsr.m = 3;
            gpsr.sendTimes = "42.6s";
            gpsr.receivers = "node87";
            @display("p=777.78,222.22;i=misc/node_vs");
        }
        node73: GpsrRouter {
            gpsr.isMalicious = false;
            gpsr.isStation = false;
            gpsr.n = 1000;
            gpsr.m = 3;
            gpsr.sendTimes = "42.8s";
            gpsr.receivers = "node96";
            @display("p=777.78,333.33;i=misc/node_vs");
        }
        node74: GpsrRouter {
            gpsr.isMalicious = true;
            gpsr.isStation = false;
            gpsr.n = 1000;
            gpsr.m = 3;
            gpsr.sendTimes = "43.3s";
            gpsr.receivers = "node89";
            @display("p=777.78,444.44;i=misc/node_vs");
        }
        node75: GpsrRouter {
            gpsr.isMalicious = false;
            gpsr.isStation = false;
            gpsr.n = 1000;
            gpsr.m = 3;
            gpsr.sendTimes = "";
            gpsr.receivers = "";
            @display("p=777.78,555.56;i=misc/node_vs");
        }
        node76: GpsrRouter {
            gpsr.isMalicious = false;
            gpsr.isStation = false;
            gpsr.n = 1000;
            gpsr.m = 3;
            gpsr.sendTimes = "43.4s";
            gpsr.receivers = "node80";
            @display("p=777.78,666.67;i=misc/node_vs");
        }
        node77: GpsrRouter {
            gpsr.isMalicious = false;
            gpsr.isStation = false;
            gpsr.n = 1000;
            gpsr.m = 3;
            gpsr.sendTimes = "";
            gpsr.receivers = "";
            @display("p=777.78,777.78;i=misc/node_vs");
        }
        node78: GpsrRouter {
            gpsr.isMalicious = false;
            gpsr.isStation = false;
            gpsr.n = 1000;
            gpsr.m = 3;
            gpsr.sendTimes = "";
            gpsr.receivers = "";
            @display("p=777.78,888.89;i=misc/node_vs");
        }
        node79: GpsrRouter {
            gpsr.isMalicious = false;
            gpsr.isStation = false;
            gpsr.n = 1000;
            gpsr.m = 3;
            gpsr.sendTimes = "";
            gpsr.receivers = "";
            @display("p=777.78,1000.0;i=misc/node_vs");
        }
        node80: GpsrRouter {
            gpsr.isMalicious = true;
            gpsr.isStation = false;
            gpsr.n = 1000;
            gpsr.m = 3;
            gpsr.sendTimes = "";
            gpsr.receivers = "";
            @display("p=888.89,0.0;i=misc/node_vs");
        }
        node81: GpsrRouter {
            gpsr.isMalicious = false;
            gpsr.isStation = false;
            gpsr.n = 1000;
            gpsr.m = 3;
            gpsr.sendTimes = "";
            gpsr.receivers = "";
            @display("p=888.89,111.11;i=misc/node_vs");
        }
        node82: GpsrRouter {
            gpsr.isMalicious = false;
            gpsr.isStation = false;
            gpsr.n = 1000;
            gpsr.m = 3;
            gpsr.sendTimes = "";
            gpsr.receivers = "";
            @display("p=888.89,222.22;i=misc/node_vs");
        }
        node83: GpsrRouter {
            gpsr.isMalicious = false;
            gpsr.isStation = false;
            gpsr.n = 1000;
            gpsr.m = 3;
            gpsr.sendTimes = "";
            gpsr.receivers = "";
            @display("p=888.89,333.33;i=misc/node_vs");
        }
        node84: GpsrRouter {
            gpsr.isMalicious = false;
            gpsr.isStation = false;
            gpsr.n = 1000;
            gpsr.m = 3;
            gpsr.sendTimes = "";
            gpsr.receivers = "";
            @display("p=888.89,444.44;i=misc/node_vs");
        }
        node85: GpsrRouter {
            gpsr.isMalicious = true;
            gpsr.isStation = false;
            gpsr.n = 1000;
            gpsr.m = 3;
            gpsr.sendTimes = "43.6s";
            gpsr.receivers = "node97";
            @display("p=888.89,555.56;i=misc/node_vs");
        }
        node86: GpsrRouter {
            gpsr.isMalicious = false;
            gpsr.isStation = false;
            gpsr.n = 1000;
            gpsr.m = 3;
            gpsr.sendTimes = "";
            gpsr.receivers = "";
            @display("p=888.89,666.67;i=misc/node_vs");
        }
        node87: GpsrRouter {
            gpsr.isMalicious = true;
            gpsr.isStation = false;
            gpsr.n = 1000;
            gpsr.m = 3;
            gpsr.sendTimes = "44.1s 44.3s";
            gpsr.receivers = "node98 node99";
            @display("p=888.89,777.78;i=misc/node_vs");
        }
        node88: GpsrRouter {
            gpsr.isMalicious = false;
            gpsr.isStation = false;
            gpsr.n = 1000;
            gpsr.m = 3;
            gpsr.sendTimes = "44.5s";
            gpsr.receivers = "node89";
            @display("p=888.89,888.89;i=misc/node_vs");
        }
        node89: GpsrRouter {
            gpsr.isMalicious = false;
            gpsr.isStation = false;
            gpsr.n = 1000;
            gpsr.m = 3;
            gpsr.sendTimes = "";
            gpsr.receivers = "";
            @display("p=888.89,1000.0;i=misc/node_vs");
        }
        node90: GpsrRouter {
            gpsr.isMalicious = false;
            gpsr.isStation = false;
            gpsr.n = 1000;
            gpsr.m = 3;
            gpsr.sendTimes = "44.9s";
            gpsr.receivers = "node98";
            @display("p=1000.0,0.0;i=misc/node_vs");
        }
        node91: GpsrRouter {
            gpsr.isMalicious = true;
            gpsr.isStation = false;
            gpsr.n = 1000;
            gpsr.m = 3;
            gpsr.sendTimes = "";
            gpsr.receivers = "";
            @display("p=1000.0,111.11;i=misc/node_vs");
        }
        node92: GpsrRouter {
            gpsr.isMalicious = false;
            gpsr.isStation = false;
            gpsr.n = 1000;
            gpsr.m = 3;
            gpsr.sendTimes = "";
            gpsr.receivers = "";
            @display("p=1000.0,222.22;i=misc/node_vs");
        }
        node93: GpsrRouter {
            gpsr.isMalicious = false;
            gpsr.isStation = false;
            gpsr.n = 1000;
            gpsr.m = 3;
            gpsr.sendTimes = "";
            gpsr.receivers = "";
            @display("p=1000.0,333.33;i=misc/node_vs");
        }
        node94: GpsrRouter {
            gpsr.isMalicious = true;
            gpsr.isStation = false;
            gpsr.n = 1000;
            gpsr.m = 3;
            gpsr.sendTimes = "";
            gpsr.receivers = "";
            @display("p=1000.0,444.44;i=misc/node_vs");
        }
        node95: GpsrRouter {
            gpsr.isMalicious = false;
            gpsr.isStation = false;
            gpsr.n = 1000;
            gpsr.m = 3;
            gpsr.sendTimes = "";
            gpsr.receivers = "";
            @display("p=1000.0,555.56;i=misc/node_vs");
        }
        node96: GpsrRouter {
            gpsr.isMalicious = false;
            gpsr.isStation = false;
            gpsr.n = 1000;
            gpsr.m = 3;
            gpsr.sendTimes = "";
            gpsr.receivers = "";
            @display("p=1000.0,666.67;i=misc/node_vs");
        }
        node97: GpsrRouter {
            gpsr.isMalicious = false;
            gpsr.isStation = false;
            gpsr.n = 1000;
            gpsr.m = 3;
            gpsr.sendTimes = "";
            gpsr.receivers = "";
            @display("p=1000.0,777.78;i=misc/node_vs");
        }
        node98: GpsrRouter {
            gpsr.isMalicious = false;
            gpsr.isStation = false;
            gpsr.n = 1000;
            gpsr.m = 3;
            gpsr.sendTimes = "";
            gpsr.receivers = "";
            @display("p=1000.0,888.89;i=misc/node_vs");
        }
        node99: GpsrRouter {
            gpsr.isMalicious = false;
            gpsr.isStation = false;
            gpsr.n = 1000;
            gpsr.m = 3;
            gpsr.sendTimes = "";
            gpsr.receivers = "";
            @display("p=1000.0,1000.0;i=misc/node_vs");
        }
        station0: GpsrRouter {
            gpsr.isMalicious = false;
            gpsr.isStation = true;
            gpsr.n = 1000;
            gpsr.m = 3;
            @display("p=166.66666666666666,166.66666666666666;i=misc/signal");
        }
        station3: GpsrRouter {
            gpsr.isMalicious = false;
            gpsr.isStation = true;
            gpsr.n = 1000;
            gpsr.m = 3;
            @display("p=166.66666666666666,500.0;i=misc/signal");
        }
        station6: GpsrRouter {
            gpsr.isMalicious = false;
            gpsr.isStation = true;
            gpsr.n = 1000;
            gpsr.m = 3;
            @display("p=166.66666666666666,833.3333333333333;i=misc/signal");
        }
        station1: GpsrRouter {
            gpsr.isMalicious = false;
            gpsr.isStation = true;
            gpsr.n = 1000;
            gpsr.m = 3;
            @display("p=500.0,166.66666666666666;i=misc/signal");
        }
        station4: GpsrRouter {
            gpsr.isMalicious = false;
            gpsr.isStation = true;
            gpsr.n = 1000;
            gpsr.m = 3;
            @display("p=500.0,500.0;i=misc/signal");
        }
        station7: GpsrRouter {
            gpsr.isMalicious = false;
            gpsr.isStation = true;
            gpsr.n = 1000;
            gpsr.m = 3;
            @display("p=500.0,833.3333333333333;i=misc/signal");
        }
        station2: GpsrRouter {
            gpsr.isMalicious = false;
            gpsr.isStation = true;
            gpsr.n = 1000;
            gpsr.m = 3;
            @display("p=833.3333333333333,166.66666666666666;i=misc/signal");
        }
        station5: GpsrRouter {
            gpsr.isMalicious = false;
            gpsr.isStation = true;
            gpsr.n = 1000;
            gpsr.m = 3;
            @display("p=833.3333333333333,500.0;i=misc/signal");
        }
        station8: GpsrRouter {
            gpsr.isMalicious = false;
            gpsr.isStation = true;
            gpsr.n = 1000;
            gpsr.m = 3;
            @display("p=833.3333333333333,833.3333333333333;i=misc/signal");
        }
}
