package scenario1000;

import inet.common.scenario.ScenarioManager;
import inet.networklayer.configurator.ipv4.Ipv4NetworkConfigurator;
import inet.networklayer.ipv4.RoutingTableRecorder;
import GpsrRouter;
import inet.physicallayer.wireless.ieee80211.packetlevel.Ieee80211ScalarRadioMedium;
import inet.visualizer.common.IntegratedVisualizer;

network DronesNetwork1000_10
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
            gpsr.sendTimes = "5.4s 5.5s 5.6s";
            gpsr.receivers = "node9 node2 node67";
            @display("p=0.0,0.0;i=misc/node_vs");
        }
        node1: GpsrRouter {
            gpsr.isMalicious = false;
            gpsr.isStation = false;
            gpsr.n = 1000;
            gpsr.m = 3;
            gpsr.sendTimes = "6.0s 6.5s";
            gpsr.receivers = "node83 node89";
            @display("p=0.0,111.11;i=misc/node_vs");
        }
        node2: GpsrRouter {
            gpsr.isMalicious = false;
            gpsr.isStation = false;
            gpsr.n = 1000;
            gpsr.m = 3;
            gpsr.sendTimes = "7.4s 7.6s 8.8s";
            gpsr.receivers = "node59 node30 node82";
            @display("p=0.0,222.22;i=misc/node_vs");
        }
        node3: GpsrRouter {
            gpsr.isMalicious = false;
            gpsr.isStation = false;
            gpsr.n = 1000;
            gpsr.m = 3;
            gpsr.sendTimes = "9.8s";
            gpsr.receivers = "node68";
            @display("p=0.0,333.33;i=misc/node_vs");
        }
        node4: GpsrRouter {
            gpsr.isMalicious = false;
            gpsr.isStation = false;
            gpsr.n = 1000;
            gpsr.m = 3;
            gpsr.sendTimes = "11.1s 11.4s";
            gpsr.receivers = "node36 node11";
            @display("p=0.0,444.44;i=misc/node_vs");
        }
        node5: GpsrRouter {
            gpsr.isMalicious = false;
            gpsr.isStation = false;
            gpsr.n = 1000;
            gpsr.m = 3;
            gpsr.sendTimes = "11.5s";
            gpsr.receivers = "node37";
            @display("p=0.0,555.56;i=misc/node_vs");
        }
        node6: GpsrRouter {
            gpsr.isMalicious = false;
            gpsr.isStation = false;
            gpsr.n = 1000;
            gpsr.m = 3;
            gpsr.sendTimes = "12.5s 12.7s";
            gpsr.receivers = "node33 node87";
            @display("p=0.0,666.67;i=misc/node_vs");
        }
        node7: GpsrRouter {
            gpsr.isMalicious = false;
            gpsr.isStation = false;
            gpsr.n = 1000;
            gpsr.m = 3;
            gpsr.sendTimes = "12.9s";
            gpsr.receivers = "node91";
            @display("p=0.0,777.78;i=misc/node_vs");
        }
        node8: GpsrRouter {
            gpsr.isMalicious = false;
            gpsr.isStation = false;
            gpsr.n = 1000;
            gpsr.m = 3;
            gpsr.sendTimes = "13.3s";
            gpsr.receivers = "node12";
            @display("p=0.0,888.89;i=misc/node_vs");
        }
        node9: GpsrRouter {
            gpsr.isMalicious = false;
            gpsr.isStation = false;
            gpsr.n = 1000;
            gpsr.m = 3;
            gpsr.sendTimes = "13.8s 13.9s";
            gpsr.receivers = "node40 node70";
            @display("p=0.0,1000.0;i=misc/node_vs");
        }
        node10: GpsrRouter {
            gpsr.isMalicious = false;
            gpsr.isStation = false;
            gpsr.n = 1000;
            gpsr.m = 3;
            gpsr.sendTimes = "14.0s 14.9s 15.1s";
            gpsr.receivers = "node54 node12 node84";
            @display("p=111.11,0.0;i=misc/node_vs");
        }
        node11: GpsrRouter {
            gpsr.isMalicious = false;
            gpsr.isStation = false;
            gpsr.n = 1000;
            gpsr.m = 3;
            gpsr.sendTimes = "15.5s 15.7s 16.0s";
            gpsr.receivers = "node87 node82 node65";
            @display("p=111.11,111.11;i=misc/node_vs");
        }
        node12: GpsrRouter {
            gpsr.isMalicious = false;
            gpsr.isStation = false;
            gpsr.n = 1000;
            gpsr.m = 3;
            gpsr.sendTimes = "16.7s 17.7s 18.0s";
            gpsr.receivers = "node93 node69 node53";
            @display("p=111.11,222.22;i=misc/node_vs");
        }
        node13: GpsrRouter {
            gpsr.isMalicious = false;
            gpsr.isStation = false;
            gpsr.n = 1000;
            gpsr.m = 3;
            gpsr.sendTimes = "18.7s";
            gpsr.receivers = "node67";
            @display("p=111.11,333.33;i=misc/node_vs");
        }
        node14: GpsrRouter {
            gpsr.isMalicious = false;
            gpsr.isStation = false;
            gpsr.n = 1000;
            gpsr.m = 3;
            gpsr.sendTimes = "18.9s";
            gpsr.receivers = "node89";
            @display("p=111.11,444.44;i=misc/node_vs");
        }
        node15: GpsrRouter {
            gpsr.isMalicious = false;
            gpsr.isStation = false;
            gpsr.n = 1000;
            gpsr.m = 3;
            gpsr.sendTimes = "19.1s 19.2s";
            gpsr.receivers = "node30 node51";
            @display("p=111.11,555.56;i=misc/node_vs");
        }
        node16: GpsrRouter {
            gpsr.isMalicious = false;
            gpsr.isStation = false;
            gpsr.n = 1000;
            gpsr.m = 3;
            gpsr.sendTimes = "";
            gpsr.receivers = "";
            @display("p=111.11,666.67;i=misc/node_vs");
        }
        node17: GpsrRouter {
            gpsr.isMalicious = false;
            gpsr.isStation = false;
            gpsr.n = 1000;
            gpsr.m = 3;
            gpsr.sendTimes = "19.4s";
            gpsr.receivers = "node38";
            @display("p=111.11,777.78;i=misc/node_vs");
        }
        node18: GpsrRouter {
            gpsr.isMalicious = false;
            gpsr.isStation = false;
            gpsr.n = 1000;
            gpsr.m = 3;
            gpsr.sendTimes = "";
            gpsr.receivers = "";
            @display("p=111.11,888.89;i=misc/node_vs");
        }
        node19: GpsrRouter {
            gpsr.isMalicious = false;
            gpsr.isStation = false;
            gpsr.n = 1000;
            gpsr.m = 3;
            gpsr.sendTimes = "";
            gpsr.receivers = "";
            @display("p=111.11,1000.0;i=misc/node_vs");
        }
        node20: GpsrRouter {
            gpsr.isMalicious = true;
            gpsr.isStation = false;
            gpsr.n = 1000;
            gpsr.m = 3;
            gpsr.sendTimes = "";
            gpsr.receivers = "";
            @display("p=222.22,0.0;i=misc/node_vs");
        }
        node21: GpsrRouter {
            gpsr.isMalicious = false;
            gpsr.isStation = false;
            gpsr.n = 1000;
            gpsr.m = 3;
            gpsr.sendTimes = "19.8s 20.3s";
            gpsr.receivers = "node41 node34";
            @display("p=222.22,111.11;i=misc/node_vs");
        }
        node22: GpsrRouter {
            gpsr.isMalicious = false;
            gpsr.isStation = false;
            gpsr.n = 1000;
            gpsr.m = 3;
            gpsr.sendTimes = "21.6s";
            gpsr.receivers = "node96";
            @display("p=222.22,222.22;i=misc/node_vs");
        }
        node23: GpsrRouter {
            gpsr.isMalicious = false;
            gpsr.isStation = false;
            gpsr.n = 1000;
            gpsr.m = 3;
            gpsr.sendTimes = "22.6s";
            gpsr.receivers = "node63";
            @display("p=222.22,333.33;i=misc/node_vs");
        }
        node24: GpsrRouter {
            gpsr.isMalicious = false;
            gpsr.isStation = false;
            gpsr.n = 1000;
            gpsr.m = 3;
            gpsr.sendTimes = "22.7s 23.1s";
            gpsr.receivers = "node64 node29";
            @display("p=222.22,444.44;i=misc/node_vs");
        }
        node25: GpsrRouter {
            gpsr.isMalicious = false;
            gpsr.isStation = false;
            gpsr.n = 1000;
            gpsr.m = 3;
            gpsr.sendTimes = "23.6s 24.1s";
            gpsr.receivers = "node97 node64";
            @display("p=222.22,555.56;i=misc/node_vs");
        }
        node26: GpsrRouter {
            gpsr.isMalicious = false;
            gpsr.isStation = false;
            gpsr.n = 1000;
            gpsr.m = 3;
            gpsr.sendTimes = "24.7s";
            gpsr.receivers = "node31";
            @display("p=222.22,666.67;i=misc/node_vs");
        }
        node27: GpsrRouter {
            gpsr.isMalicious = false;
            gpsr.isStation = false;
            gpsr.n = 1000;
            gpsr.m = 3;
            gpsr.sendTimes = "25.6s";
            gpsr.receivers = "node84";
            @display("p=222.22,777.78;i=misc/node_vs");
        }
        node28: GpsrRouter {
            gpsr.isMalicious = false;
            gpsr.isStation = false;
            gpsr.n = 1000;
            gpsr.m = 3;
            gpsr.sendTimes = "25.9s 26.1s 26.3s";
            gpsr.receivers = "node69 node38 node97";
            @display("p=222.22,888.89;i=misc/node_vs");
        }
        node29: GpsrRouter {
            gpsr.isMalicious = false;
            gpsr.isStation = false;
            gpsr.n = 1000;
            gpsr.m = 3;
            gpsr.sendTimes = "26.7s";
            gpsr.receivers = "node84";
            @display("p=222.22,1000.0;i=misc/node_vs");
        }
        node30: GpsrRouter {
            gpsr.isMalicious = false;
            gpsr.isStation = false;
            gpsr.n = 1000;
            gpsr.m = 3;
            gpsr.sendTimes = "26.8s 27.2s 27.3s 27.9s";
            gpsr.receivers = "node54 node65 node69 node80";
            @display("p=333.33,0.0;i=misc/node_vs");
        }
        node31: GpsrRouter {
            gpsr.isMalicious = false;
            gpsr.isStation = false;
            gpsr.n = 1000;
            gpsr.m = 3;
            gpsr.sendTimes = "28.0s 28.4s";
            gpsr.receivers = "node56 node72";
            @display("p=333.33,111.11;i=misc/node_vs");
        }
        node32: GpsrRouter {
            gpsr.isMalicious = false;
            gpsr.isStation = false;
            gpsr.n = 1000;
            gpsr.m = 3;
            gpsr.sendTimes = "29.3s 29.5s";
            gpsr.receivers = "node84 node54";
            @display("p=333.33,222.22;i=misc/node_vs");
        }
        node33: GpsrRouter {
            gpsr.isMalicious = false;
            gpsr.isStation = false;
            gpsr.n = 1000;
            gpsr.m = 3;
            gpsr.sendTimes = "29.6s";
            gpsr.receivers = "node82";
            @display("p=333.33,333.33;i=misc/node_vs");
        }
        node34: GpsrRouter {
            gpsr.isMalicious = false;
            gpsr.isStation = false;
            gpsr.n = 1000;
            gpsr.m = 3;
            gpsr.sendTimes = "";
            gpsr.receivers = "";
            @display("p=333.33,444.44;i=misc/node_vs");
        }
        node35: GpsrRouter {
            gpsr.isMalicious = false;
            gpsr.isStation = false;
            gpsr.n = 1000;
            gpsr.m = 3;
            gpsr.sendTimes = "29.8s";
            gpsr.receivers = "node51";
            @display("p=333.33,555.56;i=misc/node_vs");
        }
        node36: GpsrRouter {
            gpsr.isMalicious = false;
            gpsr.isStation = false;
            gpsr.n = 1000;
            gpsr.m = 3;
            gpsr.sendTimes = "30.0s 30.2s";
            gpsr.receivers = "node42 node90";
            @display("p=333.33,666.67;i=misc/node_vs");
        }
        node37: GpsrRouter {
            gpsr.isMalicious = false;
            gpsr.isStation = false;
            gpsr.n = 1000;
            gpsr.m = 3;
            gpsr.sendTimes = "30.4s 31.3s 31.7s 32.6s";
            gpsr.receivers = "node88 node76 node41 node89";
            @display("p=333.33,777.78;i=misc/node_vs");
        }
        node38: GpsrRouter {
            gpsr.isMalicious = false;
            gpsr.isStation = false;
            gpsr.n = 1000;
            gpsr.m = 3;
            gpsr.sendTimes = "33.0s";
            gpsr.receivers = "node95";
            @display("p=333.33,888.89;i=misc/node_vs");
        }
        node39: GpsrRouter {
            gpsr.isMalicious = false;
            gpsr.isStation = false;
            gpsr.n = 1000;
            gpsr.m = 3;
            gpsr.sendTimes = "33.3s";
            gpsr.receivers = "node62";
            @display("p=333.33,1000.0;i=misc/node_vs");
        }
        node40: GpsrRouter {
            gpsr.isMalicious = false;
            gpsr.isStation = false;
            gpsr.n = 1000;
            gpsr.m = 3;
            gpsr.sendTimes = "33.4s";
            gpsr.receivers = "node85";
            @display("p=444.44,0.0;i=misc/node_vs");
        }
        node41: GpsrRouter {
            gpsr.isMalicious = false;
            gpsr.isStation = false;
            gpsr.n = 1000;
            gpsr.m = 3;
            gpsr.sendTimes = "33.6s 34.9s";
            gpsr.receivers = "node90 node53";
            @display("p=444.44,111.11;i=misc/node_vs");
        }
        node42: GpsrRouter {
            gpsr.isMalicious = false;
            gpsr.isStation = false;
            gpsr.n = 1000;
            gpsr.m = 3;
            gpsr.sendTimes = "";
            gpsr.receivers = "";
            @display("p=444.44,222.22;i=misc/node_vs");
        }
        node43: GpsrRouter {
            gpsr.isMalicious = false;
            gpsr.isStation = false;
            gpsr.n = 1000;
            gpsr.m = 3;
            gpsr.sendTimes = "35.2s";
            gpsr.receivers = "node50";
            @display("p=444.44,333.33;i=misc/node_vs");
        }
        node44: GpsrRouter {
            gpsr.isMalicious = false;
            gpsr.isStation = false;
            gpsr.n = 1000;
            gpsr.m = 3;
            gpsr.sendTimes = "35.6s";
            gpsr.receivers = "node81";
            @display("p=444.44,444.44;i=misc/node_vs");
        }
        node45: GpsrRouter {
            gpsr.isMalicious = true;
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
            gpsr.sendTimes = "35.9s";
            gpsr.receivers = "node96";
            @display("p=444.44,666.67;i=misc/node_vs");
        }
        node47: GpsrRouter {
            gpsr.isMalicious = false;
            gpsr.isStation = false;
            gpsr.n = 1000;
            gpsr.m = 3;
            gpsr.sendTimes = "";
            gpsr.receivers = "";
            @display("p=444.44,777.78;i=misc/node_vs");
        }
        node48: GpsrRouter {
            gpsr.isMalicious = false;
            gpsr.isStation = false;
            gpsr.n = 1000;
            gpsr.m = 3;
            gpsr.sendTimes = "36.1s 36.2s";
            gpsr.receivers = "node55 node66";
            @display("p=444.44,888.89;i=misc/node_vs");
        }
        node49: GpsrRouter {
            gpsr.isMalicious = false;
            gpsr.isStation = false;
            gpsr.n = 1000;
            gpsr.m = 3;
            gpsr.sendTimes = "36.4s";
            gpsr.receivers = "node55";
            @display("p=444.44,1000.0;i=misc/node_vs");
        }
        node50: GpsrRouter {
            gpsr.isMalicious = false;
            gpsr.isStation = false;
            gpsr.n = 1000;
            gpsr.m = 3;
            gpsr.sendTimes = "";
            gpsr.receivers = "";
            @display("p=555.56,0.0;i=misc/node_vs");
        }
        node51: GpsrRouter {
            gpsr.isMalicious = false;
            gpsr.isStation = false;
            gpsr.n = 1000;
            gpsr.m = 3;
            gpsr.sendTimes = "36.5s 36.6s 36.7s";
            gpsr.receivers = "node52 node63 node96";
            @display("p=555.56,111.11;i=misc/node_vs");
        }
        node52: GpsrRouter {
            gpsr.isMalicious = true;
            gpsr.isStation = false;
            gpsr.n = 1000;
            gpsr.m = 3;
            gpsr.sendTimes = "";
            gpsr.receivers = "";
            @display("p=555.56,222.22;i=misc/node_vs");
        }
        node53: GpsrRouter {
            gpsr.isMalicious = false;
            gpsr.isStation = false;
            gpsr.n = 1000;
            gpsr.m = 3;
            gpsr.sendTimes = "37.6s 38.2s 39.2s";
            gpsr.receivers = "node82 node75 node58";
            @display("p=555.56,333.33;i=misc/node_vs");
        }
        node54: GpsrRouter {
            gpsr.isMalicious = true;
            gpsr.isStation = false;
            gpsr.n = 1000;
            gpsr.m = 3;
            gpsr.sendTimes = "39.3s 39.6s";
            gpsr.receivers = "node90 node70";
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
            gpsr.sendTimes = "39.8s 39.9s";
            gpsr.receivers = "node91 node93";
            @display("p=555.56,777.78;i=misc/node_vs");
        }
        node58: GpsrRouter {
            gpsr.isMalicious = false;
            gpsr.isStation = false;
            gpsr.n = 1000;
            gpsr.m = 3;
            gpsr.sendTimes = "40.1s";
            gpsr.receivers = "node88";
            @display("p=555.56,888.89;i=misc/node_vs");
        }
        node59: GpsrRouter {
            gpsr.isMalicious = false;
            gpsr.isStation = false;
            gpsr.n = 1000;
            gpsr.m = 3;
            gpsr.sendTimes = "";
            gpsr.receivers = "";
            @display("p=555.56,1000.0;i=misc/node_vs");
        }
        node60: GpsrRouter {
            gpsr.isMalicious = false;
            gpsr.isStation = false;
            gpsr.n = 1000;
            gpsr.m = 3;
            gpsr.sendTimes = "";
            gpsr.receivers = "";
            @display("p=666.67,0.0;i=misc/node_vs");
        }
        node61: GpsrRouter {
            gpsr.isMalicious = false;
            gpsr.isStation = false;
            gpsr.n = 1000;
            gpsr.m = 3;
            gpsr.sendTimes = "";
            gpsr.receivers = "";
            @display("p=666.67,111.11;i=misc/node_vs");
        }
        node62: GpsrRouter {
            gpsr.isMalicious = false;
            gpsr.isStation = false;
            gpsr.n = 1000;
            gpsr.m = 3;
            gpsr.sendTimes = "40.4s";
            gpsr.receivers = "node77";
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
            gpsr.isMalicious = true;
            gpsr.isStation = false;
            gpsr.n = 1000;
            gpsr.m = 3;
            gpsr.sendTimes = "40.7s 41.1s 42.0s";
            gpsr.receivers = "node93 node81 node74";
            @display("p=666.67,444.44;i=misc/node_vs");
        }
        node65: GpsrRouter {
            gpsr.isMalicious = true;
            gpsr.isStation = false;
            gpsr.n = 1000;
            gpsr.m = 3;
            gpsr.sendTimes = "42.2s";
            gpsr.receivers = "node98";
            @display("p=666.67,555.56;i=misc/node_vs");
        }
        node66: GpsrRouter {
            gpsr.isMalicious = false;
            gpsr.isStation = false;
            gpsr.n = 1000;
            gpsr.m = 3;
            gpsr.sendTimes = "42.5s";
            gpsr.receivers = "node79";
            @display("p=666.67,666.67;i=misc/node_vs");
        }
        node67: GpsrRouter {
            gpsr.isMalicious = false;
            gpsr.isStation = false;
            gpsr.n = 1000;
            gpsr.m = 3;
            gpsr.sendTimes = "";
            gpsr.receivers = "";
            @display("p=666.67,777.78;i=misc/node_vs");
        }
        node68: GpsrRouter {
            gpsr.isMalicious = false;
            gpsr.isStation = false;
            gpsr.n = 1000;
            gpsr.m = 3;
            gpsr.sendTimes = "";
            gpsr.receivers = "";
            @display("p=666.67,888.89;i=misc/node_vs");
        }
        node69: GpsrRouter {
            gpsr.isMalicious = false;
            gpsr.isStation = false;
            gpsr.n = 1000;
            gpsr.m = 3;
            gpsr.sendTimes = "";
            gpsr.receivers = "";
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
            gpsr.sendTimes = "42.6s";
            gpsr.receivers = "node78";
            @display("p=777.78,111.11;i=misc/node_vs");
        }
        node72: GpsrRouter {
            gpsr.isMalicious = false;
            gpsr.isStation = false;
            gpsr.n = 1000;
            gpsr.m = 3;
            gpsr.sendTimes = "42.7s";
            gpsr.receivers = "node76";
            @display("p=777.78,222.22;i=misc/node_vs");
        }
        node73: GpsrRouter {
            gpsr.isMalicious = false;
            gpsr.isStation = false;
            gpsr.n = 1000;
            gpsr.m = 3;
            gpsr.sendTimes = "43.0s";
            gpsr.receivers = "node83";
            @display("p=777.78,333.33;i=misc/node_vs");
        }
        node74: GpsrRouter {
            gpsr.isMalicious = false;
            gpsr.isStation = false;
            gpsr.n = 1000;
            gpsr.m = 3;
            gpsr.sendTimes = "";
            gpsr.receivers = "";
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
            gpsr.sendTimes = "43.1s";
            gpsr.receivers = "node88";
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
            gpsr.sendTimes = "43.3s 43.8s";
            gpsr.receivers = "node94 node98";
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
            gpsr.isMalicious = false;
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
            gpsr.sendTimes = "44.0s";
            gpsr.receivers = "node89";
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
            gpsr.isMalicious = false;
            gpsr.isStation = false;
            gpsr.n = 1000;
            gpsr.m = 3;
            gpsr.sendTimes = "";
            gpsr.receivers = "";
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
            gpsr.isMalicious = false;
            gpsr.isStation = false;
            gpsr.n = 1000;
            gpsr.m = 3;
            gpsr.sendTimes = "44.1s";
            gpsr.receivers = "node88";
            @display("p=888.89,777.78;i=misc/node_vs");
        }
        node88: GpsrRouter {
            gpsr.isMalicious = false;
            gpsr.isStation = false;
            gpsr.n = 1000;
            gpsr.m = 3;
            gpsr.sendTimes = "";
            gpsr.receivers = "";
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
            gpsr.sendTimes = "";
            gpsr.receivers = "";
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
            gpsr.sendTimes = "44.2s";
            gpsr.receivers = "node96";
            @display("p=1000.0,222.22;i=misc/node_vs");
        }
        node93: GpsrRouter {
            gpsr.isMalicious = true;
            gpsr.isStation = false;
            gpsr.n = 1000;
            gpsr.m = 3;
            gpsr.sendTimes = "";
            gpsr.receivers = "";
            @display("p=1000.0,333.33;i=misc/node_vs");
        }
        node94: GpsrRouter {
            gpsr.isMalicious = false;
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
            gpsr.sendTimes = "45.0s";
            gpsr.receivers = "node99";
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
            gpsr.isMalicious = true;
            gpsr.isStation = false;
            gpsr.n = 1000;
            gpsr.m = 3;
            gpsr.sendTimes = "";
            gpsr.receivers = "";
            @display("p=1000.0,777.78;i=misc/node_vs");
        }
        node98: GpsrRouter {
            gpsr.isMalicious = true;
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
