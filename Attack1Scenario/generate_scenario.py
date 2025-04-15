#!/usr/bin/python3

import random
import os
import argparse
import math

def get_station_index(x, y, n, m):
    l = n / m
 
    qX = int(math.floor(x / l))
    qY = int(math.floor(y / l))
 
    if qX >= m:
        qX = m - 1
    if qY >= m:
        qY = m - 1
 
    q = qY * m + qX
 
    return q

def generate_list(N, length):
    upper_limit = (3/4) * N
    # Usa un set per garantire l'unicità dei numeri
    numbers = set()
    while len(numbers) < length:
        number = round(random.uniform(5, upper_limit), 1)
        numbers.add(number)
    # Ordina i numeri in ordine crescente
    sorted_numbers = sorted(numbers)
    return sorted_numbers

def compute_station_position(area_size, M):
    
    quadrant_size = area_size / M
    quadrants_centers = []

    for i in range(M):
        for j in range(M):
            center_x = (i + 0.5) * quadrant_size
            center_y = (j + 0.5) * quadrant_size
            quadrants_centers.append((center_x, center_y))
    
    return quadrants_centers

def create_ini_file(path_ini, sim_time, minInterval, maxInterval, sendMode):
    ini_file = open(path_ini, "w")
    content = f"""
[General]
sim-time-limit = {sim_time}s

# channel physical parameters
**.wlan[*].radio.typename = "Ieee80211ScalarRadio"
**.wlan[*].bitrate = 24Mbps
**.wlan[*].mac.dcf.channelAccess.cwMin = 15
**.wlan[*].radio.transmitter.power = 2mW
**.beaconInterval = 1.0s
**.planarizationMode = "GG"

*.node*.mobility.typename = "RandomWaypointMobility"
*.node*.mobility.speed = uniform(0mps, 10mps)
*.node*.mobility.waitTime = 2s

**.statistic-recording = false
**.scalar-recording = false
**.vector-recording = false

**.sendMode = "{sendMode}"
**.minInterval = {minInterval}s
**.maxInterval = {maxInterval}s
*.node*.hasStatus = true

#power
*.node*.energyStorage.typename = "SimpleEpEnergyStorage"
*.node*.energyStorage.nominalCapacity =  55.5Wh
*.node*.energyStorage.initialCapacity =  uniform(47.175Wh, this.nominalCapacity)
*.node*.wlan[*].radio.energyConsumer.typename = "CustomStateBasedEpEnergyConsumer"
*.node*.energyConsumer.signPowerConsumption = 2mW #ECDSA 2mW; EDDSA 1.5 mW
*.node*.energyConsumer.verifyPowerConsumption = 1.5mW #ECDSA 1.5mW; EDDSA 1 mW
*.node*.energyConsumer.hoveringPowerConsumption = 150W
*.node*.energyConsumer.maxMovingPowerConsumption = 250W
    """
    ini_file.write(content)
    ini_file.close()

def create_ned_file(grid_size, area_size, malicious_percentage, output_dir, package_name, path_ini, num_msg, m, sim_time, minInterval, maxInterval, sendMode):
    distance = area_size / (grid_size - 1)
    output_file = f"{output_dir}/DronesNetwork{area_size}_{int(malicious_percentage*100)}.ned"

    total_nodes = grid_size * grid_size
    num_malicious_nodes = int(malicious_percentage * total_nodes)

    malicious_nodes = random.sample(list(range(total_nodes)), num_malicious_nodes)
    
    good_nodes = [node for node in list(range(total_nodes)) if node not in malicious_nodes]
    
    
    
    pairs = set()
    while len(pairs) < num_msg:
        pair = tuple(sorted(random.sample(list(range(total_nodes)), 2)))  
        pairs.add(pair)
        
    
    
    times = generate_list(sim_time, num_msg)
    
            
    ned_file = open(output_file, "w")

    ned_file.write(f"package {package_name};\n")
    ned_file.write("\n")
    ned_file.write("import inet.common.scenario.ScenarioManager;\n")
    ned_file.write("import inet.networklayer.configurator.ipv4.Ipv4NetworkConfigurator;\n")
    ned_file.write("import inet.networklayer.ipv4.RoutingTableRecorder;\n")
    ned_file.write("import GpsrRouter;\n")
    ned_file.write("import inet.physicallayer.wireless.ieee80211.packetlevel.Ieee80211ScalarRadioMedium;\n")
    ned_file.write("import inet.visualizer.common.IntegratedVisualizer;\n\n")

    ned_file.write(f"network DronesNetwork{area_size}_{int(malicious_percentage * 100)}\n")
    ned_file.write("{\n")
    ned_file.write("    parameters:\n")
    ned_file.write(f"        @display(\"bgb={area_size},{area_size}\");\n")
    ned_file.write("    submodules:\n")
    ned_file.write("        radioMedium: Ieee80211ScalarRadioMedium {\n")
    ned_file.write("            parameters:\n")
    ned_file.write(f"                @display(\"p={area_size + 200},211;is=s\");\n")
    ned_file.write("        }\n")
    ned_file.write("        visualizer: IntegratedVisualizer {\n")
    ned_file.write(f"            @display(\"p={area_size + 200},512\");\n")
    ned_file.write("        }\n")
    ned_file.write("        configurator: Ipv4NetworkConfigurator {\n")
    ned_file.write("            parameters:\n")
    ned_file.write(
        "                config = xml(\"<config><interface hosts='*' address='145.236.x.x' netmask='255.255.0.0'/></config>\");\n")
    ned_file.write(f"                @display(\"p={area_size + 200},110;is=s\");\n")
    ned_file.write("        }\n")
    ned_file.write("        routingTableRecorder: RoutingTableRecorder {\n")
    ned_file.write(f"            @display(\"p={area_size + 200},311\");\n")
    ned_file.write("        }\n")
    ned_file.write("        scenarioManager: ScenarioManager {\n")
    ned_file.write("            parameters:\n")
    ned_file.write(
        "                script = default(xml(\"<scenario/>\"));\n")
    ned_file.write(f"                @display(\"p={area_size + 200},412;is=s\");\n")
    ned_file.write("        }\n")

    # Configura i nodi con posizioni in griglia e i ruoli specificati
    times_index = 0
    for i in range(total_nodes):
        x = round((i // grid_size) * distance, 2)
        y = round((i % grid_size) * distance, 2)

        malicious_flag = "true" if i in malicious_nodes else "false"
        ned_file.write(f"        node{i}: GpsrRouter {{\n")
        ned_file.write(f"            gpsr.isMalicious = {malicious_flag};\n")
        ned_file.write(f"            gpsr.isStation = false;\n")
        ned_file.write(f"            gpsr.n = {area_size};\n")
        ned_file.write(f"            gpsr.m = {m};\n")
        
        
        rec = ""
        t = ""
        for pair in pairs:
            if pair[0] == i:
                rec += f"node{pair[1]} "
                t += f"{times[times_index]}s "
                times_index += 1
        
        rec = rec.strip()
        t = t.strip()
        ned_file.write(f"            gpsr.sendTimes = \"{t}\";\n")
        ned_file.write(f"            gpsr.receivers = \"{rec}\";\n")
        ned_file.write(f"            @display(\"p={x},{y};i=misc/node_vs\");\n")
        ned_file.write("        }\n")
        
    quadrants_centers = compute_station_position(area_size, m)
    
    for i, (x, y) in enumerate(quadrants_centers):
        index = get_station_index(x, y, area_size, m)
        ned_file.write(f"        station{index}: GpsrRouter {{\n")
        ned_file.write(f"            gpsr.isMalicious = false;\n")
        ned_file.write(f"            gpsr.isStation = true;\n")
        ned_file.write(f"            gpsr.n = {area_size};\n")
        ned_file.write(f"            gpsr.m = {m};\n")
        ned_file.write(f"            @display(\"p={x},{y};i=misc/signal\");\n")
        ned_file.write("        }\n")

    ned_file.write("}\n")
    
    

    print(f"File {output_file} generato con successo.")

    ned_file.close()
    
    ini_file = open(path_ini, "a")

    content = f"""
[Config Scenario{area_size}_{int(malicious_percentage*100)}]
network = {package_name}.DronesNetwork{area_size}_{int(malicious_percentage*100)}

# mobility
*.node*.mobility.constraintAreaMinX = 0m
*.node*.mobility.constraintAreaMinY = 0m
*.node*.mobility.constraintAreaMinZ = 100m
*.node*.mobility.constraintAreaMaxX = {area_size}m
*.node*.mobility.constraintAreaMaxY = {area_size}m
*.node*.mobility.constraintAreaMaxZ = 120m

*.station*.mobility.constraintAreaMinX = 0m
*.station*.mobility.constraintAreaMinY = 0m
*.station*.mobility.constraintAreaMinZ = 0m
*.station*.mobility.constraintAreaMaxX = {area_size}m
*.station*.mobility.constraintAreaMaxY = {area_size}m
*.station*.mobility.constraintAreaMaxZ = 0m


    """
    ini_file.write(content)
    ini_file.close()

    print("Aggiunte all'ini con successo")



def main():
    parser = argparse.ArgumentParser(description="Script setup scenari")
    parser.add_argument('--droneGridSize', type=int, required=True, help='Radice quadrata del numero di nodi (intero)')
    parser.add_argument('--areaMin', type=int, required=True, help='Area minima')
    parser.add_argument('--areaMax', type=int, required=True, help='Area massima')
    parser.add_argument('--malMin', type=int, required=True, help='Minima percentuale di nodi malevoli')
    parser.add_argument('--malMax', type=int, required=True, help='Massima percentuale di nodi malevoli')
    parser.add_argument('--nQuad', type=int, required=True, help='Numero di quadranti')
    parser.add_argument('--simTime', type=int, required=True, help='Durata simulazione')
    parser.add_argument('--nMsg', type=int, required=True, help='Numero messaggi')
    parser.add_argument('--minInterval', type=int, required=True, help='Intervallo minimo (Interval mode)')
    parser.add_argument('--maxInterval', type=int, required=True, help='Intervallo massimo (Interval mode)')
    parser.add_argument('--sendMode', type=str, required=True, help='Send mode (FixedMaximum o Interval)')
    args = parser.parse_args()

    grid_size = args.droneGridSize
    area_min = args.areaMin
    area_max = args.areaMax
    min_mal = int(args.malMin/10)
    max_mal = int(args.malMax/10)
    n_quad = args.nQuad
    m = math.sqrt(n_quad)
    if m.is_integer():
        m = int(m)
    else:
        print("nQuad must be perfect square")
        return
    sim_time = args.simTime
    n_msg = args.nMsg
    minInterval = args.minInterval
    maxInterval = args.maxInterval
    sendMode = args.sendMode

    create_ini_file("./omnetpp.ini", sim_time, minInterval, maxInterval, sendMode)

    for area in range(area_min, area_max+1, 1000):
        directory_name = f"scenario{area}"
        if not os.path.exists(directory_name):
            os.makedirs(directory_name)
        for j in range(min_mal, max_mal+1):
            mal_per = round(j * 0.1, 1)
            create_ned_file(grid_size, area, mal_per, f"./{directory_name}", directory_name, "./omnetpp.ini", n_msg, m, sim_time, minInterval, maxInterval, sendMode)


if __name__ == "__main__":
    main()




