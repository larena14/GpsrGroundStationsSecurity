# 📡 GPSR with ground stations: attacks and mitigations

> Simulation and analysis of secure extensions to the GPSR protocol in high-mobility environments, with a focus on location spoofing attacks and mitigations using ground stations.

## 👥 Authors

- [@ringnelson](https://github.com/ringnelson)
- [@larena14](https://github.com/larena14)
- [@VincenzoNasso](https://github.com/VincenzoNasso)
- [@ginmy00](https://github.com/ginmy00)

---

## 📖 Overview

This project extends the **GPSR (Greedy Perimeter Stateless Routing)** protocol by introducing a hybrid infrastructure with **ground stations** in order to **detect and mitigate attacks** on the network.

Main features include:

- Collaborative management of node positions via fixed ground stations
- Message authenticity verification (ECDSA/EDDSA)
- Distributed trust evaluation system (`trustness`)
- Advanced mechanisms for malicious behavior detection

> For more details, see `Docs/Relazione_Aloia_Arena_Nasso_Saporito.pdf`

---

## 🧠 Attack Scenarios and Mitigations

| Scenario           | Description                                                                 | Folder                      |
|--------------------|-----------------------------------------------------------------------------|-----------------------------|
| **Attack 1**        | Spoofing location in beacon messages                                        | `Attack1Scenario/`          |
| **Mitigation 1**    | Beacon signature                                                            | `Mitigation1Scenario/`      |
| **Attack 2**        | Spoofing location to both nodes and stations                                | `Attack2Scenario/`          |
| **Mitigation 2**    | Location verification based on estimated distance (RSSI)                   | `Mitigation2Scenario/`      |
| **Mitigation 2+**   | Location verification via distance (RSSI) + triangulation with anchors      | `Mitigation2PosScenario/`   |
| **Attack 3**        | Unreliable or compromised stations                                          | `Attack3Scenario/`          |
| **Mitigation 3**    | Signed beacons and inter-station communication checks                       | `Mitigation3Scenario/`      |

---

## 🗂️ Repository Structure

```plaintext
.
├── Attack1Scenario/           # Attack 1: beacon position spoofing
├── Attack2Scenario/           # Attack 2: spoofing to both nodes and stations
├── Attack3Scenario/           # Attack 3: compromised or unreliable stations
│
├── Mitigation1Scenario/       # Mitigation for Attack 1: digital signatures
├── Mitigation2Scenario/       # Mitigation for Attack 2: distance checks (no triangulation)
├── Mitigation2PosScenario/    # Mitigation for Attack 2: distance + triangulation (Gauss-Newton)
├── Mitigation3Scenario/       # Mitigation for Attack 3: inter-station verification
│
├── Docs/                      # PDF Documentation
│   ├── Relazione_Aloia_Arena_Nasso_Saporito.pdf
│   └── Report_tecnico_Aloia_Arena_Nasso_Saporito.pdf
```

---

## ⚙️ Requirements

Recommended operating system: **Ubuntu 22.04 LTS**

### Dependencies

- [OMNeT++ 6.0.3](https://omnetpp.org/)
- INET Framework 4.5
- Eigen 3.4.0
- Libsodium 1.0.18 (for EDDSA support)
- OpenSSL 3.0.2 (for ECDSA support)

> For detailed installation instructions, please refer to `Docs/Report_tecnico_Aloia_Arena_Nasso_Saporito.pdf`

---

## 🧪 Running the Simulations

### 🔹 Using OMNeT++

1. Open **OMNeT++**
2. Import one of the scenario folders (e.g., `Attack1Scenario/`)
3. Build the project
4. Launch the simulation by running the `omnetpp.ini` configuration file

### 🔹 From the Terminal

Before running the simulations, make sure the project is compiled.

```bash
cd directory_attack_or_mitigation/
./launch_simulations.sh --areaMin <min_area_size> \
  --areaMax <max_area_size> --malMin <min_malicious_nodes> \
  --malMax <max_malicious_nodes> --nQuad <number_of_quadrants> \
  --nSim <number_of_simulations> --gridSize <sqrt_num_drones> \
  --nMsg <messages_to_generate> --simTime <simulation_duration> \
  --minInterval <min_msg_interval> --maxInterval <max_msg_interval> \
  --sendMode <sending_mode> --signAlgorithm <signature_algorithm>
```
> For more details, please refer to `Docs/Report_tecnico_Aloia_Arena_Nasso_Saporito.pdf`

---

## 📊 Metrics Collected
Each simulation produces a JSON file that contains:

- Total time spent in Greedy and Perimeter routing modes
- Total number of messages sent and delivered
- Average number of hops
- Energy consumption per node
- Undelivered messages

These statistics can be used to evaluate the impact of attacks and the effectiveness of mitigation mechanisms.
