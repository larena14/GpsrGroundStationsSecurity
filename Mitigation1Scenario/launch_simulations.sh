#!/bin/bash


if [ $# -lt 26 ]; then
    echo "Utilizzo: $0 --areaMin <valore> --areaMax <valore> --malMin <valore> --malMax <valore> --nQuad <valore> --nSim <valore> --gridSize <valore> --nMsg <valore> --simTime <valore> --minInterval <valore> --maxInterval <valore> --sendMode <valore> --signAlgorithm <valore>"
    exit 1
fi

areaMin=0
areaMax=0
malMin=0
malMax=0
nQuad=0
nSim=0
gridSize=0
nMsg=0
simTime=0
minInterval=0
maxInterval=0
sendMode=""
signAlgorithm=""

# Parsing dei parametri
while [[ $# -gt 0 ]]; do
    case "$1" in
        --areaMin)
            areaMin="$2"
            shift 2
            ;;
        --areaMax)
            areaMax="$2"
            shift 2
            ;;
        --malMin)
            malMin="$2"
            shift 2
            ;;
        --malMax)
            malMax="$2"
            shift 2
            ;;
        --nQuad)
            nQuad="$2"
            shift 2
            ;;
        --nSim)
            nSim="$2"
            shift 2
            ;;
        --gridSize)
            gridSize="$2"
            shift 2
            ;;
        --nMsg)
            nMsg="$2"
            shift 2
            ;;
        --simTime)
            simTime="$2"
            shift 2
            ;;
        --minInterval)
            minInterval="$2"
            shift 2
            ;;
        --maxInterval)
            maxInterval="$2"
            shift 2
            ;;
        --sendMode)
            sendMode="$2"
            shift 2
            ;;
        --signAlgorithm)
            signAlgorithm="$2"
            shift 2
            ;;
        *)
            echo "Parametro sconosciuto: $1"
            exit 1
            ;;
    esac
done

echo "Valore minimo dell'area: $areaMin"
echo "Valore massimo dell'area: $areaMax"
echo "Numero minimo di nodi malevoli: $malMin"
echo "Numero massimo di nodi malevoli: $malMax"
echo "Numero di quadranti: $nQuad"
echo "Numero di simulazioni: $nSim"
echo "Dimensione della griglia: $gridSize"
echo "Numero di messaggi: $nMsg"
echo "Tempo di simulazione: $simTime"
echo "Intervallo minimo di invio: $minInterval"
echo "Intervallo massimo di invio: $maxInterval"
echo "Modalità di invio: $sendMode"
echo "Algoritmo di firma: $signAlgorithm"


cd /home/user/omnetpp-6.0.3 && source setenv

cd /home/user/omnetpp-6.0.3/samples/Mitigation1Scenario


for ((j=areaMin; j<=areaMax; j+=1000)); do
	for ((i=1; i<=nSim; i++)); do
		echo "Iteration ${i}"

		echo "Generating scenario..."

		./generate_scenario.py --droneGridSize ${gridSize} --areaMin ${areaMin} --areaMax ${j} --malMin ${malMin} \
		        --malMax ${malMax} --nQuad ${nQuad} --simTime ${simTime} --nMsg ${nMsg} --minInterval ${minInterval} \
		        --maxInterval ${maxInterval} --sendMode ${sendMode} --signAlgorithm ${signAlgorithm}

		echo "Done."

		sleep 3

	    echo "Launching simulations..."

	    for ((p=malMin; p<=malMax; p+=10)); do 
            while true; do
               ./Mitigation1Scenario -r 0 -m -u Cmdenv -c Scenario${j}_${p} -n .:../inet4.5/examples:../inet4.5/showcases:../inet4.5/src:../inet4.5/tests/validation:../inet4.5/tests/networks:../inet4.5/tutorials -x inet.common.selfdoc\;inet.linklayer.configurator.gatescheduling.z3\;inet.emulation\;inet.examples.emulation\;inet.showcases.emulation\;inet.transportlayer.tcp_lwip\;inet.applications.voipstream\;inet.examples.voipstream --image-path=../inet4.5/images -l ../inet4.5/src/INET omnetpp.ini

                if [ $? -eq 0 ]; then
                    sleep 2 
                    mv results/json/output.json results/json/Scenario${j}_${p}_${i}_${sendMode}_${signAlgorithm}.json
                    break
                else
                    echo "Simulation failed for Scenario${j}_${p}. Retrying..."
                    sleep 2
                    ./generate_scenario.py --droneGridSize ${gridSize} --areaMin ${j} --areaMax ${j} --malMin ${p} --malMax ${p} \
                                           --nQuad ${nQuad} --simTime ${simTime} --nMsg ${nMsg} \
                                           --minInterval ${minInterval} \
		        						   --maxInterval ${maxInterval} --sendMode ${sendMode} --signAlgorithm ${signAlgorithm}
                    sleep 3
                fi
            done
	    done

		echo "Done."

		sleep 3

done

done

<<comment
Example: 
./launch_simulations.sh --areaMin 1000 --areaMax 1000 --malMin 0 --malMax 40 --nQuad 9 --nSim 40 --gridSize 10 --nMsg 100 --simTime 60 --minInterval 5 --maxInterval 10 --sendMode "Interval" --signAlgorithm "ECDSA"

./launch_simulations.sh --areaMin 1000 --areaMax 1000 --malMin 0 --malMax 40 --nQuad 9 --nSim 40 --gridSize 10 --nMsg 100 --simTime 60 --minInterval 5 --maxInterval 10 --sendMode "Interval" --signAlgorithm "EDDSA"

./launch_simulations.sh --areaMin 1000 --areaMax 1000 --malMin 0 --malMax 40 --nQuad 9 --nSim 40 --gridSize 10 --nMsg 100 --simTime 60 --minInterval 5 --maxInterval 10 --sendMode "FixedMaximum" --signAlgorithm "ECDSA"

./launch_simulations.sh --areaMin 1000 --areaMax 1000 --malMin 0 --malMax 40 --nQuad 9 --nSim 40 --gridSize 10 --nMsg 100 --simTime 60 --minInterval 5 --maxInterval 10 --sendMode "FixedMaximum" --signAlgorithm "EDDSA"
comment



