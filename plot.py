#!/usr/bin/python3

import json
from cProfile import label
import matplotlib.pyplot as plt
import numpy as np
from scipy.interpolate import interp1d
import matplotlib.lines as mlines
from itertools import chain
import argparse


def compare_arrivi(area_size, mal_per_min, mal_per_max, n_sim, sim_time, sendMode):
    plt.figure(figsize=(10, 6))
    x_common = np.linspace(0, sim_time, 100)
    line_styles = ['-', '--', '-.', ':']
    colors = plt.cm.Dark2(np.linspace(0, 1, n_sim))
    max_msg = 0
    for mal in range(mal_per_min, mal_per_max, 10):
        cumulative_y = []
        for i in range(1, n_sim + 1):
            try:
                with open(
                        f"./results/json/Scenario{area_size}_{mal}_{i}_{sendMode}.json",
                        'r') as f:
                    data = json.load(f)
                total_msg = data['totalMsg']
                if (total_msg > max_msg):
                    max_msg = total_msg
                messages = data['messages']
                tempi_arrivo = sorted( float(m['time']) for m in messages )

                y_cumulative = np.arange(1, len(tempi_arrivo) + 1)

                if len(tempi_arrivo) > 1:
                    f_interp = interp1d(tempi_arrivo, y_cumulative, bounds_error=False,
                                        fill_value=(0, y_cumulative[-1]))
                    y_interp = f_interp(x_common)
                    cumulative_y.append(y_interp)
                    
                else:
                    cumulative_y.append(np.zeros_like(x_common))
            except:
                print(f"compare_arrivi - {area_size}_{mal}_{i}_{sendMode} - Error")
                continue



        cumulative_y_mean = np.mean(cumulative_y, axis=0)
        line_style = line_styles[int(mal/10) % len(line_styles)]
        color = colors[int(mal/10) % len(colors)]
        plt.plot(x_common, cumulative_y_mean, label=f'{mal}%', linewidth=3, color=color, linestyle=line_style)

    plt.xlabel("Tempo (s)")
    plt.ylabel("Arrivi cumulativi")
    plt.title(f"Arrivi cumulativi nel tempo - Scenario {area_size}_{sendMode}")
    plt.legend()
    plt.grid(True)

    plt.ylim(0, max_msg)
    plt.yticks(np.arange(0, max_msg+1, 30))

    plt.savefig(f"./plots/comparison_arrive_plot_{area_size}_{sendMode}.png", format="png", dpi=300)
    print(f"Saved figure: ./plots/comparison_arrive_plot_{area_size}_{sendMode}.png")


def plot_arrivi_cumulativi_scenario(area_size, mal_per, n_sim, sim_time, sendMode):
    x_common = np.linspace(0, sim_time, 100)
    cumulative_y = []

    plt.figure(figsize=(10, 6))

    line_styles = ['-', '--', '-.', ':']
    colors = plt.cm.Dark2(np.linspace(0, 1, n_sim))

    max_msg = 0

    for i in range(1, n_sim + 1):
        try:
            with open(
                    f"./results/json/Scenario{area_size}_{mal_per}_{i}_{sendMode}.json",
                    'r') as f:
                data = json.load(f)
                
            total_msg = data['totalMsg']
            if (total_msg > max_msg):
                max_msg = total_msg
                
            messages = data['messages']
            tempi_arrivo = sorted(float(m['time']) for m in messages)
        except:
            print(f"plot_arrivi_cumulativi - {area_size}_{mal_per}_{i}_{sendMode} - Error")
            continue

        y_cumulative = np.arange(1, len(tempi_arrivo) + 1)

        if len(tempi_arrivo) > 1:
            f_interp = interp1d(tempi_arrivo, y_cumulative, bounds_error=False, fill_value=(0, y_cumulative[-1]))
            y_interp = f_interp(x_common)
            cumulative_y.append(y_interp)

            line_style = line_styles[i % len(line_styles)]
            color = colors[i % len(colors)]
            plt.plot(x_common, y_interp, color=color, linestyle=line_style, alpha=0.7,
                     linewidth=1.5)
        else:
            cumulative_y.append(np.zeros_like(x_common))

    cumulative_y_mean = np.mean(cumulative_y, axis=0)

    plt.plot(x_common, cumulative_y_mean, label=f'Media', color='firebrick',
             linewidth=3)
    plt.xlabel("Tempo (s)")
    plt.ylabel("Arrivi cumulativi")
    plt.title(f"Arrivi cumulativi nel tempo - Scenario {area_size}_{mal_per}_{sendMode}")
    plt.legend()
    plt.grid(True)

    plt.ylim(0, max_msg)
    plt.yticks(np.arange(0, max_msg+1, 30))

    plt.savefig(f"./plots/arrive_plot_{area_size}_{mal_per}_{sendMode}.png", format="png", dpi=300)
    print(f"Saved figure: ./plots/arrive_plot_{area_size}_{mal_per}_{sendMode}.png")


def get_times(area_size, mal_per_min, mal_per_max, n_sim, sendMode):
    plt.figure(figsize=(10, 6))
    x =[]
    y_gt = []
    y_pt = []
    for mal in range(mal_per_min, mal_per_max, 10):
        x.append(mal)
        pt = []
        gt = []
        for i in range(1, n_sim + 1):
            try:
                with open(
                        f"./results/json/Scenario{area_size}_{mal}_{i}_{sendMode}.json",
                        'r') as f:
                    data = json.load(f)
                pt.append(float(data['PerimeterTime']))
                gt.append(float(data['GreedyTime']))

            except:
                print(f"get_times - {area_size}_{mal}_{i}_{sendMode} - Error")
                continue

        pt_mean = np.mean(pt)
        gt_mean = np.mean(gt)
        y_gt.append(gt_mean)
        y_pt.append(pt_mean)

    plt.plot(x, y_gt, label=f'Tempo medio in Greedy')
    plt.plot(x, y_pt, label=f'Tempo medio in Perimeter')
    plt.scatter(x, y_gt)
    plt.scatter(x, y_pt)

    plt.xlabel('Percentuale nodi malevoli (%)')
    plt.ylabel('Tempo medio')
    plt.title('Tempi medi in Perimeter e Greedy')
    plt.legend()
    plt.savefig(f"./plots/mean_time_plot_{area_size}_{sendMode}.png", format="png", dpi=300)
    print(f"Saved figure: ./plots/mean_time_plot_{area_size}_{sendMode}.png")

def mean_hops(area_size, mal_per_min, mal_per_max, n_sim, sendMode):
    plt.figure(figsize=(10, 6))
    x = []
    y = []
    for mal in range(mal_per_min, mal_per_max, 10):
        x.append(mal)
        mean = []
        for i in range(1, n_sim+1):
            try:
                with open(f"./results/json/Scenario{area_size}_{mal}_{i}_{sendMode}.json",
                        'r') as f:
                    data = json.load(f)
                    messages = data['messages']
                    hops = []
                    for m in messages:
                        hops.append(int(m['hops']))
                    mean.append(np.mean(hops))
            except:
                print(f"mean_hops - {area_size}_{mal}_{i}_{sendMode} - Error")
                continue

        y.append(np.mean(mean))

    plt.plot(x, y)
    plt.scatter(x, y)

    plt.xlabel('Percentuale nodi malevoli (%)')
    plt.ylabel('hops medi')
    plt.title('hops medi (arrived)')
    plt.savefig(f"./plots/mean_hops_plot_{area_size}_{sendMode}.png", format="png", dpi=300)
    print(f"Saved figure: ./plots/mean_hops_plot_{area_size}_{sendMode}.png")
    
def mean_hops_not_arrived(area_size, mal_per_min, mal_per_max, n_sim, sendMode):
    plt.figure(figsize=(10, 6))
    x = []
    y = []
    for mal in range(mal_per_min, mal_per_max, 10):
        x.append(mal)
        mean = []
        for i in range(1, n_sim+1):
            try:
                with open(f"./results/json/Scenario{area_size}_{mal}_{i}_{sendMode}.json",
                        'r') as f:
                    data = json.load(f)
                    messages = data['notArrivedMessages']
                    hops = []
                    for m in messages:
                        hops.append(int(m['hops']))
                    mean.append(np.mean(hops))
            except:
                print(f"mean_hops_not_arrived - {area_size}_{mal}_{i}_{sendMode} - Error")
                continue

        y.append(np.mean(mean))

    plt.plot(x, y)
    plt.scatter(x, y)

    plt.xlabel('Percentuale nodi malevoli (%)')
    plt.ylabel('hops medi')
    plt.title('hops medi (not arrived)')
    plt.savefig(f"./plots/mean_hops_plot_not_arrived_{area_size}_{sendMode}.png", format="png", dpi=300)
    print(f"Saved figure: ./plots/mean_hops_plot_not_arrived_{area_size}_{sendMode}.png")

def main():
    parser = argparse.ArgumentParser(description="Script per plottare")
    parser.add_argument('--areaMin', type=int, required=True, help='Area minima')
    parser.add_argument('--areaMax', type=int, required=True, help='Area massima')
    parser.add_argument('--malMin', type=int, required=True, help='Minima percentuale di nodi malevoli')
    parser.add_argument('--malMax', type=int, required=True, help='Massima percentuale di nodi malevoli')
    parser.add_argument('--nSim', type=int, required=True, help='Numero di simulazioni')
    parser.add_argument('--simTime', type=int, required=True, help='Durata simulazione')
    parser.add_argument('--sendMode', type=str, required=True, help='Send mode (FixedMaximum o Interval)')

    args = parser.parse_args()

    area_min = args.areaMin
    area_max = args.areaMax
    min_mal = args.malMin
    max_mal = args.malMax
    n_sim = args.nSim
    sim_time = args.simTime
    sendMode = args.sendMode

    for i in range(area_min, area_max+1, 1000):
        for j in range(min_mal, max_mal+1, 10):
            plot_arrivi_cumulativi_scenario(i, j, n_sim, sim_time, sendMode)

    for i in range(area_min, area_max+1, 1000):
        get_times(i, min_mal, max_mal+1, n_sim, sendMode)

    for i in range(area_min, area_max+1, 1000):
        compare_arrivi(i, min_mal, max_mal+1, n_sim, sim_time, sendMode)

    for i in range(area_min, area_max+1, 1000):
        mean_hops(i, min_mal, max_mal+1, n_sim, sendMode)
    
    for i in range(area_min, area_max+1, 1000):
        mean_hops_not_arrived(i, min_mal, max_mal+1, n_sim, sendMode)


if __name__ == "__main__":
    main()



