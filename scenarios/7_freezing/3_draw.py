#!/usr/bin/env python

import numpy as np
import glob
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker

xlabel_fontsize = 15
ylabel_fontsize = 15

def read_expFile():
    fileNameList = glob.glob("out/freezing_part_*_aero_time.txt")
    timeList = []
    tempList = []
    ice_ratio_mean = []
    ice_ratio_max = []
    ice_ratio_min = []
    ice_ratio = []

    for i_file, fileName in enumerate(fileNameList):
        ice_ratio_i_file = []
        with open(fileName, "r") as fr:
            for iline in fr:
                time, temp, ff = iline.strip().split()
                if i_file == 0:
                    timeList.append(float(time))
                    tempList.append(float(temp))
                ice_ratio_i_file.append(float(ff))
        ice_ratio.append(ice_ratio_i_file)

    timeList = np.array(timeList)
    tempList = np.array(tempList)
    ice_ratio = np.array(ice_ratio)

    ice_ratio_mean = ice_ratio.mean(axis = 0)
    ice_ratio_max = ice_ratio.max(axis = 0)
    ice_ratio_min = ice_ratio.min(axis = 0)
    return {
        "timeList": timeList,
        "tempList": tempList,
        "ice_ratio_mean": ice_ratio_mean,
        "ice_ratio_max": ice_ratio_max,
        "ice_ratio_min": ice_ratio_min
    }

def draw(fig, ax):
    OutDir = "output"
    output = read_expFile()
    timeList = output["timeList"]
    temperature  = output["tempList"]
    axt = ax.twinx()
    axt.plot(timeList / 60, temperature - 273.15, color = "red")
    axt.set_yticks([-40, -30, -20, -10])
    axt.set_ylim([-102.5, -5])
    axt.grid(linestyle = "--")
    ice_ratio_mean = output["ice_ratio_mean"]
    ice_ratio_max =  output["ice_ratio_max"]
    ice_ratio_min =  output["ice_ratio_min"]
    ax.plot(
        timeList / 60,
        ice_ratio_mean * 100,
        color = "blue",
        linewidth = 0.7
    )
    ax.fill_between(
        timeList / 60,
        ice_ratio_min * 100,
        ice_ratio_max * 100,
        color = "blue",
        alpha = 0.3,
        edgecolor = None
    )
    ax.set_yticks([0, 20, 40, 60, 80])
    ax.set_ylim([-15, 140])
    ax.grid()
    ax.xaxis.set_major_locator(ticker.MultipleLocator(1))
    for label in axt.get_yticklabels():
        label.set_fontsize(12)
    axt.set_ylabel("Temperature (˚C)", loc = "top", fontsize = ylabel_fontsize)
    for label in ax.get_yticklabels():
        label.set_fontsize(12)
    for label in ax.get_xticklabels():
        label.set_fontsize(12)
    ax.set_ylabel("Frozen fraction (%)", loc = "bottom", fontsize = ylabel_fontsize)
    return fig, ax


if __name__ == "__main__":    

    fig = plt.figure(figsize = (10, 5))
    ax = fig.add_subplot(1, 1, 1)
    fig, ax = draw(fig, ax)
    ax.set_xlabel("Time (min)", fontsize = xlabel_fontsize)

    plt.savefig("out/TSs.pdf")
    print("The figure has been saved to out/TSs.pdf")
