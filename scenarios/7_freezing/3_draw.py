#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import matplotlib.ticker as ticker

cmap = plt.cm.Blues
colors = cmap(np.linspace(0.4, 1, 10))  
chi_list =  np.array([0, 0.01, 0.05, 0.1, 0.15, 0.2, 0.4, 0.6, 0.8, 1])
chi_colors = cmap(0.3 + 0.7 * chi_list)
chi_colors = chi_colors[::-1]

xlabel_fontsize = 15
ylabel_fontsize = 15

def read_expFile(caseName):
    fileName = "out/freezing_part_data_" + caseName + ".txt"
    timeList = []
    tempList = []
    ice_ratio_mean = []
    ice_ratio_max = []
    ice_ratio_min = []
    with open(fileName, "r") as fr:
        for iline in fr:
            time, temp, ff_mean, ff_max, ff_min = iline.strip().split()
            timeList.append(float(time))
            tempList.append(float(temp))
            ice_ratio_mean.append(float(ff_mean))
            ice_ratio_max.append(float(ff_max))
            ice_ratio_min.append(float(ff_min))
    timeList = np.array(timeList)
    tempList = np.array(tempList)
    ice_ratio_mean = np.array(ice_ratio_mean)
    ice_ratio_max = np.array(ice_ratio_max)
    ice_ratio_min = np.array(ice_ratio_min)
    return timeList, tempList, ice_ratio_mean, ice_ratio_max, ice_ratio_min

        


def draw(fig, ax, casesName, ax_label):
    OutDir = "output"
    casesLabel = ["100% $Illite$", "100% $Fe_2O_3$", "Internal Mixture", "External Mixture"]
    casesColor = ["green", "grey", chi_colors[0], chi_colors[-1]]

    sampleName = casesName[0]
    casesDict = {}

    for caseName in casesName:
        timeList, tempList, ice_ratio_mean, ice_ratio_max, ice_ratio_min = read_expFile(caseName)
        CaseDict = {"timeList": timeList, "tempList": tempList, "ice_ratio_mean": ice_ratio_mean,
                    "ice_ratio_max": ice_ratio_max, "ice_ratio_min": ice_ratio_min}
        casesDict[caseName] = CaseDict

    timeList = casesDict[sampleName]["timeList"]
    timeList = timeList / 60
    temperature  = casesDict[sampleName]["tempList"]

    axt = ax.twinx()
    axt.plot(timeList, temperature - 273.15, color = "red")
    axt.set_yticks([-40, -30, -20, -10])
    axt.set_ylim([-102.5, -5])
    axt.grid(linestyle = "--")


    for ind, caseName in enumerate(casesName):
        ice_ratio_mean = casesDict[caseName]["ice_ratio_mean"]
        ice_ratio_max =  casesDict[caseName]["ice_ratio_max"]
        ice_ratio_min =  casesDict[caseName]["ice_ratio_min"]
        ax.plot(timeList, ice_ratio_mean * 100, label = casesLabel[ind], color = casesColor[ind], linewidth = 0.7)
        ax.fill_between(timeList, ice_ratio_min * 100, ice_ratio_max * 100, color = casesColor[ind], alpha = 0.3, edgecolor = None)
        #print(caseName, ice_ratio_mean[-1] * 100)

    ax.set_yticks([0, 20, 40, 60, 80])
    ax.set_ylim([-15, 140])
    ax.grid()

    ax.xaxis.set_major_locator(ticker.MultipleLocator(1))
    for label in axt.get_yticklabels():
        label.set_fontsize(12)
        label.set_fontname("serif")
    axt.set_ylabel("Temperature (˚C)", loc = "top", fontsize = ylabel_fontsize)
    for label in ax.get_yticklabels():
        label.set_fontsize(12)
        label.set_fontname("serif")
    for label in ax.get_xticklabels():
        label.set_fontsize(12)
        label.set_fontname("serif")


    ax.set_ylabel("Frozen fraction (%)", loc = "bottom", fontsize = ylabel_fontsize)
    
    ax.text(-0.01, 0.96, ax_label, horizontalalignment = "right", verticalalignment = "top", transform = ax.transAxes, fontsize = 15)
    return fig, ax


if __name__ == "__main__":    

    fig, axes = plt.subplots(2, 1, figsize = (10, 7), sharex = True)
    plt.subplots_adjust(hspace = 0.1)
    ax1, ax2 = axes[0], axes[1]
    casesName = ["exp1", "exp2", "exp3", "exp4"]
    fig, ax1 = draw(fig, ax1, casesName, ax_label = "(a)")

    casesName = ["exp5", "exp6", "exp7", "exp8"]
    fig, ax2 = draw(fig, ax2, casesName, ax_label = "(b)")
    ax2.set_xlabel("Time (min)", fontsize = xlabel_fontsize)

    #plt.show()
    plt.savefig("out/TSs.png", dpi = 500)
    print("The figure has been saved to out/TSs.png")
