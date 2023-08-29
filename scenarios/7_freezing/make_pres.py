#!/usr/bin/env python

import numpy as np
from pdb import set_trace
from scipy.interpolate import interp1d
import pandas as pd

min_SI = 60
hour_SI = 3600
day_SI = 24 * 3600


#total_hour = 24
#timeList = np.arange(0, 3600 * (total_hour + 1), 3600)
total_time = 3600
dtime = 60
timeList = np.arange(0, total_time + dtime, dtime)
Nt = len(timeList)

def change_linear(start_T, end_T, Nt):
    return np.linspace(start_T, end_T, Nt)

def change_linear_stages(tList, tempList):
    f = interp1d(tList, tempList)
    temperature = f(timeList)
    return temperature

def change_curve_file(fileName):
    df = pd.read_csv(fileName, header = None)
    x = df[0].values
    y = df[1].values
    f = interp1d(x, y)
    data = f(timeList)
    return data

def make_file(timeList, TList, fileName = "pressure.dat"):
    with open(fileName, "w") as fw:
        fw.write("#   time    (s)\n")
        fw.write("# pressure (Pa)\n")
        fw.write("time ")
        for time in timeList:
            fw.write("%d    " % (time))
        fw.write("\n")
        fw.write("pressure ")
        for T in TList:
            fw.write("%.3f    " % (T))
        fw.write("\n")


if __name__ == "__main__":
    print(Nt)
    TLists = [
        #change_linear(288, 240, 25),
        #change_linear(240, 288, 24),
        #change_linear(240, 240, 5),
        #change_linear(240, 278, 19),
        #change_linear(240, 240, 24),

        #change_linear(265, 225, 61),
        #change_linear(250, 250, total_hour + 1 - 13)
        change_linear_stages(
            #tList = [0,          5 * min_SI,  15 * min_SI, 20 * min_SI,  30 * min_SI, 35 * min_SI, 45 * min_SI, 50 * min_SI, total_time],
            #tempList = [280 , 237.5 , 280 ,  237.5 , 280 ,  237.5 , 280 , 237.5 , 280 ],
            #tList = [0, 14 * min_SI, 15 * min_SI, total_time],
            #tempList  = [270 , 270 ,  225 ,  225],
            tList = [0, total_time],
            tempList  = [100000, 100000],

        )
        #change_curve_file("/data/keeling/a/wenhant2/Scripts/analysis/3_freezing_BenchPlots/Mohler2008_data/Mohler2008_fig3_P_black.csv") * 100
    ]
    TList = np.concatenate(TLists)
    assert len(timeList) == len(TList)
    make_file(timeList, TList)





