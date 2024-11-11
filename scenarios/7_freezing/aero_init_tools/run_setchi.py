#!/usr/bin/env python

from set_chi import classChi, classChi_bins
import numpy as np
import sys

def read_mass_ratio_dict(compFile):
    with open(compFile) as f:
        lines = []
        for line in f:
            lines.append(line)
    mass_ratio_dict = {}
    for line in lines:
        if line.strip()[0] == "#":
            continue
        infoList = line.strip().split() 
        if len(infoList) == 2:
            mass_ratio_dict[infoList[0]] = float(infoList[1])
    return mass_ratio_dict

#compFile = "../aero_init_comp.dat"
#mass_ratio_dict = read_mass_ratio_dict(compFile)
#print(mass_ratio_dict)
#exit()

chi = float(sys.argv[1])
init_fileName = sys.argv[2]
isbin = sys.argv[3]
assert isbin in ["T", "F"]
if isbin == "T":
    isbin = True
else:
    isbin = False

#compFile = sys.argv[3]
print("===============")
print(chi, init_fileName, ("Bin mode" if isbin else "Global mode"))#, compFile)
print("===============")
#mass_ratio_dict = read_mass_ratio_dict(compFile)
if isbin:
    objChi = classChi_bins(init_fileName = init_fileName, Nbins = 100)
else:
    objChi = classChi(init_fileName = init_fileName)
if isbin:
    objChi.set_chi(chi = chi)
    objChi.show_chi_bins()
else:
    while True:
        objChi.set_chi(chi = chi)#, mass_ratio_dict = mass_ratio_dict)
        if np.abs(objChi.chi - chi) < 0.01:
            print("Set chi = ", objChi.chi, ".")
            break
        else:
            print("Chi expected = ", chi, "; chi received = ", objChi.chi, ", redo ...")
objChi.output(outFile = "restart.nc")
    
