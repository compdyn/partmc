#!/usr/bin/env python

from set_chi import classChi
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
compFile = sys.argv[3]
print("===============")
print(chi, init_fileName, compFile)
print("===============")
mass_ratio_dict = read_mass_ratio_dict(compFile)
objChi = classChi(init_fileName = init_fileName)
while True:
    objChi.set_chi(chi = chi)#, mass_ratio_dict = mass_ratio_dict)
    if np.abs(objChi.chi - chi) < 0.01:
        print("Set chi = ", objChi.chi, ".")
        break
    else:
        print("Chi expected = ", chi, "; chi received = ", objChi.chi, ", redo ...")
objChi.output(outFile = "restart.nc")
    
