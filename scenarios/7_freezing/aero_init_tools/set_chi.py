#!/usr/bin/env python

import numpy as np
import netCDF4 as nc
import matplotlib.pyplot as plt
from pdb import set_trace
import os

class classChi(object):

    """
    def __init__(self, mean = 1e-6, logstd = 0.326, Npart = 1000, species = ["OIN", "ILT", "KLN"]):

        self.Npart = Npart
        self.species = species
        self.Ns = len(species)

        self.aero_dens = {
            "H2O": 1000,
            "OIN": 2600,
            "ILT": 2750,
            "KLN": 2650,
        }
        self.species_dens = []
        for spec in self.species:
            self.species_dens.append(self.aero_dens[spec])
        self.species_dens = np.array(self.species_dens)
        self.species_dens = self.species_dens.reshape(self.Ns, 1)
        self.logDp = np.random.normal(np.log10(mean), logstd, size = (self.Npart))
        self.Dp = 10 **self.logDp
        self.Vp = 1.0/6 * np.pi * self.Dp **3
    """

    def __init__(self, init_fileName):
        self.init_fileName = init_fileName
        self.init_dir = "/".join(init_fileName.split("/")[:-1])
        ncf = nc.Dataset(init_fileName)
        aero_particle_mass = ncf.variables["aero_particle_mass"][:].filled(np.nan)
        self.species_dens = ncf.variables["aero_density"][:].filled(np.nan)
        self.species_dens = self.species_dens[:, np.newaxis]
        self.species_name = ncf.variables["aero_species"].names.strip().split(",")
        self.Ns_orig, _ = aero_particle_mass.shape
        H2O_ind = np.nan
        exc_H2O_ind = []
        for ind, specName in enumerate(self.species_name):
            if specName.lower() != "h2o":
                exc_H2O_ind.append(ind)
            else:
                assert np.isnan(H2O_ind), "Double \"H2O\" found."
                H2O_ind = ind
        self.H2O_ind = H2O_ind
        self.exc_H2O_ind = exc_H2O_ind
        if ~np.isnan(self.H2O_ind):
            self.H2O_mass = aero_particle_mass[self.H2O_ind, :]

        aero_particle_mass = aero_particle_mass[exc_H2O_ind, :]
        self.species_dens = self.species_dens[exc_H2O_ind, :]
        #set_trace()
        #self.species_name = self.species_name[exc_H2O_ind]

        self.Ns, self.Npart = aero_particle_mass.shape
        if np.isnan(self.H2O_ind):
            print("No H2O found.")
            assert self.Ns_orig == self.Ns
        self.species_volume = aero_particle_mass / self.species_dens 
        print(np.sum(self.species_volume, axis=1))
        self.Vp = np.sum(self.species_volume, axis = 0)
        #assert len(self.species_name) == len(self.species_dens)
        #self.aero_dens = {}
        #for ind, speciesName in enumerate(self.species_name):
        #    self.aero_dens[speciesName] = self.species_dens[ind, 0]
        ncf.close()




    def init_random(self):

        self.species_volume_ratio  = np.random.rand(self.Ns, self.Npart)
        self.species_volume_ratio_one_parm = np.sum(self.species_volume_ratio, axis = 0)
        self.species_volume_ratio = self.species_volume_ratio / self.species_volume_ratio_one_parm
        self.species_volume = self.species_volume_ratio * self.Vp
        self.update()

    def init_external(self, species_total_volume_ratio):
        self.species_volume_ratio = np.full((self.Ns, self.Npart), 0.0)
        species_total_volume_ratio = np.array(species_total_volume_ratio)
        species_total_volume_ratio = species_total_volume_ratio / np.sum(species_total_volume_ratio)
        aero_ind = np.arange(self.Npart)
        np.random.shuffle(aero_ind)
        species_volume_cumsum = np.cumsum(self.Vp[aero_ind])
        species_volume_cumratio = species_volume_cumsum / species_volume_cumsum[-1]
        threshold_ratio = np.cumsum(species_total_volume_ratio)
        threshold = [0]
        for ratio in threshold_ratio:
            if ratio == 0:
                 threshold.append(0)
            else:
                threshold.append(np.argmin(np.abs(species_volume_cumratio - ratio)) + 1)
        for i_spec in range(self.Ns):
            spec_index = aero_ind[threshold[i_spec]:threshold[i_spec+1]]
            self.species_volume_ratio[i_spec, spec_index] = 1.0
        self.species_volume = self.species_volume_ratio * self.Vp
        #species_total_volume = (self.species_volume_ratio * self.Vp).sum(axis=1)
        #species_total_volume_ratio = species_total_volume / np.sum(species_total_volume)
        self.update()

    def update(self):
        self.species_volume_ratio = self.species_volume / self.Vp
        self.species_mass = self.species_volume * self.species_dens
        self.aero_mass = np.sum(self.species_mass, axis = 0)
        self.species_total_mass = np.sum(self.species_mass, axis = 1)
        self.total_mass = np.sum(self.aero_mass)
        self.P_ai = self.species_mass / self.aero_mass
        self.P_i = self.aero_mass / self.total_mass
        self.P_a = self.species_total_mass / self.total_mass
        self.compute_chi()

    def compute_chi(self):
        self.Hi = -np.sum(np.where(self.P_ai == 0.0, 0.0, self.P_ai * np.log(self.P_ai)), axis = 0)
        self.Ha = np.sum(self.P_i * self.Hi)
        self.Hy = -np.sum(np.where(self.P_a == 0.0, 0.0, self.P_a * np.log(self.P_a)))
        self.Di = np.exp(self.Hi)
        self.Da = np.exp(self.Ha)
        self.Dy = np.exp(self.Hy)
        self.Db = self.Dy / self.Da
        self.chi = (self.Da - 1) / (self.Dy - 1)

    def mixing(self, Npars, chi_obj = None):
        assert Npars * 2 <= self.Npart
        aero_index = np.arange(self.Npart)
        np.random.shuffle(aero_index)
        parA_ind = aero_index[:Npars]
        parB_ind = aero_index[-Npars:]
        parA_vol = self.Vp[parA_ind]
        parB_vol = self.Vp[parB_ind]
        min_vol = np.min(np.array([parA_vol, parB_vol]), axis = 0)
        exchange_ratio = np.random.rand(len(min_vol)) 
        if not(chi_obj is None):
            exchange_ratio *= np.abs(self.chi - chi_obj) / np.abs(chi_obj)
        exchange_vol = exchange_ratio * min_vol
        exchange_species_vol_A = self.species_volume_ratio[:, parA_ind] * exchange_vol
        exchange_species_vol_B = self.species_volume_ratio[:, parB_ind] * exchange_vol
        self.species_volume[:, parB_ind] = self.species_volume[:, parB_ind] - exchange_species_vol_B + exchange_species_vol_A
        self.species_volume[:, parA_ind] = self.species_volume[:, parA_ind] - exchange_species_vol_A + exchange_species_vol_B
        self.update()
        #set_trace()
        #print(np.sum(self.species_volume, axis = 0))

    def set_chi(self, chi):
        #volume_ratio = []
        #for specName in self.aero_dens:
        #    if specName in mass_ratio_dict:
        #        volume_ratio.append(mass_ratio_dict[specName] / self.aero_dens[specName])
        #    else:
        #        volume_ratio.append(0)
        #volume_ratio = np.array(volume_ratio)
        #volume_ratio /= np.sum(volume_ratio)
        volume_ratio = np.sum(self.species_volume,axis=1)/np.sum(self.species_volume)
        print("volume_ratio = ", volume_ratio)
        self.init_external(volume_ratio)
        chi_list = []
        while True:
            #self.mixing(Npars = int(self.Npart / 200) if int(self.Npart / 200) > 1 else 1)#, chi_obj = chi_obj)
            self.mixing(Npars = 2)#, chi_obj = chi_obj)
    #a.compute_chi()
            print("chi = ", self.chi)
            print(np.sum(self.species_volume, axis=1))
            #print(np.sum(self.species_mass, axis=1))
            #print(np.sum(self.species_volume, axis=0)[500])
            chi_list.append(self.chi)
            if chi <= self.chi:
                break
        print(self.species_mass)

    def output(self, outFile = None):
        if outFile is None:
            ncf = nc.Dataset(self.init_fileName, "r+")
        else:
            os.system("cp -p " + self.init_fileName + " " + self.init_dir + "/" + outFile)
            ncf = nc.Dataset(self.init_dir + "/" + outFile, "r+")
            print("Create " + self.init_dir + "/" + outFile)
        aero_particle_mass = np.full((self.Ns_orig, self.Npart), 0.0)
        if ~np.isnan(self.H2O_ind):
            aero_particle_mass[self.H2O_ind, :] = self.H2O_mass
            aero_particle_mass[self.exc_H2O_ind, :] = self.species_mass
        else:
            aero_particle_mass = self.species_mass
        ncf.variables["aero_particle_mass"][:, :] = aero_particle_mass
        ncf.close()





if __name__ == "__main__":
    chi = 0.9
    mass_ratio_dict = {
        #"OIN": 0.3,
        "ILT": 0.5,
        #"KLN": 0.2,
        "NVF": 0.5,
    }
    aa = classChi(init_fileName = "/data/keeling/a/wenhant2/modeling/partmc/scenarios/7_freezing/output/chiexp_0.2/freezing_part_0001_00000001.nc")
    while True:
        aa.set_chi(chi = chi, mass_ratio_dict = mass_ratio_dict)
        if np.abs(aa.chi - chi) < 0.01:
            break
        else:
            print("Chi expected = ", chi, "; chi received = ", aa.chi, ", rerun ...")
    aa.output(outFile = "freezing_part_0001_restart.nc")
