#!/usr/bin/env python

import numpy as np
import netCDF4 as nc
import matplotlib.pyplot as plt
from pdb import set_trace
from copy import deepcopy
import os
import time

#np.random.seed(100)
seed = int(int(time.time() * 1e6)%1e8)
print("Set random seed to ", seed)
np.random.seed(seed)

class classChi_bins(object):
    def __init__(self, init_fileName, Nbins = 50):

        self.Nbins = Nbins
        self.init_fileName = init_fileName
        self.init_dir = "/".join(init_fileName.split("/")[:-1])
        ncf = nc.Dataset(init_fileName)
        Npart = ncf.dimensions["aero_particle"].size
       
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
        #print(np.sum(self.species_volume, axis=1))
        self.Vp_core = np.sum(self.species_volume, axis = 0)
        #self.Vp_core = np.sum(self.species_volume[self.exc_H2O_ind, :], axis = 0)
        #assert len(self.species_name) == len(self.species_dens)
        #self.aero_dens = {}
        #for ind, speciesName in enumerate(self.species_name):
        #    self.aero_dens[speciesName] = self.species_dens[ind, 0]
        ncf.close()
        self.log10_Vp_core = np.log10(self.Vp_core)
        _, bin_edges = np.histogram(self.log10_Vp_core, bins = self.Nbins)
        self.bin_edges = bin_edges
        self.index_dict = {}
        for ind, (edge1, edge2) in enumerate(zip(bin_edges[:-1], bin_edges[1:])):
            if ind == 0:
                choosed_ind = np.where( (self.log10_Vp_core >= edge1) & (self.log10_Vp_core <= edge2) )[0]
            else:
                choosed_ind = np.where( (self.log10_Vp_core > edge1) & (self.log10_Vp_core <= edge2) )[0]
            #print(edge1, edge2)
            #print(choosed_ind)
            self.index_dict[ind] = choosed_ind.copy()
        count = 0
        for ind in self.index_dict:
            count += len(self.index_dict[ind])
        print("Npart = ", self.Npart, "; counted = ", count)

    def set_chi(self, chi):
        orig_species_total_volume = self.species_volume.sum(axis = 1)
        for ibin in self.index_dict:
            part_index = self.index_dict[ibin]
            if len(part_index) < 2:
                continue
            species_volume_ibin = self.set_chi_ibin(chi, ibin)
            #set_trace()
            self.species_volume[:, part_index] = species_volume_ibin
        self.update()
        new_species_total_volume = self.species_volume.sum(axis = 1)
        print("=" * 40)
        print("Species total volume difference:")
        print("orig = ", orig_species_total_volume)
        print("new = ", new_species_total_volume)
        print("=" * 40)
        
    def set_chi_ibin(self, chi, ibin):
        
        part_index = self.index_dict[ibin]
        

        objChi = classChi(init_fileName = self.init_fileName, specific_part_ind = part_index)
        count = 0
        max_count = 200
        best_objChi = None
        min_diff = np.inf
        while True:
            count += 1
            #objChi.set_chi(chi = chi, mass_ratio_dict = mass_ratio_dict)
            objChi.set_chi(chi = chi, Npars = 1)
            diff = np.abs(objChi.chi - chi)
            if diff < min_diff:
                min_diff = diff
                best_objChi = deepcopy(objChi)
            if np.abs(objChi.chi - chi) < 0.01:
                break
            if count >= max_count:
                print("Warning, exceed the maximum times of iteration, force to stop")
                print("Current best chi = ", best_objChi.chi, "; diff = ", np.abs(best_objChi.chi - chi))
                objChi = best_objChi
                break
            else:
                print("Chi expected = ", chi, "; chi received = ", objChi.chi, ", rerun ...")
        assert np.abs(objChi.chi - chi) == min_diff
        print("Final chi = ", objChi.chi)
        return objChi.species_volume

    def update(self):
        self.species_volume_ratio = self.species_volume / self.Vp_core
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

    def show_chi_bins(self):
        for ibin in self.index_dict:
            part_index = self.index_dict[ibin]
            if len(part_index) == 0:
                print("bin #", ibin, ": no particle")
            else:
                chi_ibin = self.compute_chi_partial(part_index)
                print("bin #", ibin, ": len = ", len(part_index), "; chi = ", chi_ibin)

    def compute_chi_partial(self, index_list):

        species_mass = self.species_mass[:, index_list]
        aero_mass = self.aero_mass[index_list]
        species_total_mass = np.sum(species_mass, axis = 1)
        total_mass = np.sum(aero_mass)

        P_ai = species_mass / aero_mass
        P_i = aero_mass / total_mass
        P_a = species_total_mass / total_mass

        Hi = -np.sum(np.where(P_ai == 0.0, 0.0, P_ai * np.log(P_ai)), axis = 0)
        Ha = np.sum(P_i * Hi)
        Hy = -np.sum(np.where(P_a == 0.0, 0.0, P_a * np.log(P_a)))
        Di = np.exp(Hi)
        Da = np.exp(Ha)
        Dy = np.exp(Hy)
        Db = Dy / Da
        chi = (Da - 1) / (Dy - 1)
        return chi

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

    def __init__(self, init_fileName, specific_part_ind = None):
        if specific_part_ind is None:
            self.reginal_mode = False
        else:
            self.reginal_mode = True
        self.init_fileName = init_fileName
        self.init_dir = "/".join(init_fileName.split("/")[:-1])
        ncf = nc.Dataset(init_fileName)
        Npart = ncf.dimensions["aero_particle"].size
        if not(self.reginal_mode):
            specific_part_ind = np.arange(Npart)
        aero_particle_mass = ncf.variables["aero_particle_mass"][:, specific_part_ind].filled(np.nan)
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
        #print(np.sum(self.species_volume, axis=1))
        self.Vp_core = np.sum(self.species_volume, axis = 0)
        #self.Vp_core = np.sum(self.species_volume[self.exc_H2O_ind, :], axis = 0)
        #assert len(self.species_name) == len(self.species_dens)
        #self.aero_dens = {}
        #for ind, speciesName in enumerate(self.species_name):
        #    self.aero_dens[speciesName] = self.species_dens[ind, 0]
        ncf.close()




    def init_random(self):

        self.species_volume_ratio  = np.random.rand(self.Ns, self.Npart)
        self.species_volume_ratio_one_parm = np.sum(self.species_volume_ratio, axis = 0)
        self.species_volume_ratio = self.species_volume_ratio / self.species_volume_ratio_one_parm
        self.species_volume = self.species_volume_ratio * self.Vp_core
        self.update()

    def init_external(self, species_total_volume_ratio):
        self.species_volume_ratio = np.full((self.Ns, self.Npart), 0.0)
        species_total_volume_ratio = np.array(species_total_volume_ratio)
        species_total_volume_ratio = species_total_volume_ratio / np.sum(species_total_volume_ratio)
        aero_ind = np.arange(self.Npart)
        np.random.shuffle(aero_ind)
        species_volume_cumsum = np.cumsum(self.Vp_core[aero_ind])
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
        self.species_volume = self.species_volume_ratio * self.Vp_core
        orig_species_total_volume = species_total_volume_ratio * np.sum(self.Vp_core)
        new_species_total_volume = self.species_volume.sum(axis = 1)

        #print("species_total_volume difference after init_external:")
        #for orig, new in zip(orig_species_total_volume, new_species_total_volume):
        #    print(orig, new)
        #print("="*30)
        #print("Np = ", len(self.Vp_core))
        #print(orig_species_total_volume.sum(), new_species_total_volume.sum())

        #species_total_volume = (self.species_volume_ratio * self.Vp).sum(axis=1)
        #species_total_volume_ratio = species_total_volume / np.sum(species_total_volume)
        self.update()

    def update(self):
        self.species_volume_ratio = self.species_volume / self.Vp_core
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

    """
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
    """

    def mixing(self, Npars, chi_obj = None):
        assert Npars * 2 <= self.Npart
        aero_index = np.arange(self.Npart)
        np.random.shuffle(aero_index)
        parA_ind = aero_index[:Npars]
        parB_ind = aero_index[-Npars:]
        parA_vol = self.Vp_core[parA_ind]
        parB_vol = self.Vp_core[parB_ind]
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

    def set_chi(self, chi, Npars = 10):
        #volume_ratio = []
        #for specName in self.aero_dens:
        #    if specName in mass_ratio_dict:
        #        volume_ratio.append(mass_ratio_dict[specName] / self.aero_dens[specName])
        #    else:
        #        volume_ratio.append(0)
        #volume_ratio = np.array(volume_ratio)
        #volume_ratio /= np.sum(volume_ratio)
        volume_ratio = np.sum(self.species_volume,axis=1)/np.sum(self.species_volume)
        #print("volume_ratio = ", volume_ratio)
        self.init_external(volume_ratio)
        chi_list = []
        while True:
            #self.mixing(Npars = int(self.Npart / 200) if int(self.Npart / 200) > 1 else 1)#, chi_obj = chi_obj)
            self.mixing(Npars = Npars)#, chi_obj = chi_obj)
    #a.compute_chi()
            print("chi = ", self.chi)
            #print(np.sum(self.species_volume, axis=1))
            #print(np.sum(self.species_mass, axis=1))
            #print(np.sum(self.species_volume, axis=0)[500])
            chi_list.append(self.chi)
            if chi <= self.chi:
                break
        #print(self.species_mass)

    def output(self, outFile = None):
        assert self.reginal_mode == False
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





#if __name__ == "__main__":
if __name__ == "__1__":
    chi = 0.2
    mass_ratio_dict = {
        #"OIN": 0.3,
        "ILT": 0.5,
        #"KLN": 0.2,
        "NVF": 0.5,
    }
    aa = classChi(init_fileName = "/data/keeling/a/wenhant2/modeling/partmc/scenarios/7_freezing/output/chiexp_indsize_0.8/freezing_part_0001_00000001.nc", specific_part_ind = np.arange(900, 902))
    count = 0
    max_count = 200
    best_aa = None
    min_diff = np.inf
    while True:
        count += 1
        #aa.set_chi(chi = chi, mass_ratio_dict = mass_ratio_dict)
        aa.set_chi(chi = chi)
        diff = np.abs(aa.chi - chi)
        if diff < min_diff:
            min_diff = diff
            best_aa = deepcopy(aa)
        if np.abs(aa.chi - chi) < 0.01:
            break
        if count >= max_count:
            print("Warning, exceed the maximum times of iteration, force to stop")
            print("Current best chi = ", best_aa.chi, "; diff = ", np.abs(best_aa.chi - chi))
            aa = best_aa
            break
        else:
            print("Chi expected = ", chi, "; chi received = ", aa.chi, ", rerun ...")
    assert np.abs(aa.chi - chi) == min_diff
    print("Final chi = ", aa.chi)
    #aa.output(outFile = "freezing_part_0001_restart.nc")
if __name__ == "__main__":
    aa = classChi_bins(init_fileName = "./freezing_part_0001_00000001.nc", Nbins = 100)
    aa.set_chi(chi = 0.6)
    #print(aa.chi)
    aa.show_chi_bins()
    print("Overall chi = ", aa.chi)
    #aa.output(outFile = "freezing_part_0001_restart.nc")
    
