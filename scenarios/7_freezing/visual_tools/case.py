#!/usr/bin/env python

import numpy as np
import netCDF4 as nc
from glob import glob
import os

from pdb import set_trace

class PartMC_Case(object):

    def __init__(self, caseName, OutDir = ".", prefix = "freezing_part", disable_check_time = False):
        
        #self.prefix = "freezing_part"
        self.prefix = prefix
        self.caseName = caseName
        if os.path.exists(OutDir + "/" + caseName):
            self.case_output = OutDir + "/" + caseName
        elif os.path.exists(OutDir + "/out_" + caseName):
            self.case_output = OutDir + "/out_" + caseName
        else:
            assert False, caseName + " couln'd be found!"
        self.fileList = glob(self.case_output + "/" + self.prefix + "_*.nc")
        self.fileList.sort()
        self.sample_file = self.fileList[0]
        self.ensemble_nameList = []
        for fileName in self.fileList:
            ensemble_name = fileName.split("_")[-2]
            if not(ensemble_name in self.ensemble_nameList):
                self.ensemble_nameList.append(ensemble_name)

        self.N_ensemble = len(self.ensemble_nameList)
        self.ensemble_fileList = {}
        self.N_file_per_ensemble = None
        self.default_ensemble_name = self.ensemble_nameList[0]

        for ensemble_name in self.ensemble_nameList:
            fileList = glob(self.case_output + "/" + self.prefix + "_" + ensemble_name + "_*.nc")
            fileList.sort()
            self.ensemble_fileList[ensemble_name] = fileList
            if self.N_file_per_ensemble is None:
                self.N_file_per_ensemble = len(fileList)
            else:
                assert self.N_file_per_ensemble == len(fileList)

        if disable_check_time:
            timeList = []
            for fileName in self.ensemble_fileList[self.ensemble_nameList[0]]:
                ncf = nc.Dataset(fileName)
                time = float(ncf.variables["time"][0].filled(np.nan))
                ncf.close()
                timeList.append(time)
            self.timeList = np.array(timeList)
        else:
            self.read_check_time()
        self.nTime = len(self.timeList)
        self.read_dimensions()
        
        self.database = {}
        for ensemble_name in self.ensemble_nameList:
            self.database[ensemble_name] = {}

        print("PartMC case created.")
        print("case name: " + self.caseName)
        print("ensemble number: " + str(self.N_ensemble))
        print()

    def clean_database(self):
        del(self.database)
        self.database = {}
        for ensemble_name in self.ensemble_nameList:
            self.database[ensemble_name] = {}


    def add_to_database(self, ensemble_name, varName, dimensions, data, info = {}):
        if varName in self.database[ensemble_name][varName]:
            print("Warning! " + varName + " already exists in the database, which will be covered.")
        self.database[ensemble_name][varName] = (data, dimensions, info)

    def read_from_database(self, ensemble_name, varName):
        assert ensemble_name in self.database, "Error! Couldn't find " + ensemble_name + " in " + self.caseName + "."
        assert varName in self.database[ensemble_name], "Error! Couldn't find " + varName + " in " + self.caseName + ":" + ensemble_name + "."
        data, dimensions, info = self.database[ensemble_name][varName]
        return data, dimensions, info

    def compute_RH_ice(self, ensemble_name = None):
        if ensemble_name is None:
            ensemble_name = self.default_ensmeble_name
        if "RH_ice" in self.database[ensemble_name]:
            RH_ice, dimensions, info = self.database[ensemble_name]["RH_ice"]
        else:
            RH_water, _, _ = self.get_var("relative_humidity", ensemble_name)
            temperature_tList, _, _ = self.get_var("temperature", ensemble_name)
            RH_water = np.array(RH_water)
            temperature_tList = np.array(temperature_tList)
            es = self.compute_saturated_vapor_pressure_water(temperature_tList)
            esi = self.compute_saturated_vapor_pressure_ice(temperature_tList)
            ratio = es/esi
            RH_ice = RH_water * ratio
            dimensions = ("time",)
            info = {
                "units": "1",
            }
            self.database[ensemble_name]["RH_ice"] = (RH_ice, dimensions, info)
        return RH_ice, dimensions, info

    def compute_ice_mixing_ratio(self, ensemble_name = None):

        if ensemble_name is None:
            ensemble_name = self.default_ensemble_name

        if "ice_mixing_ratio" in self.database[ensemble_name]:
            mixing_ratio_tList, dimensions, info = self.database[ensemble_name]["ice_mixing_ratio"]
        else:
            temperature_tList, _, _ = self.get_var("temperature", ensemble_name)
            pressure_tList, _, _ = self.get_var("pressure", ensemble_name)
            aero_num_conc_tList, _, _ = self.get_var("aero_num_conc", ensemble_name)
            part_frozen_tList, _, _ = self.get_var("aero_frozen", ensemble_name)
            particle_mass_tList, _, _ = self.get_var("aero_particle_mass", ensemble_name)

            mixing_ratio_tList = np.full((self.nTime), np.nan)
            for i_time in range(self.nTime):
                particle_mass = particle_mass_tList[i_time].sum(axis = 0)
                air_density = pressure_tList[i_time] / (286 * temperature_tList[i_time])
                freezing_mass = np.sum(part_frozen_tList[i_time] * aero_num_conc_tList[i_time] * particle_mass)
                mixing_ratio = freezing_mass / (air_density * 1) 
                #mixing_ratio *= 1000 # g/g -> g/Kg
                mixing_ratio_tList[i_time] = mixing_ratio

            dimensions = ("time",)
            info = {
                "units": "g/g",
            }
            self.database[ensemble_name]["ice_mixing_ratio"] = (mixing_ratio_tList, dimensions, info)

        return mixing_ratio_tList, dimensions, info
    
    def compute_ice_num_conc(self, ensemble_name = None):

        if ensemble_name is None:
            ensemble_name = self.default_ensemble_name
        if "ice_num_conc" in self.database[ensemble_name]:
            ice_num_conc_tList, dimensions, info = self.database[ensemble_name]["ice_num_conc"]
        else:
            part_frozen_tList, _, _ = self.get_var("aero_frozen", ensemble_name)
            aero_num_conc_tList, _, _ = self.get_var("aero_num_conc", ensemble_name)
            
            ice_num_conc_tList = np.full((self.nTime), np.nan)
            for i_time in range(self.nTime):
                ice_num_conc = np.sum(part_frozen_tList[i_time] * aero_num_conc_tList[i_time])
                ice_num_conc_tList[i_time] = ice_num_conc
            dimensions = ("time", )
            info = {
                "units": "m^-1",
            }
            self.database[ensemble_name]["ice_num_conc"] = (ice_num_conc_tList, dimensions, info)

        return ice_num_conc_tList, dimensions, info


    def compute_ice_ratio(self, ensemble_name = None):
        
        if ensemble_name is None:
            ensemble_name = self.default_ensemble_name

        if "ice_ratio" in self.database[ensemble_name]:
            ice_ratio, dimensions, info = self.database[ensemble_name]["ice_ratio"]
        else:
            part_frozen_tList, _, _ = self.get_var("aero_frozen", ensemble_name)
            aero_num_conc_tList, _, _ = self.get_var("aero_num_conc", ensemble_name)

            ice_ratio_tList = np.full((self.nTime), np.nan)
            for i_time in range(self.nTime):
                ice_ratio = np.sum(part_frozen_tList[i_time] * aero_num_conc_tList[i_time]) / np.sum(aero_num_conc_tList[i_time])
                ice_ratio_tList[i_time] = ice_ratio
            dimensions = ("time",)
            info = {
                #"units": "unitless",
                "units": "1",
            }
            self.database[ensemble_name]["ice_ratio"] = (ice_ratio_tList, dimensions, info)

        return ice_ratio_tList, dimensions, info

    def compute_dry_diameter(self, ensemble_name = None):

        if ensemble_name is None:
            ensemble_name = self.default_ensemble_name

        if "dry_diameter" in self.database[ensemble_name]:
            dry_diameter_tList, dimensions, info = self.database[ensemble_name]["dry_diameter"]
        else:
            aero_particle_mass, dimensions, _ = self.get_var("aero_particle_mass", ensemble_name = ensemble_name)
            aero_density, _, _ = self.get_var("aero_density", ensemble_name = ensemble_name) 
            aero_species_ind, _, aero_species_info = self.get_var("aero_species", ensemble_name = ensemble_name)# - 1

            #_, _, aero_species_info = self.get_var("aero_species", ensemble_name = ensemble_name)
            aero_species_names = aero_species_info["names"].split(",")

            dry_diameter_tList = []
            for i_time in range(self.nTime):
                
                aero_species_ind2name = {ind:name for ind, name in zip(aero_species_ind[i_time] - 1, aero_species_names)}
                aero_species_name2ind = {name:ind for ind, name in zip(aero_species_ind[i_time] - 1, aero_species_names)}
                #diameter = self.aero_particle_diameter(aero_particle_mass[i_time], aero_density[i_time], aero_species_name2ind)
                dry_diameter = self.aero_particle_dry_diameter(aero_particle_mass[i_time], aero_density[i_time], aero_species_name2ind)
                dry_diameter_tList.append(dry_diameter)

            dimensions = ("time", ) + (dimensions[-1], )
            info = {
                "units": "m",
            }
            self.database[ensemble_name]["dry_diameter"] = (dry_diameter_tList, dimensions, info)
        return dry_diameter_tList, dimensions, info

    def convert_to_ice_volume(self, ensemble_name = None):

        if ensemble_name is None:
            ensemble_name = self.default_ensemble_name
        timeList = self.timeList
        aero_mass_list, _, mass_info = self.get_var("aero_particle_mass", ensemble_name)
        aero_density_list, _, _ = self.get_var("aero_density", ensemble_name = ensemble_name)
        aero_species_ind_list, _, aero_species_info = self.get_var("aero_species", ensemble_name = ensemble_name)
        den_ice_list, dim, info = self.get_var("aero_ice_density", ensemble_name)
        frozen_list, dim, info = self.get_var("aero_frozen", ensemble_name)
        aero_volume_list = []
        for i_time in range(self.nTime):
            aero_mass = aero_mass_list[i_time]
            aero_density = aero_density_list[i_time]
            aero_species_ind = aero_species_ind_list[i_time]
            den_ice = den_ice_list[i_time]
            frozen = frozen_list[i_time]
            aero_species_names = aero_species_info["names"].split(",")
            aero_species_ind2name = {ind:name for ind, name in zip(aero_species_ind - 1, aero_species_names)}
            aero_species_name2ind = {name:ind for ind, name in zip(aero_species_ind - 1, aero_species_names)}
            i_water = aero_species_name2ind["H2O"]
            aero_volume = np.full(aero_mass.shape, np.nan)
            for i_spec in range(len(aero_density)):
                if i_spec == i_water:
                    aero_volume[i_spec, :] = np.where(frozen == 1, aero_mass[i_spec, :] / den_ice, aero_mass[i_spec, :] / aero_density[i_spec])
                else:
                    aero_volume[i_spec, :] = aero_mass[i_spec, :] / aero_density[i_spec]
            aero_volume_list.append(aero_volume)
        #return aero_volume_list, aero_mass_list, den_ice_list, i_water
        return aero_volume_list

    def compute_ice_diameter(self, ensemble_name = None):
        if ensemble_name is None:
            ensemble_name = self.default_ensemble_name
        aero_volume_list = self.convert_to_ice_volume(ensemble_name)
        aero_volume_list = np.array(aero_volume_list)
        Dp_ice = (6 / np.pi * aero_volume_list.sum(axis = 1)) ** (1/3.0)
        return Dp_ice

    def compute_diameter(self, ensemble_name = None):

        if ensemble_name is None:
            ensemble_name = self.default_ensemble_name

        if "diameter" in self.database[ensemble_name]:
            diameter_tList, dimensions, info = self.database[ensemble_name]["diameter"]
        else:
            aero_particle_mass, dimensions, _ = self.get_var("aero_particle_mass", ensemble_name = ensemble_name)
            aero_density, _, _ = self.get_var("aero_density", ensemble_name = ensemble_name) 
            aero_species_ind, _, aero_species_info = self.get_var("aero_species", ensemble_name = ensemble_name)# - 1

            #_, _, aero_species_info = self.get_var("aero_species", ensemble_name = ensemble_name)
            aero_species_names = aero_species_info["names"].split(",")

            diameter_tList = []
            for i_time in range(self.nTime):
                
                aero_species_ind2name = {ind:name for ind, name in zip(aero_species_ind[i_time] - 1, aero_species_names)}
                aero_species_name2ind = {name:ind for ind, name in zip(aero_species_ind[i_time] - 1, aero_species_names)}
                diameter = self.aero_particle_diameter(aero_particle_mass[i_time], aero_density[i_time], aero_species_name2ind)
                #dry_diameter = self.aero_particle_dry_diameter(aero_particle_mass[i_time], aero_density[i_time], aero_species_name2ind)
                diameter_tList.append(diameter)

            dimensions = ("time", ) + (dimensions[-1],)
            info = {
                "units": "m",
            }
            self.database[ensemble_name]["diameter"] = (diameter_tList, dimensions, info)
        return diameter_tList, dimensions, info

    def compute_dwetddry(self, ensemble_name = None):
        if ensemble_name is None:
            ensemble_name = self.default_ensemble_name
        if "dwetddry" in self.database[ensemble_name]:
            dwetddry_tList, dimensions, info = self.database[ensemble_name]["dwetddry"]
        else:
            dwet_tList, dimensions, _ = self.compute_diameter(ensemble_name = ensemble_name)
            ddry_tList, _, _ = self.compute_dry_diameter(ensemble_name = ensemble_name)
            dwetddry_tList = []

            for i_time in range(self.nTime):
                dwetddry = dwet_tList[i_time] / ddry_tList[i_time]
                dwetddry_tList.append(dwetddry)
            info = {
                "units": "1",
            }
            self.database[ensemble_name]["dwetddry"] = (dwetddry_tList, dimensions, info)
        return dwetddry_tList, dimensions, info

    def compute_droplets_ratio(self, droplets_threshold = 1.4, ensemble_name = None):

        if ensemble_name is None:
            ensemble_name = self.default_ensemble_name
        if "droplets_ratio" in self.database[ensemble_name]:
            droplets_ratio_tList, dimensions, info = self.database[ensemble_name]["droplets_ratio"]
        else:
            dwetddry_tList, dimensions, info = self.compute_dwetddry(ensemble_name = ensemble_name)
            droplets_ratio_tList = []
            for i_time in range(self.nTime):
                droplets_ratio = np.sum(dwetddry_tList[i_time] >= droplets_threshold) / len(dwetddry_tList[i_time])
                droplets_ratio_tList.append(droplets_ratio)

            droplets_ratio_tList = np.array(droplets_ratio_tList)
            dimensions = ("time",)
            info = {
                "units": "1",
            }
            self.database[ensemble_name]["droplets_ratio"] = droplets_ratio_tList, dimensions, info
        return droplets_ratio_tList, dimensions, info

     

    
    def get_var(self, varName, ensemble_name = None):
        if ensemble_name is None:
            ensemble_name = self.default_ensemble_name
            print("Warning! Undefine ensemble_name, using " + ensemble_name + " as a default value.")
        if varName in self.database[ensemble_name]:
            data_array, dimensions, info = self.database[ensemble_name][varName]
        else:
            print("Couldn't find \"" + varName + "\" from the database: \"" + ensemble_name + "\", search it from the files.")
            data_array = []
            #shape = None
            dimensions = None
            info = None
            for fileName in self.ensemble_fileList[ensemble_name]:
                ncf = nc.Dataset(fileName)
                data = ncf.variables[varName][:].filled(np.nan)
                #len(data.shape) == len(ncf.dimensions
                #if shape is None:
                #    shape = data.shape
                #else:
                #    assert shape == data.shape
                if dimensions is None:
                    dimensions = ncf.variables[varName].dimensions
                else:
                    assert dimensions == ncf.variables[varName].dimensions
                if info is None:
                    info = {}
                    for attr in ncf.variables[varName].ncattrs():
                        info[attr] =  ncf.variables[varName].getncattr(attr)
                #assert len(shape) == len(dimensions)
                data_array.append(data)
            #data_array = np.array(data_array)
            dimensions = ("time",) + dimensions
            self.database[ensemble_name][varName] = (data_array, dimensions, info)
        return data_array, dimensions, info

    def read_check_time(self):
        self.timeList = None
        for ensemble_name in self.ensemble_nameList:
            timeList = []
            for fileName in self.ensemble_fileList[ensemble_name]:
                ncf = nc.Dataset(fileName)
                time = float(ncf.variables["time"][0].filled(np.nan))
                ncf.close()
                timeList.append(time)
            timeList = np.array(timeList)
            if self.timeList is None:
                self.timeList = timeList
            else:
                assert np.max(np.abs(self.timeList - timeList)) < 10e-8

    def read_dimensions(self):
        ncf = nc.Dataset(self.sample_file)
        self.dimensions = {}
        for dimension in ncf.dimensions:
            dimension_name = ncf.dimensions[dimension].name
            dimension_size = ncf.dimensions[dimension].size
            self.dimensions[dimension_name] = dimension_size
        ncf.close()
            

                
    def compute_ensemble_mean_std(self, varName, ensembles = None, getFun = None, getFun_kwargs = {}, return_info = False):

        if ensembles is None:
            ensembles = self.ensemble_nameList
        mean_list = []
        std_list = []
        for i_time in range(self.nTime):
            combine_array = []
            shape_check = None
            for ensemble_name in ensembles:
                if not(varName in self.database[ensemble_name]):
                    if getFun is None:
                        self.get_var(varName, ensemble_name)
                    else:
                        getFun(ensemble_name = ensemble_name, **getFun_kwargs)

                assert varName in self.database[ensemble_name]
                data = self.database[ensemble_name][varName][0][i_time]
                if shape_check is None:
                    shape_check = data.shape
                else:
                    assert shape_check == data.shape
                combine_array.append(data)
            combine_array = np.array(combine_array)
            mean_value = np.nanmean(combine_array, axis = 0)
            std_value = np.nanstd(combine_array, axis = 0)
            mean_list.append(mean_value)
            std_list.append(std_value)
        mean_list = np.array(mean_list)
        std_list = np.array(std_list)
        if return_info:
            info_dict = self.database[ensemble_name][varName][2]
            return mean_list, std_list, info_dict
        else:
            return mean_list, std_list

    def compute_ensemble_max_min(self, varName, ensembles = None, getFun = None, getFun_kwargs = {}, return_info = False):
        if ensembles is None:
            ensembles = self.ensemble_nameList
        max_list = []
        min_list = []
        for i_time in range(self.nTime):
            combine_array = []
            shape_check = None
            for ensemble_name in ensembles:
                if not(varName in self.database[ensemble_name]):
                    if getFun is None:
                        self.get_var(varName, ensemble_name)
                    else:
                        getFun(ensemble_name = ensemble_name, **getFun_kwargs)

                assert varName in self.database[ensemble_name]
                data = self.database[ensemble_name][varName][0][i_time]
                if shape_check is None:
                    shape_check = data.shape
                else:
                    assert shape_check == data.shape
                combine_array.append(data)
            combine_array = np.array(combine_array)
            max_value = np.nanmax(combine_array, axis = 0)
            min_value = np.nanmin(combine_array, axis = 0)
            max_list.append(max_value)
            min_list.append(min_value)
        max_list = np.array(max_list)
        min_list = np.array(min_list)
        if return_info:
            info_dict = self.database[ensemble_name][varName][2]
            return max_list, min_list, info_dict
        else:
            return max_list, min_list
                
            
            
    def vol2rad(self, vol):
        return (3 * vol / 4 / np.pi) ** (1/3.0)

    def aero_particle_volume(self, aero_particle_mass, aero_density, aero_species_name2ind):
        vol = 0
        for spec_name in aero_species_name2ind:
            vol += aero_particle_mass[aero_species_name2ind[spec_name], :] / aero_density[aero_species_name2ind[spec_name]]
        return vol

    def aero_particle_dry_volume(self, aero_particle_mass, aero_density, aero_species_name2ind):
        vol = 0
        for spec_name in aero_species_name2ind:
            if spec_name == "H2O":
                continue
            vol += aero_particle_mass[aero_species_name2ind[spec_name], :] / aero_density[aero_species_name2ind[spec_name]]
        return vol

    def aero_particle_diameter(self, aero_particle_mass, aero_density, aero_species_name2ind):
        vol = self.aero_particle_volume(aero_particle_mass, aero_density, aero_species_name2ind)
        radius = self.vol2rad(vol)
        diameter = 2 * radius
        return diameter

    def aero_particle_dry_diameter(self, aero_particle_mass, aero_density, aero_species_name2ind):
        vol = self.aero_particle_dry_volume(aero_particle_mass, aero_density, aero_species_name2ind)
        radius = self.vol2rad(vol)
        diameter = 2 * radius
        return diameter

    def compute_saturated_vapor_pressure_water(self, T):
        tmp = 54.842763 \
            - 6763.22 / T \
            - 4.210 * np.log(T) \
            + 0.000367 * T \
            + np.tanh( 0.0415 * (T - 218.8)) \
                * (53.878 - 1331.22 / T - 9.44523 * np.log(T) + 0.014025 * T)
        es = np.exp(tmp)
        return es

    def compute_saturated_vapor_pressure_ice(self, T):
        tmp = 9.550426 \
            - 5723.265 / T \
            + 3.53068 * np.log(T) \
            - 0.00728332 * T
        esi = np.exp(tmp)
        return esi


    def compute_histogram(self, time, ensemble_name = None, bins = 100):

        if ensemble_name is None:
            ensemble_name = self.default_ensemble_name

        i_time = np.where(time - self.timeList >= 0 )[0][-1]

        aero_particle_mass, _, _ = self.get_var("aero_particle_mass", ensemble_name = ensemble_name) 
        aero_particle_mass = aero_particle_mass[i_time]# (aero_species, aero_particle)

        aero_density, _, _ = self.get_var("aero_density", ensemble_name = ensemble_name) 
        aero_density = aero_density[i_time] #(aero_species) 

        aero_num_conc, _, _ = self.get_var("aero_num_conc", ensemble_name = ensemble_name) 
        aero_num_conc = aero_num_conc[i_time] # (aero_particle)

        aero_particle_vol = (aero_particle_mass.T / aero_density).T

        aero_species_ind, _, _ = self.get_var("aero_species", ensemble_name = ensemble_name)# - 1
        aero_species_ind = aero_species_ind[i_time] - 1

        _, _, aero_species_info = self.get_var("aero_species", ensemble_name = ensemble_name)
        aero_species_names = aero_species_info["names"].split(",")

        aero_species_ind2name = {ind:name for ind, name in zip(aero_species_ind, aero_species_names)}
        aero_species_name2ind = {name:ind for ind, name in zip(aero_species_ind, aero_species_names)}

        diameter = self.aero_particle_diameter(aero_particle_mass, aero_density, aero_species_name2ind)
        dry_diameter = self.aero_particle_dry_diameter(aero_particle_mass, aero_density, aero_species_name2ind)

        part_frozen_tList, _, _ = self.get_var("aero_frozen", ensemble_name)
        part_frozen = part_frozen_tList[i_time]

        diameter *= 10**6
        dry_diameter *= 10**6

        hists_dict = {}
    
        hist, bin_edges = np.histogram(np.log10(diameter), weights=aero_num_conc, bins = bins)
        hists_dict["diameter_log10"] = (hist, bin_edges)
        hist, bin_edges = np.histogram(np.log10(dry_diameter), weights=aero_num_conc, bins = bins)
        hists_dict["dry_diameter_log10"] = (hist, bin_edges)
        hist, bin_edges = np.histogram(np.log10(diameter), weights = aero_num_conc * part_frozen, bins = bins)
        hists_dict["ice_diameter_log10"] = (hist, bin_edges)
        hist, bin_edges = np.histogram(np.log10(dry_diameter), weights = aero_num_conc * part_frozen, bins = bins)
        hists_dict["ice_dry_diameter_log10"] = (hist, bin_edges)

        hist, bin_edges = np.histogram(diameter, weights=aero_num_conc, bins = bins)
        hists_dict["diameter"] = (hist, bin_edges)
        hist, bin_edges = np.histogram(dry_diameter, weights=aero_num_conc, bins = bins)
        hists_dict["dry_diameter"] = (hist, bin_edges)
        hist, bin_edges = np.histogram(diameter, weights=aero_num_conc * part_frozen, bins = bins)
        hists_dict["ice_diameter"] = (hist, bin_edges)
        hist, bin_edges = np.histogram(dry_diameter, weights=aero_num_conc * part_frozen, bins = bins)
        hists_dict["ice_dry_diameter"] = (hist, bin_edges)

        return hists_dict






            




