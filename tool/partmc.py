#!/usr/bin/env python
# Copyright (C) 2007-2009 Matthew West
# Licensed under the GNU General Public License version 2 or (at your
# option) any later version. See the file COPYING for details.

import os
import sys
import re
import textwrap
import math
import numpy
import random
import scipy
import Scientific.IO.NetCDF

class aero_data_t(object):

    """Stores the physical constants describing aerosol species. All
    data atributes are 1D arrays, with one entry per species. Thus
    names[i] and densities[i] give the name and density of the i-th
    aerosol species, respectively.

    All data attribute arrays have the same length, so the number of
    species can be found, for example, with:
    
    >>> n_aero_species = len(aero_data.names)

    The data attributes are:
    
    names - names of the aerosol species (strings)
    mosaic_indices - species index numbers in MOSAIC (integers)
    densities - density of each species (kg/m^3)
    num_ions - number of ions for dissociation (integers)
    solubilities - solubility factors (dimensionless)
    molec_weights - molecular weights (kg/mole)
    kappas - hydroscopicity parameters (dimensionless)"""

    def __init__(self, ncf = None, n_species = None):
        
        """ Creates an aero_data_t object. If the ncf parameter is
        passed a NetCDFFile object output from PartMC then the
        aero_data information will be loaded from the file. Otherwise
        the n_species parameter must be passed and an empty
        aero_data_t object will be created. For example:

        >>> ncf = Scientific.IO.NetCDF.NetCDFFile('filename.nc')
        >>> aero_data = partmc.aero_data_t(ncf)

        or

        >>> aero_data = partmc.aero_data_t(n_species = 10)"""

        if ncf is not None:
            if "aero_species" not in ncf.variables.keys():
                raise Exception("aero_species variable not found in NetCDF file")
            if "names" not in dir(ncf.variables["aero_species"]):
                raise Exception("aero_species variable does not have 'names' attribute")
            self.names = ncf.variables["aero_species"].names.split(",")

            for (ncf_var, self_var) in [
                ("aero_mosaic_index", "mosaic_indices"),
                ("aero_density", "densities"),
                ("aero_num_ions", "num_ions"),
                ("aero_solubility", "solubilities"),
                ("aero_molec_weight", "molec_weights"),
                ("aero_kappa", "kappas"),
                ]:
                if ncf_var not in ncf.variables.keys():
                    raise Exception("%s variable not found in NetCDF file" % ncf_var)
                self.__dict__[self_var] = numpy.asarray(ncf.variables[ncf_var].getValue())
                
        if n_species is not None:
            if ncf is not None:
                raise Exception("ncf and n_species arguments cannot both be specified")
            self.names = ["" for i in range(n_species)]
            self.mosaic_indices = numpy.zeros([n_species], float)
            self.densities = numpy.zeros([n_species], float)
            self.num_ions = numpy.zeros([n_species], int)
            self.solubilities = numpy.zeros([n_species], float)
            self.molec_weights = numpy.zeros([n_species], float)
            self.kappas = numpy.zeros([n_species], float)

        if ncf is None and n_species is None:
            raise Exception("either ncf or n_species parameter must be specified")

class env_state_t(object):

    """Stores the environment state. All data attributes are scalars,
    giving the state of the variable at a single point in time and
    space. The data attributes are:

    temperature - air temperature (K)
    relative_humidity - air relative humidity (dimensionless)
    pressure - air pressure (Pa)
    longitude - location longitude (degrees)
    latitude - location latitude (degrees)
    altitude - location altitude (m)
    start_time_of_day - time of day when simulation started (s)
    start_day_of_year - day number of year when the simulation
        started (integer)
    elapsed_time - elapsed simulation time since the simulation start (s)
    height - mixing layer height (m)
    """
    
    def __init__(self, ncf = None):
        
        """ Creates an env_state_t object. If the ncf parameter is
        passed a NetCDFFile object output from PartMC then the
        env_state information will be loaded from the file. Otherwise
        an empty env_state_t object will be created. For example:

        >>> ncf = Scientific.IO.NetCDF.NetCDFFile('filename.nc')
        >>> env_state = partmc.env_state_t(ncf)

        or

        >>> env_state = partmc.env_state_t()"""

        if ncf is not None:
            for (ncf_var, self_var, var_type) in [
                ("temperature", "temperature", float),
                ("relative_humidity", "relative_humidity", float),
                ("pressure", "pressure", float),
                ("longitude", "longitude", float),
                ("latitude", "latitude", float),
                ("altitude", "altitude", float),
                ("start_time_of_day", "start_time_of_day", float),
                ("start_day_of_year", "start_day_of_year", int),
                ("elapsed_time", "elapsed_time", float),
                ("height", "height", float),
                ]:
                if ncf_var not in ncf.variables.keys():
                    raise Exception("%s variable not found in NetCDF file" % ncf_var)
                self.__dict__[self_var] = var_type(ncf.variables[ncf_var].getValue())
        else:
            self.temperature = 0.0
            self.relative_humidity = 0.0
            self.pressure = 0.0
            self.longitude = 0.0
            self.latitude = 0.0
            self.altitude = 0.0
            self.start_time_of_day = 0.0
            self.start_day_of_year = 0
            self.elapsed_time = 0.0
            self.height = 0.0

    def A(self, constants):

        """ Computes the A parameter (for condensation
        calculations). For example:

        >>> constants = partmc.constants_t('path_to_src/constants.f90')
        >>> A = env_state.A(constants)"""
        
        return (4.0 * constants.water_surf_eng * constants.water_molec_weight
                / (constants.univ_gas_constants * self.temperature *
                   constants.water_density))

class gas_data_t(object):

    """Stores the physical constants describing gas species. All data
    atributes are 1D arrays, with one entry per species. Thus names[i]
    and molec_weights[i] give the name and molecular weight of the
    i-th gas species, respectively.

    All data attribute arrays have the same length, so the number of
    species can be found, for example, with:
    
    >>> n_gas_species = len(gas_data.names)

    The data attributes are:
    
    names - names of the gas species (strings)
    mosaic_indices - species index numbers in MOSAIC (integers)
    molec_weights - molecular weights (kg/mole)"""
    
    def __init__(self, ncf = None, n_species = None):
        
        """ Creates a gas_data_t object. If the ncf parameter is
        passed a NetCDFFile object output from PartMC then the
        gas_data information will be loaded from the file. Otherwise
        the n_species parameter must be passed and an empty gas_data_t
        object will be created. For example:

        >>> ncf = Scientific.IO.NetCDF.NetCDFFile('filename.nc')
        >>> gas_data = partmc.gas_data_t(ncf)

        or

        >>> gas_data = partmc.gas_data_t(n_species = 60)"""

        if ncf is not None:
            if "gas_species" not in ncf.variables.keys():
                raise Exception("gas_species variable not found in NetCDF file")
            if "names" not in dir(ncf.variables["gas_species"]):
                raise Exception("gas_species variable does not have 'names' attribute")
            self.names = ncf.variables["gas_species"].names.split(",")

            for (ncf_var, self_var) in [
                ("gas_mosaic_index", "mosaic_indices"),
                ("gas_molec_weight", "molec_weights"),
                ]:
                if ncf_var not in ncf.variables.keys():
                    raise Exception("%s variable not found in NetCDF file" % ncf_var)
                self.__dict__[self_var] = numpy.asarray(ncf.variables[ncf_var].getValue())
                
        if n_species is not None:
            if ncf is not None:
                raise Exception("ncf and n_species arguments cannot both be specified")
            self.names = ["" for i in range(n_species)]
            self.mosaic_indices = numpy.zeros([n_species], float)
            self.molec_weights = numpy.zeros([n_species], float)

        if ncf is None and n_species is None:
            raise Exception("either ncf or n_species parameter must be specified")

class gas_state_t(object):

    """Stores the gas state at a single point in time and space. The
    data attributes are:

    raw_mixing_ratios - mixing ratios of each gas species (dimensionless)
    gas_data - object of type gas_data_t with per-species physical data

    The mixing ratio of a gas species should be obtained by using the
    mixing_ratio() method, not by accessing the raw_mixing_ratios
    attribute directly, to avoid confusion with species numbers. For
    example:

    >>> ozone = gas_state.mixing_ratio('O3')"""
    
    def __init__(self, ncf = None, gas_data = None):
        
        """Creates a gas_state_t object. If the ncf parameter is
        passed a NetCDFFile object output from PartMC then the
        gas_state information will be loaded from the file. Otherwise
        the gas_data parameter must be passed and an empty gas_state
        will be created with the appropriate number of species. For
        example:

        >>> ncf = Scientific.IO.NetCDF.NetCDFFile('filename.nc')
        >>> gas_state = partmc.gas_state_t(ncf)

        or

        >>> ncf = Scientific.IO.NetCDF.NetCDFFile('filename.nc')
        >>> gas_data = partmc.gas_data_t(ncf)
        >>> gas_state = partmc.gas_state_t(gas_data = gas_data)"""

        if ncf is not None:
            if gas_data is not None:
                raise Exception("cannot specify both ncf and gas_data parameters")
            self.gas_data = gas_data_t(ncf)
            if "gas_mixing_ratio" not in ncf.variables.keys():
                raise Exception("gas_mixing_ratio variable not found in NetCDF file")
            self.raw_mixing_ratios = numpy.asarray(ncf.variables["gas_mixing_ratio"])
        else:
            if gas_data is None:
                raise Exception("must specify either ncf or gas_data parameters")
            self.gas_data = gas_data
            self.raw_mixing_ratios = numpy.zeros([len(gas_data.names)], float)

    def mixing_ratio(self, species):

        """Return the mixing ratio of the given species name. For
        example, the mixing ratio of ozone is given by:

        >>> ozone = gas_state.mixing_ratio('O3')"""
        
        if species not in self.gas_data.name:
            raise Exception("unknown species: %s" % species)
        index = self.gas_data.name.index(species)
        return self.raw_mixing_ratio[index]

class aero_particle_array_t(object):

    """Stores the particles at a single point in time and space. All
    data attributes (except aero_data) are 1D or 2D arrays with one
    entry (or column if 2D) per aerosol particle. The data attributes
    are (assuming S species and N particles):

    aero_data - object of type aero_data_t with per-species physical data
    raw_masses - S x N array of the mass (kg) of species i in particle j
    n_orig_parts - length N array with the number of original particles
        that coagulated to form each particle
    absorb_cross_sects - length N array with absorbion cross section (m^2
        of each particle
    scatter_cross_sects - length N array with scattering cross section
        (m^2) of each particle
    asymmetries - length N array with asymmetry parameter of each particle
    refract_shell_reals - length N array with the real part of the
        refractive index of the shell of each particle
    refract_shell_imags - length N array with the imaginary part of the
        refractive index of the shell of each particle
    refract_core_reals - length N array with the real part of the
        refractive index of the core of each particle
    refract_core_imags - length N array with the imaginary part of the
        refractive index of the core of each particle
    core_vols - length N array with the volume of the core (m^3) of each
        particle
    water_hyst_legs - length N array with which leg of the water
        hysteresis curve the particle is on
    comp_vols - length N array with the computational volume (m^3)
        associated with each particle
    ids - length N array with the ID number of each particle
    least_create_times - length N array with the earliest creation
        time (s) of any component of each particle
    greatest_create_times - length N array with the latest creation
        time (s) of any component of each particle

    Most of these attributes should be accessed directly, except for
    the raw_masses array. To prevent confusion, this should be
    accessed by the masses() method.

    To get the number of particles we can use, for example:
    >>> n_particles = len(aero_particle_array.ids)"""

    def __init__(self, ncf = None, n_particles = None, aero_data = None,
                 include_ids = None, exclude_ids = None):

        """Creates an aero_particle_array_t object. If the ncf
        parameter is passed a NetCDFFile object output from PartMC
        then the aero_particle_array information will be loaded from
        the file. The include_ids and exclude_ids parameters can be
        given lists of particle IDs to include or exclude when loading
        from the file.

        If the ncf parameter is not specified then the aero_data and
        n_particles parameters must be given and an empty
        aero_particle_array_t will be created with the appropriate
        number of particles.

        For example:

        >>> ncf = Scientific.IO.NetCDF.NetCDFFile('filename.nc')
        >>> aero_particle_array = partmc.aero_particle_array_t(ncf)

        or

        >>> ncf = Scientific.IO.NetCDF.NetCDFFile('filename.nc')
        >>> aero_particle_array = partmc.aero_particle_array_t(ncf,
                include_ids = [45, 182, 7281])

        or

        >>> ncf = Scientific.IO.NetCDF.NetCDFFile('filename.nc')
        >>> aero_data = partmc.aero_data_t(ncf)
        >>> aero_particle_array \\
                = partmc.aero_particle_array_t(n_particles = 100,
                                               aero_data = aero_data)"""
        
        if ncf == None:
            if aero_data is None:
                raise Exception("must pass aero_data when ncf parameter is not given")
            if n_particles is None:
                raise Exception("must pass n_particles when ncf parameter is not given")
            if include_ids is not None:
                raise Exception("cannot provide include_ids without the ncf parameter")
            if exclude_ids is not None:
                raise Exception("cannot provide include_ids without the ncf parameter")
            self.aero_data = aero_data
            self.raw_masses = zeros([len(aero_data.name), n_particles])
            self.n_orig_parts = zeros(n_particles, int)
            self.absorb_cross_sects = zeros(n_particles)
            self.scatter_cross_sects = zeros(n_particles)
            self.asymmetries = zeros(n_particles)
            self.refract_shell_reals = zeros(n_particles)
            self.refract_shell_imags = zeros(n_particles)
            self.refract_core_reals = zeros(n_particles)
            self.refract_core_imags = zeros(n_particles)
            self.core_vols = zeros(n_particles)
            self.water_hyst_legs = zeros(n_particles, int)
            self.comp_vols = zeros(n_particles)
            self.ids = zeros(n_particles, int)
            self.least_create_times = zeros(n_particles)
            self.greatest_create_times = zeros(n_particles)
            return

        if aero_data is not None:
            raise Exception("cannot specify aero_data when ncf parameter is given")
        if n_particles is not None:
            raise Exception("cannot specify n_particles when ncf parameter is given")

        self.aero_data = aero_data_t(ncf)

        for (ncf_var, self_var) in [
            ("aero_particle_mass", "raw_masses"),
            ("aero_n_orig_part", "n_orig_parts"),
            ("aero_absorb_cross_sect", "absorb_cross_sects"),
            ("aero_scatter_cross_sect", "scatter_cross_sects"),
            ("aero_asymmetry", "asymmetries"),
            ("aero_refract_shell_real", "refract_shell_reals"),
            ("aero_refract_shell_imag", "refract_shell_imags"),
            ("aero_refract_core_real", "refract_core_reals"),
            ("aero_refract_core_imag", "refract_core_imags"),
            ("aero_core_vol", "core_vols"),
            ("aero_water_hyst_leg", "water_hyst_legs"),
            ("aero_comp_vol", "comp_vols"),
            ("aero_id", "ids"),
            ("aero_least_create_time", "least_create_times"),
            ("aero_greatest_create_time", "greatest_create_times"),
            ]:
            if ncf_var not in ncf.variables.keys():
                raise Exception("%s variable not found in NetCDF file" % ncf_var)
            self.__dict__[self_var] = asarray(ncf.variables[ncf_var].getValue())

        if include_ids != None or exclude_ids != None:
            keep_indexes = [i for i in range(size(self.ids)) \
                            if (include_ids != None and self.ids[i] in include_ids) \
                            or (exclude_ids != None and self.ids[i] not in exclude_ids)]
            self.masses = self.masses[:, keep_indexes]
            self.n_orig_parts = self.n_orig_parts[keep_indexes]
            self.absorb_cross_sects = self.absorb_cross_sects[keep_indexes]
            self.scatter_cross_sects = self.scatter_cross_sects[keep_indexes]
            self.asymmetries = self.asymmetries[keep_indexes]
            self.refract_shell_reals = self.refract_shell_reals[keep_indexes]
            self.refract_shell_imags = self.refract_shell_imags[keep_indexes]
            self.refract_core_reals = self.refract_core_reals[keep_indexes]
            self.refract_core_imags = self.refract_core_imags[keep_indexes]
            self.core_vols = self.core_vols[keep_indexes]
            self.water_hyst_legs = self.water_hyst_legs[keep_indexes]
            self.comp_vols = self.comp_vols[keep_indexes]
            self.ids = self.ids[keep_indexes]
            self.least_create_times = self.least_create_times[keep_indexes]
            self.greatest_create_times = self.greatest_create_times[keep_indexes]

    def sum_masses_weighted(self, include = None, exclude = None,
                            species_weights = None):

        """Computes the weighted sum of the component masses of each
        particle, including and excluding the given species lists by
        name. This function should generally not be called directly,
        but instead the masses(), volumes(), etc. functions should be
        used.

        Example usage:
        >>> dry_volumes = aero_particle_array.sum_masses_weighted(
                exclude = ['H2O'], species_weights
                = 1 / aero_particle_array.aero_data.density)"""
        
        if include != None:
            for species in include:
                if species not in self.aero_data.name:
                    raise Exception("unknown species: %s" % species)
            species_list = set(include)
        else:
            species_list = set(self.aero_data.name)
        if exclude != None:
            for species in exclude:
                if species not in self.aero_data.name:
                    raise Exception("unknown species: %s" % species)
            species_list -= set(exclude)
        species_list = list(species_list)
        if len(species_list) == 0:
            raise Exception("no species left to sum over")
        index = self.aero_data.name.index(species_list[0])
        if species_weights != None:
            val = self.masses[index,:].copy() * species_weights[index]
        else:
            val = array(self.masses[index,:].copy())
        for i in range(len(species_list) - 1):
            index = self.aero_data.name.index(species_list[i + 1])
            if species_weights != None:
                val += self.masses[index,:] * species_weights[index]
            else:
                val += array(self.masses[index,:])
        return val
    
    def masses(self, include = None, exclude = None):

        """Return the total mass (kg) of each particle as an array,
        including or excluding the given species by name. Examples:

        >>> total_masses = aero_particle_array.masses()
        >>> water_masses = aero_particle_array.masses(include = ['H2O'])
        >>> dry_masses = aero_particle_array.masses(exclude = ['H2O'])
        >>> carbon_masses = aero_particle_array.masses(include = ['BC',
                'OC'])"""
        
        return self.sum_masses_weighted(include = include, exclude = exclude)

    def volumes(self, include = None, exclude = None):

        """Return the total volume (m^3) of each particle as an array,
        including or excluding the given species by name. See the
        masses() method for examples of usage."""
        
        species_weights = 1.0 / array(self.aero_data.density)
        return self.sum_masses_weighted(include = include, exclude = exclude,
                                        species_weights = species_weights)

    def moles(self, include = None, exclude = None):

        """Return the total moles (dimensionless) in each particle as
        an array, including or excluding the given species by
        name. See the masses() method for examples of usage."""
        
        species_weights = self.aero_data.molec_weight \
                          / self.aero_data.density
        return self.sum_masses_weighted(include = include, exclude = exclude,
                                        species_weights = species_weights)

    def radii(self):

        """Return the radius (m) of each particle as an array."""
        
        return (self.volume() * 3.0/4.0 / math.pi)**(1.0/3.0)

    def dry_radii(self):

        """Return the dry radius (m) of each particle as an array."""
        
        return (self.volume(exclude = ["H2O"]) * 3.0/4.0 / math.pi)**(1.0/3.0)

    def diameters(self):

        """Return the diameter (m) of each particle as an array."""
        
        return 2.0 * self.radius()

    def dry_diameters(self):

        """Return the dry diameter (m) of each particle as an array."""
        
        return 2.0 * self.dry_radius()

    def surface_areas(self):

        """Return the surface area (m^2) of each particle as an array."""
        
        return 4.0 * math.pi * self.radius()**2

    def dry_surface_areas(self):

        """Return the dry surface area (m^2) of each particle as an array."""
        
        return 4.0 * math.pi * self.dry_radius()**2

    def kappas(self):

        """Return the total kappa (dimensionless hydroscopicity
        parameter) of each particle as an array. This is computed as a
        volume-weighted sum of the per-species kappa values. Each
        species kappa can be specified directly in aero_data.kappa, or
        if this is zero it is calculated from aero_data.num_ions."""
        
        if "H2O" not in self.aero_data.names:
            raise Exception("unable to find water species index by name 'H2O'")
        i_water = self.aero_data.name.index("H2O")
        M_w = self.aero_data.molec_weight[i_water]
        rho_w = self.aero_data.density[i_water]
        species_weights = zeros([len(self.aero_data.name)])
        for i_spec in range(size(species_weights)):
            if i_spec == self.aero_data.name == "H2O":
                continue
            if self.aero_data.num_ions[i_spec] > 0:
                if self.aero_data.kappa[i_spec] != 0:
                    raise Exception("species has nonzero num_ions and kappa: %s" % self.name[i_spec])
                M_a = self.aero_data.molec_weight[i_spec]
                rho_a = self.aero_data.density[i_spec]
                species_weights[i_spec] = M_w * rho_a / (M_a * rho_w) \
                                          * self.aero_data.num_ions[i_spec]
            else:
                species_weights[i_spec] = self.aero_data.kappa[i_spec]
        species_weights /= self.aero_data.density
        volume_kappa = self.sum_masses_weighted(exclude = ["H2O"],
                                                species_weights = species_weights)
        dry_volume = self.volume(exclude = ["H2O"])
        return volume_kappa / dry_volume

    def critical_rel_humids_approx(self, env_state, constants):

        """Compute the critical relative humidities (dimensionless) of
        each particle as an array, using a fast approximate method."""
        
        A = env_state.A(constants)
        C = sqrt(4.0 * A**3 / 27.0)
        dry_diameters = self.dry_diameters()
        kappas = self.solute_kappas()
        S = C / sqrt(kappas * dry_diameters**3) + 1.0
        return S

    def critical_rel_humids(self, env_state, constants):

        """Compute the critical relative humidities (dimensionless) of
        each particle as an array."""
        
        kappas = self.solute_kappas()
        dry_diameters = self.dry_diameters()
        return critical_rel_humids(env_state, constants, kappas, dry_diameters)

    def critical_diameters(self, env_state, constants):

        """Compute the critical diameters (m) of each particle as an
        array."""
        
        kappas = self.solute_kappas()
        dry_diameters = self.dry_diameters()
        return critical_diameters(env_state, constants, kappas, dry_diameters)

    def bin_average(self, diameter_axis, dry_diameter = True):

        """Return a new aero_particle_array_t object with one particle
        per grid cell which is an average of all the original
        particles within that grid cell (by diameter)."""
        
        averaged_particles = aero_particle_array_t(n_particles = diameter_axis.n_bin,
                                                   aero_data = self.aero_data)
        if dry_diameter:
            diameter = self.dry_diameter()
        else:
            diameter = self.diameter()
        diameter_bin = diameter_axis.find(diameter)
        num_conc = zeros(diameter_axis.n_bin)
        masses_conc = zeros([self.masses.shape[0], diameter_axis.n_bin])
        for i in range(self.n_particles):
            b = diameter_bin[i]
            if diameter_axis.valid_bin(b):
                num_conc[b] += 1.0 / self.comp_vols[i]
                masses_conc[:,b] += self.masses[:,i] / self.comp_vols[i]
        for b in range(averaged_particles.n_particles):
            averaged_particles.comp_vols[b] = 1.0 / num_conc[b]
            averaged_particles.masses[:,b] = masses_conc[:,b] * averaged_particles.comp_vols[b]
        return averaged_particles

def critical_rel_humids(env_state, constants, kappas, dry_diameters):

    """Compute the critical relative humidity (dimensionless) for each
    kappa and dry_diameter.

    The kappas and dry_diameters parameters should be 1D arrays of the
    same length N, and the return value will be a 1D array of length N
    with each entry i being the critical RH for a particle with
    average kappa given by kappas[i] and a dry diameter of
    dry_diameters[i].

    The env_state and constants parameters should be objects of type
    env_state_t and constants_t, respectively.

    Example:
    >>> ncf = Scientific.IO.NetCDF.NetCDFFile('filename.nc')
    >>> env_state = partmc.env_state_t(ncf)
    >>> constants = partmc.constants_t(ncf)
    >>> kappas = numpy.array([0.5, 0.2])
    >>> dry_diameters = numpy.array([1e-8, 5e-8])
    >>> crit_rhs = partmc.critical_rel_humids(env_state, constants,
            kappas, dry_diameters)"""
    
    A = env_state.A(constants)
    critical_diameters = critical_diameters(env_state, constants, kappas, dry_diameters)
    return (critical_diameters**3 - dry_diameters**3) \
        / (critical_diameters**3 - dry_diameters**3 * (1 - kappas)) \
        * exp(A / critical_diameters)

def critical_diameters(env_state, constants, kappas, dry_diameters):

    """Compute the critical diameters (m) for each kappa and
    dry_diameter.

    The kappas and dry_diameters parameters should be 1D arrays of the
    same length N, and the return value will be a 1D array of length N
    with each entry i being the critical diameter for a particle with
    average kappa given by kappas[i] and a dry diameter of
    dry_diameters[i].

    The env_state and constants parameters should be objects of type
    env_state_t and constants_t, respectively.

    Example:
    >>> ncf = Scientific.IO.NetCDF.NetCDFFile('filename.nc')
    >>> env_state = partmc.env_state_t(ncf)
    >>> constants = partmc.constants_t(ncf)
    >>> kappas = numpy.array([0.5, 0.2])
    >>> dry_diameters = numpy.array([1e-8, 5e-8])
    >>> crit_diams = partmc.critical_diameters(env_state, constants,
            kappas, dry_diameters)"""
    
    A = env_state.A(constants)
    c4 = - 3.0 * dry_diameters**3 * kappas / A
    c3 = - dry_diameters**3 * (2.0 - kappas)
    c0 = dry_diameters**6 * (1.0 - kappas)
    dc = zeros_like(d1)
    for i in range(len(kappas)):
        def f(d):
            return d**6 + c4[i] * d**4 + c3[i] * d**3 + c0[i]
        d1 = dry_diameters[i]
        if not (f(d1) < 0):
            raise Exception("initialization failure for d1")
        d2 = 2 * d1
        for iteration in range(100):
            if f(d2) > 0:
                break
            d2 *= 2
        else:
            raise Exception("intialization failure for d2")
        dc[i] = scipy.optimize.brentq(f, d1, d2)
    return dc

class aero_removed_info_t(object):

    """Stores information about the particles removed from the aerosol
    particle population at a single point in time and space. All data
    attributes are arrays with one entry per removed aerosol
    particle. The data attributes are:

    ids - particle ID numbers of the removed particles
    actions - action code of each removed particle
    other_ids - associated particle ID numbers for the removals (or 0
        if there is no associated ID)

    The action codes are integers with the following values:
    aero_removed_info_t.AERO_INFO_NONE - No information.
    aero_removed_info_t.AERO_INFO_DILUTION - The particle was removed
        by diluting out of the parcel.
    aero_removed_info_t.AERO_INFO_COAG - The particle was removed due
        to coagulating with another particle. The ID of the particle
        coagulated with is given in the other_ids array.
    aero_removed_info_t.AERO_INFO_HALVED - The particle was removed
        due to a halving of the particle population.

    To obtain the number of removed particles we can use, for example:
    >>> n_removed_particles = len(aero_removed_info.ids)

    Example:
    >>> ncf = Scientific.IO.NetCDF.NetCDFFile('filename.nc')
    >>> aero_removed_info = partmc.aero_removed_info_t(ncf)
    >>> n_removed_particles = len(aero_removed_info.ids)
    >>> for i in range(n_removed_particles):
    >>>     print 'removed particle %d' % i
    >>>     print '    id = %d' % aero_removed_info.ids[i]
    >>>     action = aero_removed_info.actions[i]
    >>>     if action == aero_removed_info.AERO_INFO_NONE:
    >>>         print '    no information'
    >>>     elif action == aero_removed_info.AERO_INFO_DILUTION:
    >>>         print '    removed due to dilution'
    >>>     elif action == aero_removed_info.AERO_INFO_COAG:
    >>>         print '    removed due to coagulation' \\
    >>>             ' with particle id %d' \\
    >>>             % aero_removed_info.other_ids[i]
    >>>     elif action == aero_removed_info.AERO_INFO_HALVING:
    >>>         print '    removed due to halving'"""

    AERO_INFO_NONE = 0
    AERO_INFO_DILUTION = 1
    AERO_INFO_COAG = 2
    AERO_INFO_HALVED = 3

    def __init__(self, ncf):

        """Creates an aero_remove_info_t object by reading data from a
        NetCDF file output from PartMC.
        
        For example:
        >>> ncf = Scientific.IO.NetCDF.NetCDFFile('filename.nc')
        >>> aero_removed_info = partmc.aero_removed_info_t(ncf)"""
        
        for (ncf_var, self_var) in [
            ("aero_removed_id", "ids"),
            ("aero_removed_action", "actions"),
            ("aero_removed_other_id", "other_ids"),
            ]:
            if ncf_var not in ncf.variables.keys():
                raise Exception("%s variable not found in NetCDF file" % ncf_var)
            self.__dict__[self_var] = asarray(ncf.variables[ncf_var].getValue())

        if (len(self.aero_removed_id) == 1) and (self.aero_removed_id[0] == 0):
            self.id = array([],'int32')
            self.action = array([],'int32')
            self.other_id = array([],'int32')

class grid(object):

    """Base class for 1D grids. See linear_grid and log_grid for
    specific grid types."""

    def __init__(self):

        """Do not call this. Instead call linear_grid() or log_grid()
        to make a specific type of grid."""
        
        raise NotImplementedError()

    def find_clipped(self, values):

        """Find the bins for each entry of values, clipped to
        [0, n_bin - 1].

        The parameter values should be a 1D array of values to locate
        within the grid, and the return value is a 1D array of
        integers between 0 (the first bin) and n_bin - 1 (the last
        bin). If a value is below the first bin then 0 is returned
        while if it is above the last bin then n_bin - 1 is
        returned."""
        
        indices = self.find(values)
        indices = indices.clip(0, self.n_bin - 1)
        return indices

    def find_clipped_outer(self, values):

        """Find the bins for each entry of values, clipped to
        [-1, n_bin].

        The parameter values should be a 1D array of values to locate
        within the grid, and the return value is a 1D array of
        integers between -1 and n_bin. If each value is within a bin
        then a number within [0, n_bin - 1] is returned. If a value is
        below the first bin then -1 is returned while if it is above
        the last bin then n_bin is returned."""
        
        indices = self.find(values)
        indices = indices.clip(-1, self.n_bin)
        return indices

    def closest_edge(self, value):

        """Find the closest bin edge to the given value.

        Value should be a single scalar and the return value is an
        integer between 0 and n_bin giving the number of the edge
        closest to value."""
        
        i = self.find_clipped(value)
        lower_edge = self.edge(i)
        upper_edge = self.edge(i + 1)
        if abs(value - lower_edge) < abs(value - upper_edge):
            return i
        else:
            return i + 1

    def edges(self):

        """Return a length (n_bin + 1) array of the bin edges in the
        grid."""
        
        return array([self.edge(i) for i in range(self.n_bin + 1)])

    def centers(self):

        """Return a length n_bin array of the bin centers in the
        grid."""
        
        return array([self.center(i) for i in range(self.n_bin)])

class linear_grid(grid):

    def __init__(self, min, max, n_bin):
        self.min = float(min)
        self.max = float(max)
        self.n_bin = n_bin

    def scale(self, factor):
        self.min = self.min * factor
        self.max = self.max * factor

    def grid_size(self, index):
        return (self.max - self.min) / float(self.n_bin)

    def valid_bin(self, bin):
        if (bin >= 0) and (bin < self.n_bin):
            return True
        return False

    def find(self, values):
        indices = (floor((values - self.min) * self.n_bin
                         / (self.max - self.min))).astype(int)
        return indices

    def edge(self, index):
        if (index < 0) or (index > self.n_bin):
            raise Exception("index out of range: %d" % index)
        if index == self.n_bin:
            return self.max
        elif index == 0:
            return self.min
        else:
            return float(index) / float(self.n_bin) * (self.max - self.min) \
                   + self.min

    def center(self, index):
        if (index < 0) or (index >= self.n_bin):
            raise Exception("index out of range: %d" % index)
        return (float(index) + 0.5) / float(self.n_bin) \
               * (self.max - self.min) + self.min

    def half_sample(self):
        if self.n_bin % 2 != 0:
            raise Exception("n_bin must be an even number")
        return linear_grid(min = self.min, max = self.max,
                               n_bin = self.n_bin / 2)
        
class log_grid(grid):

    def __init__(self, min, max, n_bin):
        self.min = float(min)
        self.max = float(max)
        self.n_bin = n_bin

    def scale(self, factor):
        self.min = self.min * factor
        self.max = self.max * factor

    def grid_size(self, index):
        return (log10(self.max) - log10(self.min)) / float(self.n_bin)

    def valid_bin(self, bin):
        if (bin >= 0) and (bin < self.n_bin):
            return True
        return False

    def find(self, value):
        indices = (floor((log(value) - log(self.min)) * self.n_bin
                   / (log(self.max) - log(self.min)))).astype(int)
        return indices

    def edge(self, index):
        if (index < 0) or (index > self.n_bin):
            raise Exception("index out of range: %d" % index)
        if index == self.n_bin:
            return self.max
        elif index == 0:
            return self.min
        else:
            return math.exp(float(index) / float(self.n_bin)
                            * (log(self.max) - log(self.min))
                            + log(self.min))
        
    def center(self, index):
        if (index < 0) or (index >= self.n_bin):
            raise Exception("index out of range: %d" % index)
        return math.exp((float(index) + 0.5) / float(self.n_bin)
                        * (log(self.max) - log(self.min))
                        + log(self.min))

    def half_sample(self):
        if self.n_bin % 2 != 0:
            raise Exception("n_bin must be an even number")
        return log_grid(min = self.min, max = self.max,
                            n_bin = self.n_bin / 2)

def histogram_1d(x_values, x_axis, weights = None):
    """Make a 1D histogram.

    The histogram is of points at positions x_values[i] for each i.
    Example:
    >>> x_axis = partmc.log_grid(min = 1e-8, max = 1e-5, n_bin = 70)
    >>> hist = partmc.histogram_1d(diam, x_axis, weights = 1 / particles.comp_vols)
    >>> plt.semilogx(x_axis.centers(), hist)
    """
    if weights is not None:
        if len(x_values) != len(weights):
            raise Exception("x_values and weights have different lengths")
    x_bins = x_axis.find(x_values)
    hist = numpy.zeros([x_axis.n_bin])
    for i in range(len(x_values)):
        if x_axis.valid_bin(x_bins[i]):
            value = 1.0 / x_axis.grid_size(x_bins[i])
            if weights is not None:
                value *= weights[i]
            hist[x_bins[i]] += value
    return hist

def histogram_2d(x_values, y_values, x_axis, y_axis, weights = None, only_positive = True):
    """Make a 2D histogram.

    The histogram is of points at positions (x_values[i], y_values[i])
    for each i. Example:
    >>> x_axis = partmc.log_grid(min = 1e-8, max = 1e-5, n_bin = 70)
    >>> y_axis = partmc.linear_grid(min = 0, max = 1, n_bin = 50)
    >>> hist = partmc.histogram_2d(diam, bc_frac, x_axis, y_axis, weights = 1 / particles.comp_vols)
    >>> plt.pcolor(x_axis.edges(), y_axis.edges(), hist.transpose(),
                   norm = matplotlib.colors.LogNorm(), linewidths = 0.1)
    """
    if len(x_values) != len(y_values):
        raise Exception("x_values and y_values have different lengths")
    if weights is not None:
        if len(x_values) != len(weights):
            raise Exception("x_values and weights have different lengths")
    x_bins = x_axis.find(x_values)
    y_bins = y_axis.find(y_values)
    hist = numpy.zeros([x_axis.n_bin, y_axis.n_bin])
    for i in range(len(x_values)):
        if x_axis.valid_bin(x_bins[i]) and y_axis.valid_bin(y_bins[i]):
            value = 1.0 / (x_axis.grid_size(x_bins[i]) * y_axis.grid_size(y_bins[i]))
            if weights is not None:
                value *= weights[i]
            hist[x_bins[i], y_bins[i]] += value
    if only_positive:
        mask = numpy.ma.make_mask(hist <= 0.0)
        hist = numpy.ma.array(hist, mask = mask)
    return hist

def multival_2d(x_values, y_values, z_values, x_axis, y_axis, rand_arrange = True):
    """Make a 2D matrix with 0%/33%/66%/100% percentile values.

    The returned matrix represents z_values[i] at position
    (x_values[i], y_values[i]) for each i. Example:
    >>> x_axis = partmc.log_grid(min = 1e-8, max = 1e-5, n_bin = 140)
    >>> y_axis = partmc.linear_grid(min = 0, max = 1, n_bin = 100)
    >>> vals = partmc.multival_2d(diam, bc_frac, h2o, x_axis, y_axis)
    >>> plt.pcolor(x_axis.edges(), y_axis.edges(), vals.transpose(),
                   norm = matplotlib.colors.LogNorm(), linewidths = 0.1)
    """
    if len(x_values) != len(y_values):
        raise Exception("x_values and y_values have different lengths")
    if len(x_values) != len(z_values):
        raise Exception("x_values and z_values have different lengths")

    low_x_axis = x_axis.half_sample()
    low_y_axis = y_axis.half_sample()
    x_bins = low_x_axis.find(x_values)
    y_bins = low_y_axis.find(y_values)
    z = [[[] for j in range(low_y_axis.n_bin)]
         for i in range(low_x_axis.n_bin)]
    for i in range(len(x_values)):
        if low_x_axis.valid_bin(x_bins[i]) and low_y_axis.valid_bin(y_bins[i]):
            z[x_bins[i]][y_bins[i]].append(z_values[i])
    for x_bin in range(low_x_axis.n_bin):
        for y_bin in range(low_y_axis.n_bin):
            z[x_bin][y_bin].sort()
    grid = numpy.zeros([x_axis.n_bin, y_axis.n_bin])
    mask = numpy.zeros([x_axis.n_bin, y_axis.n_bin], bool)
    for x_bin in range(low_x_axis.n_bin):
        for y_bin in range(low_y_axis.n_bin):
            if len(z[x_bin][y_bin]) > 0:
                subs = [(0,0),(0,1),(1,0),(1,1)]
                if rand_arrange:
                    random.shuffle(subs)
                (sub_min, sub_max, sub_low, sub_high) = subs
                val_min = min(z[x_bin][y_bin])
                val_max = max(z[x_bin][y_bin])
                val_low = percentile(z[x_bin][y_bin], 0.3333)
                val_high = percentile(z[x_bin][y_bin], 0.6666)
                for (sub, val) in [(sub_min, val_min),
                                   (sub_max, val_max),
                                   (sub_low, val_low),
                                   (sub_high, val_high)]:
                    sub_i, sub_j = sub
                    i = x_bin * 2 + sub_i
                    j = y_bin * 2 + sub_j
                    grid[i,j] = val
                    mask[i,j] = True
    mask = logical_not(mask)
    vals = numpy.ma.array(grid, mask = mask)
    return vals

def pmc_histogram_2d(array, x_axis, y_axis, mask = None, inv_mask = None):
    data = []
    for i in range(x_axis.n_bin):
        for j in range(y_axis.n_bin):
            if ((mask == None) and (inv_mask == None) and (array[i,j] > 0.0)) \
                   or ((mask != None) and (mask[i,j] != 0)) \
                   or ((inv_mask != None) and (inv_mask[i,j] == 0)):
                data.append([x_axis.edge(i), x_axis.edge(i + 1),
                             y_axis.edge(j), y_axis.edge(j + 1),
                             array[i,j]])
    return data

def pmc_histogram_2d_multi(array_list, x_axis, y_axis, mask = None,
                           inv_mask = None):
    data = []
    for i in range(x_axis.n_bin):
        for j in range(y_axis.n_bin):
            if ((mask == None) and (inv_mask == None)
                and (array_list[0][i,j] > 0.0)) \
                or ((mask != None) and (mask[i,j] != 0)) \
                or ((inv_mask != None) and (inv_mask[i,j] == 0)):
                data_item = [x_axis.edge(i), x_axis.edge(i + 1),
                             y_axis.edge(j), y_axis.edge(j + 1)]
                for array in array_list:
                    data_item.append(array[i,j])
                data.append(data_item)
    return data

def time_of_day_string(time_seconds, separator = ":"):
    time_of_day = time_seconds % (24 * 3600.0)
    hours = int(time_of_day / 3600.0)
    minutes = int(time_of_day / 60.0) % 60
    seconds = int(time_of_day) % 60
    return "%02d%s%02d" % (hours, separator, minutes)

def read_history(constructor, directory, filename_pattern, print_progress = False):
    filenames = os.listdir(directory)
    filenames.sort()
    data = []
    filename_re = re.compile(filename_pattern)
    for filename in filenames:
        if filename_re.search(filename):
            if print_progress:
                print filename
            netcdf_filename = os.path.join(directory, filename)
            ncf = NetCDFFile(netcdf_filename)
            env_state = env_state_t(ncf)
            data.append([env_state.elapsed_time, constructor(ncf)])
            ncf.close()
    data.sort()
    return data

def read_any(constructor, directory, filename_pattern):
    filenames = os.listdir(directory)
    filename_re = re.compile(filename_pattern)
    for filename in filenames:
        if filename_re.search(filename):
            netcdf_filename = os.path.join(directory, filename)
            ncf = NetCDFFile(netcdf_filename)
            data = constructor(ncf)
            ncf.close()
            return data
    raise Exception("no NetCDF file found in %s matching %s"
                    % (directory, filename_pattern))

def get_filename_list(dir, file_pattern):
    filename_list = []
    filenames = os.listdir(dir)
    if len(filenames)  == 0:
        raise Exception("No files in %s match %s" % (dir, file_pattern))
    file_re = re.compile(file_pattern)
    for filename in filenames:
        match = file_re.search(filename)
        if match:
            output_key = match.group(1)
            full_filename = os.path.join(dir, filename)
            filename_list.append(full_filename)
    filename_list.sort()
    if len(filename_list) == 0:
        raise Exception("No files found in %s matching %s"
                        % (dir, file_pattern))
    return filename_list

def get_time_filename_list(dir, file_pattern):
    time_filename_list = []
    filenames = os.listdir(dir)
    if len(filenames) == 0:
        raise Exception("No files in %s match %s" % (dir, file_pattern))
    file_re = re.compile(file_pattern)
    for filename in filenames:
        match = file_re.search(filename)
        if match:
            output_key = match.group(1)
            netcdf_filename = os.path.join(dir, filename)
            ncf = NetCDFFile(netcdf_filename)
            env_state = env_state_t(ncf)
            time_filename_list.append([env_state.elapsed_time,
                                       netcdf_filename,
                                       output_key])
            ncf.close()
    time_filename_list.sort()
    if len(time_filename_list) == 0:
        raise Exception("No files found in %s matching %s"
                        % (dir, file_pattern))
    return time_filename_list

def find_nearest_time(time_indexed_data, search_time):
    min_diff = abs(search_time - time_indexed_data[0][0])
    min_i = 0
    for i in range(1,len(time_indexed_data)):
        diff = abs(search_time - time_indexed_data[i][0])
        if diff < min_diff:
            min_diff = diff
            min_i = i
    return min_i

def file_filename_at_time(time_filename_list, search_time):
    i = find_nearest_time(time_filename_list, search_time)
    return time_filename_list[i][1]

def load_constants(filename):
    constants = {}
    consts_file = open(filename) # filename of "constants.f90"
    in_const_t = False
    found_const_t = False
    start_re = re.compile("^ *type const_t *$")
    end_re = re.compile("^ *end type const_t *$")
    const_re = re.compile("^ *real[(]kind=dp[)] :: ([^ ]+) = ([-0-9.]+)d([-0-9]+) *$")
    for line in consts_file:
        if in_const_t:
            match = const_re.search(line)
            if match:
                name = match.group(1)
                mantissa = float(match.group(2))
                exponent = float(match.group(3))
                constants[name] = mantissa * 10.0**exponent
            if end_re.search(line):
                in_const_t = False
        else:
            if start_re.search(line):
                in_const_t = True
                found_const_t = True
    if not found_const_t:
        raise Exception("constants.f90 ended without finding const_t")
    if in_const_t:
        raise Exception("constants.f90 ended without finding end of const_t")
    return constants

def nan_to_value(x, value):
    values = x.copy()
    values.fill(value)
    return where(isnan(x), values, x)

def inf_to_value(x, value):
    values = x.copy()
    values.fill(value)
    return where(isinf(x), values, x)

def delta(arr):
    return (arr[1:] - arr[:-1])

def mid(arr):
    return (arr[1:] + arr[:-1]) / 2.0

def filter_inf(plot_data):
    return [[x, y] for [x, y] in plot_data if isfinite(y)]

#def sign(x):
#    if x > 0:
#        return 1
#    elif x < 0:
#        return -1
#    return 0

def mean(x):
    return float(sum(x)) / len(x)

def chop_sign_data_helper(plot_data, chopped_data):
    if len(plot_data) == 0:
        return chopped_data
    if sign(plot_data[0][1]) == 0:
        return chop_sign_data_helper(plot_data[1:], chopped_data)
    i = 0
    while (i < len(plot_data)) \
              and (sign(plot_data[i][1]) == sign(plot_data[0][1])):
        i += 1
    chopped_data.append(plot_data[0:i])
    return chop_sign_data_helper(plot_data[i:], chopped_data)

def chop_sign_data(plot_data):
    return chop_sign_data_helper(plot_data, [])

def chop_sign_data_iterative(plot_data):
    chopped_data = []
    start = 0
    while start < len(plot_data):
        while (start < len(plot_data)) and (sign(plot_data[start][1]) == 0):
            start += 1
        i = 0
        while (start + i < len(plot_data)) \
                and (sign(plot_data[start + i][1])
                     == sign(plot_data[start][1])):
            i += 1
        if i > 0:
            chopped_data.append(plot_data[start:(start + i)])
        start += i
    return chopped_data

def smooth(x,window_len=10,window='hanning'):
    """smooth the data using a window with requested size.
    
    This method is based on the convolution of a scaled window with the signal.
    The signal is prepared by introducing reflected copies of the signal 
    (with the window size) in both ends so that transient parts are minimized
    in the begining and end part of the output signal.
    
    input:
        x: the input signal 
        window_len: the dimension of the smoothing window
        window: the type of window from 'flat', 'hanning', 'hamming', 'bartlett', 'blackman'
            flat window will produce a moving average smoothing.

    output:
        the smoothed signal
        
    example:

    t=linspace(-2,2,0.1)
    x=sin(t)+randn(len(t))*0.1
    y=smooth(x)
    
    see also: 
    
    numpy.hanning, numpy.hamming, numpy.bartlett, numpy.blackman, numpy.convolve
    scipy.signal.lfilter
 
    TODO: the window parameter could be the window itself if an array instead of a string   
    """

    if x.ndim != 1:
        raise ValueError, "smooth only accepts 1 dimension arrays."

    if x.size < window_len:
        raise ValueError, "Input vector needs to be bigger than window size."

    if window_len<3:
        return x

    if not window in ['flat', 'hanning', 'hamming', 'bartlett', 'blackman']:
        raise ValueError, "Window is not one of 'flat', 'hanning', 'hamming', 'bartlett', 'blackman'"

    s=numpy.r_[2*x[0]-x[window_len:1:-1],x,2*x[-1]-x[-1:-window_len:-1]]
    if window == 'flat':
        w=ones(window_len,'d')
    else:
        w=eval('numpy.'+window+'(window_len)')

    y=numpy.convolve(w/w.sum(),s,mode='same')
    return y[window_len-1:-window_len+1]

def cumulative_plot_data(x, y_inc, start = 0.0, final = None):
    plot_data = []
    y = start
    for i in range(x.size):
        plot_data.append([x[i], y])
        y += y_inc[i]
        if (i == x.size - 1) and (final != None):
            y = final
        plot_data.append([x[i], y])
    return plot_data

def cumulative_hi_res(x, y_inc, start = 0.0, final = None,
                      min_x_step = None, min_y_step = None,
                      min_x_factor = None, min_y_factor = None):
    plot_data = []
    i = 0
    y = start
    plot_data.append([x[i], y])
    last_x, last_y = x[i], y
    for i in range(1,x.size):
        y += y_inc[i]
        if (i == x.size - 1) and (final != None):
            y = final
        if ((min_x_step != None) and (x[i] - last_x > min_x_step)) \
                or ((min_y_step != None) and (y - last_y > min_y_step)) \
                or ((min_x_factor != None) and (x[i]/last_x > min_x_factor)) \
                or ((min_y_factor != None) and (y/last_y > min_y_factor)) \
                or (i == x.size - 1):
            plot_data.append([x[i], y])
            last_x = x[i]
            last_y = y
    return plot_data

def percentile(data, p):
    # p in [0,1]
    # data must be sorted
    i = int(floor((len(data) - 1) * p + 0.5))
    return data[i]

def grid_plot_data(x, y, z, x_axis, y_axis, hi_x_axis, hi_y_axis):
    x_bin = x_axis.find(x)
    y_bin = y_axis.find(y)
    values = [[[] for j in range(y_axis.n_bin)]
              for i in range(x_axis.n_bin)]
    for i in range(len(x)):
        if x_axis.valid_bin(x_bin[i]) and y_axis.valid_bin(y_bin[i]):
            values[x_bin[i]][y_bin[i]].append(z[i])
    for x_bin in range(x_axis.n_bin):
        for y_bin in range(y_axis.n_bin):
            values[x_bin][y_bin].sort()
    grid = numpy.zeros([hi_x_axis.n_bin, hi_y_axis.n_bin])
    mask = numpy.zeros([hi_x_axis.n_bin, hi_y_axis.n_bin], int)
    for x_bin in range(x_axis.n_bin):
        for y_bin in range(y_axis.n_bin):
            if len(values[x_bin][y_bin]) > 0:
                subs = [(0,0),(0,1),(1,0),(1,1)]
                random.shuffle(subs)
                (sub_min, sub_max, sub_low, sub_high) = subs
                val_min = min(values[x_bin][y_bin])
                val_max = max(values[x_bin][y_bin])
                val_low = percentile(values[x_bin][y_bin], 0.3333)
                val_high = percentile(values[x_bin][y_bin], 0.6666)
                for (sub, val) in [(sub_min, val_min),
                                   (sub_max, val_max),
                                   (sub_low, val_low),
                                   (sub_high, val_high)]:
                    sub_i, sub_j = sub
                    i = x_bin * 2 + sub_i
                    j = y_bin * 2 + sub_j
                    grid[i,j] = val
                    mask[i,j] = 1
    plot_data = pmc_histogram_2d(grid, hi_x_axis, hi_y_axis, mask = mask)
    return plot_data
