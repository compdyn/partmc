from netCDF4 import Dataset
import math

def readArray(dataSet, var):
	madeDataSet = isinstance(dataSet, str)
	if madeDataSet:
		dataSet = Dataset(dataSet, 'r', format='NETCDF3_CLASSIC')
	result = dataSet.variables[var][:]
	if madeDataSet:
		dataSet.close()
	return result

def readMassVectors(dataSet):
	return readArray(dataSet, 'aero_particle_mass').T
	
def readDensities(dataSet):
	return readArray(dataSet, 'aero_density')

BC_index = 19 - 1
H2O_index = 20 - 1
numSpecies = 20
nonWetId = np.ones(numSpecies)
nonWetId[H2O_index] = 0
densities = readDensities('out/run0001/normal_0001_00000001.nc')
nonWetInvDensities = nonWetId/densities
	
def sphereVol2Rad(vol):
	return ((.75/math.pi)*vol)**(1.0/3.0)
	
def sphereVol2Diam(vol):
	return 2*sphereVol2Rad(vol)
	
	
def compDryVol(mass_vector):
	return (mass_vector*nonWetInvDensities).sum()
	
def compDryDiam(mass_vector):
	return sphereVol2Diam(compDryVol(mass_vector))
	
def compMass(mass_vector):
	return mass_vector.sum()
	
def compDryMass(mass_vector):
	return (mass_vector*nonWetId).sum()

def compBcDryMassFraction(mass_vector):
	return mass_vector[BC_index]/compDryMass(mass_vector)
	
def compWaterMassFraction(mass_vector):
	return mass_vector[H2O_index]/compMass(mass_vector)
	
def return1(mass_vector):
	return 1.0
