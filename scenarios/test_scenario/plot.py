from netCDF4 import Dataset
import math
import numpy as np

def compDistortion1D(exact, observed, exactWeights=None, observedWeights=None,
		useLog=False, useRel=None):
		
	if exactWeights == None:
		exactWeights = np.ones(len(exact))
	if observedWeights == None:
		observedWeights = np.ones(len(observed))
		
	if len(exact) != len(exactWeights):
		raise Exception('exact and exactWeights must have same length')
	if len(observed) != len(observedWeights):
		raise Exception('observed and observedWeights must have same length')
	if exactWeights.min() < 0.0:
		raise Exception('exactWeights must be a non-negative array')
	if observedWeights.min() < 0.0:
		raise Exception('observedWeights must be a non-negative array')
		
	if useRel == None:
		useRel = not useLog
	if useLog:
		exact = np.log(exact)
		observed = np.log(observed)
		
	netWeightExact = float(exactWeights.sum())
	netWeightObserved = float(observedWeights.sum())
	
	if netWeightExact <= 0.0 and netWeightObserved <= 0.0:
		return (0.0, 0.0)
	elif netWeightExact <= 0.0 or netWeightObserved <= 0.0:
		return (float('inf'), float('inf'))
		
	netWeightError = abs(netWeightObserved-netWeightExact)/netWeightExact
	
	permute = exact.argsort()
	exact = exact[permute]
	exactWeights = exactWeights[permute]
	
	permute = observed.argsort()
	observed = observed[permute]
	observedWeights = observedWeights[permute]
	
	exactWeights /= netWeightExact
	observedWeights /= observedWeights
	
	netNormExact = 0.0
	if useRel:
		for i in xrange(0, len(exact)):
			netNormExact += exactWeights[i]*float(exact[i])**2
		if netNormExact <= 0.0:
			return (float('inf'), netWeightError)
	else:
		netNormExact = 1.0
	
	arrA = exact
	arrB = observed
	weightsA = exactWeights
	weightsB = observedWeights
	
	indA = 0
	indB = 0
	weightA = 0.0
	weightB = weightsB[0]
	
	sum = 0.0
	while indA < len(arrA):
		weightA = weightsA[indA]
		if weightB < weightA:
			arrA, arrB = arrB, arrA
			weightsA, weightsB = weightsB, weightsA
			indA, indB = indB, indA
			weightA, weightB = weightB, weightA
		sum += weightA*float(arrA[indA]-arrB[indB])**2
		weightB -= weightA
		indA += 1
		
	return (math.sqrt(sum/netNormExact), netWeightError)

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
