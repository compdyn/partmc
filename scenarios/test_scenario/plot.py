from netCDF4 import Dataset
import math
import numpy as np
import matplotlib.pylab as pl
import numpy.random as random

charsPerOutFileNum = 4

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
		
	origExactLen = len(exact)
	origObservedLen = len(observed)
		
	if useRel == None:
		useRel = not useLog
	if useLog:
		nz = exact.nonzero()
		exact = np.log(exact[nz])
		exactWeights = exactWeights[nz]
		nz = observed.nonzero()
		observed = np.log(observed[nz])
		observedWeights = observedWeights[nz]
		
	netWeightExact = float(exactWeights.sum())
	netWeightObserved = float(observedWeights.sum())
	
	print 'exact:', len(exact), 'particles used out of', origExactLen
	print 'observed:', len(observed), 'particles used out of', origObservedLen
	
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
	observedWeights /= netWeightObserved
	
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

def numToStr(num, numChars):
	strNum = str(num)
	if len(strNum) > numChars:
		raise Exception(strNum + ' is more than ' + str(numChars) + ' chars')
	return strNum.rjust(numChars, '0')
		
def getOutDir(runNumber):
	return 'out/run' + numToStr(runNumber, charsPerOutFileNum)

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


def compParticleFeatures(numRuns, fileSuffix, featureFunctions):
	features = np.zeros((len(featureFunctions), 0))
	for run in xrange(1, numRuns+1):
		outDir = getOutDir(run)
		massVectors = readMassVectors(outDir + '/' + fileSuffix)
		values = np.zeros((len(featureFunctions), len(massVectors)))
		for j in xrange(0, len(massVectors)):
			for i in xrange(0, len(features)):
				values[i,j] = featureFunctions[i](massVectors[j])
		features = np.hstack((features, values))
	return features
	
def compPartMCDistortion(numRuns, fileSuffixNormal, fileSuffixPerturbed,
		xFeature, yFeature, useLog=False, useRel=None):
	print 'computing difference between', fileSuffixNormal, 'and', \
		fileSuffixPerturbed
	featuresA = compParticleFeatures(numRuns, fileSuffixNormal,
		(xFeature, yFeature))
	print fileSuffixNormal, 'contains', len(featuresA[0]), 'particles'
	featuresB = compParticleFeatures(numRuns, fileSuffixPerturbed,
		(xFeature, yFeature))
	print fileSuffixPerturbed, 'contains', len(featuresB[0]), 'particles'
	return compDistortion1D(featuresA[0,:], featuresB[0,:],
		featuresA[1,:]/numRuns, featuresB[1,:]/numRuns,
		useLog=useLog, useRel=useRel)
	
def plotTest1(xFeature, yFeature, useLog=False, useRel=None):
	numDataPoints = 25
	splitIndex = 1
	numRuns = 4
	distError = np.zeros(numDataPoints)
	magError = np.zeros(numDataPoints)
	for i in xrange(0, numDataPoints):
		result = compPartMCDistortion(numRuns,
			'normal_0001_' + numToStr(splitIndex + i, 8) + '.nc',
			'perturbed_0001_' + numToStr(1 + i, 8) + '.nc',
			xFeature, yFeature, useLog=useLog, useRel = useRel)
		distError[i] = result[0]
		magError[i] = result[1]
	hours = range(splitIndex, splitIndex + numDataPoints)
	pl.hold(True)
	pl.semilogy(hours, distError, label='distribution error')
	pl.semilogy(hours, magError, label='magnitude error')
	return (hours, distError, magError)
	
def plotNumParticlesDivergence():
	numDataPoints = 18
	splitIndex = 8
	numParticlesA = np.zeros(numDataPoints)
	numParticlesB = np.zeros(numDataPoints)
	for i in xrange(0, numDataPoints):
		numParticlesA[i] = len(readMassVectors(getOutDir(1) + '/' + 
			'normal_0001_' + numToStr(splitIndex + i, 8) + '.nc'))
		numParticlesB[i] = len(readMassVectors(getOutDir(1) + '/' + 
			'perturbed_0001_' + numToStr(1 + i, 8) + '.nc'))
	hours = range(splitIndex, splitIndex + numDataPoints)
	pl.hold(True)
	pl.semilogy(hours, numParticlesA, label='normal run')
	pl.semilogy(hours, numParticlesB, label='perturbed run')
	return (hours, numParticlesA, numParticlesB)
	
def gaussTest():
	xVals = np.logspace(0, 4)
	distError = np.zeros(len(xVals))
	magError = np.zeros(len(xVals))
	numTrials = 10
	for t in xrange(0, numTrials):
		for i in xrange(0, len(xVals)):
			valsA = random.normal(size = max(1,
				xVals[i]*(1 + random.normal()/math.sqrt(xVals[i]))))
			valsB = random.normal(size = max(1,
				xVals[i]*(1 + random.normal()/math.sqrt(xVals[i]))))
			result = compDistortion1D(valsA, valsB, useRel=False)
			distError[i] += result[0]
			magError[i] += result[1]
	pl.hold(True)
	pl.loglog(xVals, distError/numTrials, label='distribution error')
	pl.loglog(xVals, magError/numTrials, label='magnitude error')
	pl.xlabel('expected num samples')
	pl.ylabel('expected error')
	pl.legend()
	pl.title('1D distribution distance measure test')