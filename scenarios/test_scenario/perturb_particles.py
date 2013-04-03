from netCDF4 import Dataset
import numpy as np
import numpy.linalg as la
from scipy.cluster import vq

def readMassVectors(dataSet):
	var = dataSet.variables['aero_particle_mass']
	return var[:,:].T

def writeMassVectors(dataSet, newVectors):
	newVectors = newVectors.T
	var = dataSet.variables['aero_particle_mass']
	assert newVectors.shape == var.shape
	var[:,:] = newVectors

def perturbFile(fileName):
	data = Dataset(fileName, 'r+', format='NETCDF3_CLASSIC')
	vectors = readMassVectors(data)
	result = mergeVectorsInGrid(vectors, 2)
	writeMassVectors(data, result[0])
	data.close()
	distortions = computeDistortions(vectors, result[0])
	return (distortions, result[1])
	
def readMassVectorsFromFile(filename):
	data = Dataset(filename, 'r', format='NETCDF3_CLASSIC')
	vectors = readMassVectors(data)
	data.close()
	return vectors
	
# def printDiff(vectorA, vectorB):
# 	line = '['
# 	for i in xrange(0, len(vectorA)):
# 		if(vectorA[i] != vectorB[i]):
# 			line += str(i) + ': '
# 			#line += str(abs(vectorB[i] - vectorA[i])/abs(vectorA[i]))
# 			line += str(vectorA[i]) + ' -> ' + str(vectorB[i])
# 			line += '; '
# 	line += ']'
# 	print line

def computeLogDistortion(vectorA, vectorB):
	sum = 0.0
	for i in xrange(0, len(vectorA)):
		if vectorA[i] == 0.0:
			if vectorB[i] != 0.0:
				sum = float('inf')
		else:
			sum += (np.log(vectorA[i]) - np.log(vectorB[i]))**2
	return sum

def mergeVectorsInGrid(vectors, cellGrowthFactor):
	map = {}
	for vector in vectors:
		key = getCellKey(vector, cellGrowthFactor)
		if map.has_key(key):
			entry = map[key]
			entry[1] += 1
			alpha = 1.0/entry[1]
			entry[0] = (1.0-alpha)*entry[0] + alpha*vector
		else:
			map[key] = [np.array(vector), 1]
	
	newVectors = np.zeros(vectors.shape)
	
	for i in xrange(0, len(vectors)):
		key = getCellKey(vectors[i], cellGrowthFactor)
		newVectors[i,:] = map[key][0]
		
	superWeightMap = {}
	for entry in map.values():
		if superWeightMap.has_key(entry[1]):
			superWeightMap[entry[1]] += 1
		else:
			superWeightMap[entry[1]] = 1
	superWeights = superWeightMap.items()
	superWeights.sort()
	
	return (newVectors, superWeights)
	
def kmeansVQ(vectors, numCodevectors):
	codevectors = vq.kmeans(vectors, numCodevectors)
	codevectorWeights = np.zeros(len(codevectors))
	newVectors = np.zeros(vectors.shape)
	for i in xrange(0, len(vectors)):
		vector = vectors[i]
		codeVectorIndex = -1
		minDist = float("infinity")
		for j in xrange(0, len(codevectors)):
			dist = la.norm(codevectors[j] - vector)
			if dist < minDist:
				minDist = dist
				codeVectorIndex = j
		newVectors[i] = codevectors[codeVectorIndex]
		codevectorWeights[codeVectorIndex] += 1
		
	superWeightMap = {}
	for weight in codevectorWeights:
		if superWeightMap.has_key(weight):
			superWeightMap[weight] += 1
		else:
			superWeightMap[weight] = 1
	
	return (newVectors, superWeightMap.items())
	
def computeDistortions(oldVectors, newVectors):
	distortion = 0.0
	relDistortion = 0.0
	relDistortion2 = 0.0
	logDistortion = 0.0
	
	for i in xrange(0, len(oldVectors)):
		sqDist = la.norm(oldVectors[i] - newVectors[i])**2
		distortion += sqDist
		relDistortion += sqDist/la.norm(oldVectors[i])**2
		relDistortion2 += la.norm((oldVectors[i] - newVectors[i])
			/ (oldVectors[i]+1e-200))**2
		logDistortion += computeLogDistortion(oldVectors[i], newVectors[i])
		
	distortions = [distortion, relDistortion, relDistortion2, logDistortion]
	for i in xrange(0, len(distortions)):
		distortions[i] /= oldVectors.shape[0]*oldVectors.shape[1]
		
	return distortions
	
def getCellKey(vector, cellGrowthFactor):
	doubleArr = np.floor(np.log(vector)/np.log(cellGrowthFactor))
	return tuple(doubleArr.astype(np.int32))
