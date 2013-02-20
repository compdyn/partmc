from netCDF4 import Dataset
import numpy as np
import numpy.linalg as la

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
	writeMassVectors(data, vectors)
	data.close()
	return result
	
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
			
	distortion = 0.0
	relDistortion = 0.0
	relDistortion2 = 0.0
	logDistortion = 0.0
	
	for vector in vectors:
		key = getCellKey(vector, cellGrowthFactor)
		newVector = map[key][0]
		
		distortion += la.norm(vector - newVector)**2
		relDistortion += (la.norm(vector - newVector)/la.norm(vector))**2
		relDistortion2 += la.norm((vector - newVector)/(vector+1e-200))**2
		logDistortion += computeLogDistortion(vector, newVector)
		
		vector[:] = newVector
		
	distortions = [distortion, relDistortion, relDistortion2, logDistortion]
	for i in xrange(0, len(distortions)):
		distortions[i] /= vectors.shape[0]*vectors.shape[1]
		
	superWeightMap = {}
	for entry in map.values():
		if superWeightMap.has_key(entry[1]):
			superWeightMap[entry[1]] += 1
		else:
			superWeightMap[entry[1]] = 1
	superWeights = superWeightMap.items()
	superWeights.sort()
	return (distortions, superWeights)
	
def getCellKey(vector, cellGrowthFactor):
	doubleArr = np.floor(np.log(vector)/np.log(cellGrowthFactor))
	return tuple(doubleArr.astype(np.int32))
