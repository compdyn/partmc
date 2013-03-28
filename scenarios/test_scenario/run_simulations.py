import os
import shutil
import perturb_particles as pp
from FileData import FileData
import numpy as np

charsPerOutFileNum = 4

def numToStr(num, numChars):
	strNum = str(num)
	if len(strNum) > numChars:
		raise Exception(strNum + ' is more than ' + str(numChars) + ' chars')
	return strNum.rjust(numChars, '0')

def updateSpecFileRunNumbers(fileName, oldNum, newNum, numChanges):
	f = open(fileName, 'r+')
	text = f.read(-1)
	index = -1
	searchStr = getOutDir(oldNum)
	newNumStr = numToStr(newNum, charsPerOutFileNum)
	actNumChanges = 0
	while True:
		index = text.find(searchStr+'/', index+1)
		if index < 0:
			break
		f.seek(index + len(searchStr) - charsPerOutFileNum)
		f.write(newNumStr)
		actNumChanges += 1
	f.close()
	if numChanges != actNumChanges:
		raise Exception('Expected ' + str(numChanges)
			+ ' file writes, but performed ' + str(actNumChanges))
	
#def updateSpecFilesRunNumbers(oldNum, newNum):
#	updateSpecFileRunNumbers('urban_plume_altered.spec', oldNum, newNum, 1)
#	updateSpecFileRunNumbers('urban_plume_restart.spec', oldNum, newNum, 2)
	
def ensureDirectory(dirName):
	if os.path.exists(dirName):
		if not os.path.isdir(dirName):
			raise Exception(dirName + ' is not a directory')
	else:
		os.makedirs(dirName)
		
def getOutDir(runNumber):
	return 'out/run' + numToStr(runNumber, charsPerOutFileNum)
	
def runPartMC(specFile):
	if specFile.find('0') >= 0 or not specFile.endswith('.spec'):
		raise Exception('invalid spec file: ' + specFile)
	code = os.system('../../build/partmc ' + specFile)
	if code != 0:
		raise Exception('bad PartMC exit code: ' + str(code))
		
def trimExtension(file, extension):
	checkExtension(extension)
	return file[0:-len(extension)]

def checkExtension(file, extension):
	if not file.endswith(extension):
		raise Exception('Expected ' + extension + ' file, got ' + file)
		
def makeRestartSpecFile(baseSpecFileName, restartIndex):
	checkExtension(baseSpecFileName, '.spec')
	data = FileData()
	data.readFile(baseSpecFileName)
	outPrefix = data.getEntry('output_prefix')
	endOfOutPrefix = '/normal'
	if not outPrefix[1].endswith(endOfOutPrefix):
		raise Exception('out_prefix should end with ' + endOfOutPrefix)
	data.getEntry('restart')[1] = 'yes'
	data.addEntry(['restart_file',
		outPrefix[1][0:-len(endOfOutPrefix)] + '/perturbed_start.nc'],
		data.getIndex('restart') + 1 )
	outPrefix[1] = outPrefix[1][0:-len(endOfOutPrefix)] + '/perturbed'
	data.getEntry('t_max')[1] -= indexToTime(data, restartIndex)
	data.removeEntries('gas_data')
	data.removeEntries('gas_init')
	data.removeEntries('aerosol_data')
	data.removeEntries('aerosol_init')
	data.writeToFile('restart_' + baseSpecFileName)
	
def adjustDatFilesForRestart(fileList, restartTime):
	for file in fileList:
		checkExtension(file, '.dat')
		data = FileData()
		data.readFile(file)
		entry = data.getEntry('time')
		for i in xrange(1, len(entry)):
			entry[i] -= restartTime
		data.writeToFile(file)
		
datFileBackupPrefix = 'original_version_of_'

def copyDatFiles(fileList):
	for file in fileList:
		checkExtension(file, '.dat')
		shutil.copyfile(file, datFileBackupPrefix + file)
	
def restoreDatFiles(fileList):
	for file in fileList:
		checkExtension(file, '.dat')
		shutil.copyfile(datFileBackupPrefix + file, file)
		os.remove(datFileBackupPrefix + file)
		
def indexToTime(file, timeIndex):
	if isinstance(file, str):
		f = FileData()
		f.readFile(file)
		file = f
	result = file.getEntry('t_output')[1]*(timeIndex-1)
	if result > file.getEntry('t_max')[1]:
		raise Exception('invalid time index ' + timeIndex)
	return result
	
def getNumOutputFiles(fileName):
	f = FileData()
	f.readFile(fileName)
	interval = f.getEntry('t_output')[1]
	full = f.getEntry('t_max')[1]
	if full % interval != 0:
		raise Exception('t_output does not divide t_max')
	return (full / interval) + 1
	
def checkNumSimulations(numSimulations):
	if numSimulations <= 0:
		raise Exception('numSimulations must be positive')
	numToStr(numSimulations, charsPerOutFileNum)
	
def runSimulationsNormal(numSimulations):
	checkNumSimulations(numSimulations)
	for runNum in xrange(1, numSimulations+1):
		updateSpecFileRunNumbers('urban_plume.spec', runNum-1, runNum, 1)
		ensureDirectory(getOutDir(runNum))
		runPartMC('urban_plume.spec')
	updateSpecFileRunNumbers('urban_plume.spec', numSimulations, 0, 1)
	
def perturbSimulations(numSimulations, restartIndex):
	checkNumSimulations(numSimulations)
	for runNum in xrange(1, numSimulations+1):
		outDir = getOutDir(runNum)
		shutil.copyfile(
			outDir + '/normal_0001_' + numToStr(restartIndex, 8) + '.nc',
			outDir + '/perturbed_start.nc')
		perturbResult = pp.perturbFile(outDir + '/perturbed_start.nc')
		logFile = open(outDir + '/perturb.log', 'w')
		logFile.write('Perturbed output file ' + str(restartIndex) + '\n')
		logFile.write(str(perturbResult))
		logFile.close()
		
def getAveDistortions(numSimulations, restartIndex, methodName):
	numDataPoints = 7
	sums = np.zeros((numDataPoints, 4))
	numCodevectorSums = np.zeros(numDataPoints)
	params = getDistortParameters(methodName, numDataPoints)
	for runNum in xrange(1, numSimulations+1):
		print 'runNum: ' + str(runNum)
		vectors = pp.readMassVectorsFromFile(getOutDir(runNum)
			+ '/normal_0001_' + numToStr(restartIndex, 8) + '.nc')
		for i in xrange(0, numDataPoints):
			print 'Data point: ' + str(i)
			result = distortVectors(vectors, methodName, params[i])
			numCodevectorSums[i] += (np.array(result[1]).T)[1].sum()
			sums[i,:] += pp.computeDistortions(vectors, result[0])
	sums /= numSimulations
	numCodevectorSums /= numSimulations
	return (numCodevectorSums, sums)
				
def distortVectors(vectors, methodName, param):
	if methodName == 'grid':
		return pp.mergeVectorsInGrid(vectors, param)
	elif methodName == 'kmeans':
		return pp.kmeansVQ(vectors, param)
	else:
		raise Exception('invalid methodName: ' + str(methodName))
		
def getDistortParameters(methodName, numPoints):
	if methodName == 'grid':
		return np.logspace(.1, 1, numPoints)
	elif methodName == 'kmeans':
		return np.logspace(0, 4, numPoints)
	else:
		raise Exception('invalid methodName: ' + str(methodName))
	
def runSimulationsRestart(numSimulations, restartIndex):
	checkNumSimulations(numSimulations)
	datFilesWithTime = ['aero_back.dat', 'aero_emit.dat', 'gas_back.dat',
		'gas_emit.dat', 'height.dat', 'pres.dat', 'temp.dat']
	makeRestartSpecFile('urban_plume.spec', restartIndex)
	restartTime = indexToTime('urban_plume.spec', restartIndex)
	copyDatFiles(datFilesWithTime)
	adjustDatFilesForRestart(datFilesWithTime, restartTime)
	for runNum in xrange(1, numSimulations+1):
		updateSpecFileRunNumbers('restart_urban_plume.spec', runNum-1, runNum,
			2)
		runPartMC('restart_urban_plume.spec')
	restoreDatFiles(datFilesWithTime)
	updateSpecFileRunNumbers('restart_urban_plume.spec', numSimulations, 0, 2)

def runSimulations(numSimulations, restartIndex):
	runSimulationsNormal(numSimulations)
	perturbSimulations(numSimulations, restartIndex)
	runSimulationsRestart(numSimulations, restartIndex)

def perturbedRunToOutDir(numSimulations, restartIndex):
	os.system('rm out/*.nc')
	numOutputFiles = getNumOutputFiles('urban_plume.spec')
	for runNum in xrange(1, numSimulations+1):
		outDir = getOutDir(runNum)
		for outIndex in xrange(1, numOutputFiles+1):
			prefix = 'normal'
			srcOutIndex = outIndex
			if outIndex >= restartIndex:
				prefix = 'perturbed'
				srcOutIndex = outIndex - restartIndex + 1
			src = outDir + '/' + prefix + '_0001_'
			src += numToStr(srcOutIndex,8) + '.nc'
			dst = 'out/urban_plume_' + numToStr(runNum, 4) + '_'
			dst += numToStr(outIndex, 8) + '.nc'
			shutil.copyfile(src, dst)
		
def normalRunToOutDir(numSimulations):
	perturbedRunToOutDir(numSimulations, 1000000000)
	