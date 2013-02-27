import os
import shutil
import perturb_particles

charsPerOutFileNum = 2

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
		index = text.find(searchStr, index+1)
		if index < 0:
			break
		f.seek(index + len(searchStr) - charsPerOutFileNum)
		f.write(newNumStr)
		actNumChanges += 1
	f.close()
	if numChanges != actNumChanges:
		raise Exception('Expected ' + str(numChanges)
			+ ' file writes, but performed ' + str(actNumChanges))
	
def updateSpecFilesRunNumbers(oldNum, newNum):
	updateSpecFileRunNumbers('urban_plume.spec', oldNum, newNum, 1)
	updateSpecFileRunNumbers('urban_plume_restart.spec', oldNum, newNum, 2)
	
def ensureDirectory(dirName):
	if os.path.exists(dirName):
		if not os.path.isdir(dirName):
			raise Exception(dirName + ' is not a directory')
	else:
		os.makedirs(dirName)
		
def getOutDir(runNumber):
	return 'out/run' + numToStr(runNumber, charsPerOutFileNum)
	
def specFileToDict(specFile):
	f = open('urban_plume_restart.spec')
	lines = f.readlines(-1)
	f.close()
	dict = {}
	for line in lines:
		commentI = line.find('#')
		if commentI >= 0:
			line = line[0:commentI]
		words = line.split()
		if len(words) == 0:
			continue
		if len(words) != 2 or dict.has_key(words[0]):
			dict[words[0]] = 'UNKNOWN VALUE'
		else:
			dict[words[0]] = words[1]
	return dict
	
def findRestartIndex():
	f = open('urban_plume_restart.spec')
	text = f.read(-1)
	f.close()
	searchStr = 'restart_file out/run00/perturbed_0001_'
	i = text.find(searchStr)
	if i < 0:
		raise Exception('cannot find appropriate restart file')
	startI = i + len(searchStr)
	endI = startI + 8
	if endI > len(text):
		raise Exception('cannot find appropriate restart file')
	indexStr = text[startI:endI]
	index = int(indexStr)
	if index < 0:
		raise Exception('cannot find appropriate restart file')
	return index
	
def runPartMC(specFile):
	if specFile.find('0') >= 0 or not specFile.endswith('.spec'):
		raise Exception('invalid spec file: ' + specFile)
	code = os.system('../../build/partmc ' + specFile)
	if code != 0:
		raise Exception('bad PartMC exit code: ' + str(code))

def runSimulations(numSimulations):
	if numSimulations <= 0:
		return
	numToStr(numSimulations-1, charsPerOutFileNum)
	lastRunNum = 0
	restartIndex = findRestartIndex()
	for runNum in xrange(0, numSimulations):
		outDir = getOutDir(runNum)
		updateSpecFilesRunNumbers(lastRunNum, runNum)
		ensureDirectory(outDir)
		runPartMC('urban_plume.spec')
		for i in xrange(1, restartIndex+1):
			shutil.copyfile(
				outDir + '/normal_0001_' + numToStr(i, 8) + '.nc',
				outDir + '/perturbed_0001_' + numToStr(i, 8) + '.nc')
		perturbResult = perturb_particles.perturbFile(
			outDir + '/perturbed_0001_' + numToStr(restartIndex, 8) + '.nc')
		logFile = open(outDir + '/perturb.log', 'w')
		logFile.write('Perturbed output file ' + str(restartIndex) + '\n')
		logFile.write(str(perturbResult))
		logFile.close()
		runPartMC('urban_plume_restart.spec')
	updateSpecFilesRunNumbers(numSimulations-1, 0)
