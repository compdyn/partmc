import os
import perturb_particles

charsPerOutFileNum = 2
partmcCommand = '../../build/partmc'

def numToStr(num, numChars):
	strNum = str(num)
	if len(strNum) > numChars:
		raise Exception(strNum + ' is more than ' + str(numChars) + ' chars')
	return '0'*(numChars - len(strNum)) + strNum

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
	
def ensureOutDirectory(runNumber):
	dirName = getOutDir(runNumber)
	if os.path.exists(dirName):
		if not os.path.isdir(dirName):
			raise Exception(dirName + ' is not a directory')
	else:
		os.makedirs(dirName)
		
def getOutDir(runNumber):
	return 'out/run' + numToStr(runNumber, charsPerOutFileNum)
	
def findRestartIndexStr():
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
	return text[startI:endI]

def runSimulations(numSimulations):
	if numSimulations == 0:
		return
	numToStr(numSimulations-1, charsPerOutFileNum)
	lastRunNum = 0
	restartIndexStr = findRestartIndexStr()
	for runNum in xrange(0, numSimulations):
		outDir = getOutDir(runNum)
		updateSpecFilesRunNumbers(lastRunNum, runNum)
		ensureOutDirectory(runNum)
		os.system(partmcCommand + ' urban_plume.spec')
		os.system('cp ' + outDir + '/normal_0001_' + restartIndexStr
			+ '.nc ' + outDir + '/perturbed_0001_' + restartIndexStr
			+ '.nc')
		perturbResult = perturb_particles.perturbFile(
			outDir + '/perturbed_0001_' + restartIndexStr + '.nc')
		logFile = open(outDir + '/perturb.log', 'w')
		logFile.write(str(perturbResult))
		logFile.close()
		os.system(partmcCommand + ' urban_plume_restart.spec')
	updateSpecFilesRunNumbers(numSimulations-1, 0)
		