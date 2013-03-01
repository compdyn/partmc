class FileData:
	__dataArr = []
	
	def readFile(self, fileName):
		self.__dataArr = []
		f = open(fileName, 'r')
		lines = f.readlines(-1)
		f.close()
		for line in lines:
			commentI = line.find('#')
			if commentI >= 0:
				line = line[0:commentI]
			words = line.split()
			if len(words) > 0:
				self.addEntry(words)
				
	def copy(self, other):
		self.__dataArr = []
		for elem in other.__dataArr:
			self.__dataArr.append(elem[:])
	
	def addEntry(self, keyAndValues, index=-1):
		if index < 0:
			index = len(self.__dataArr)
		if len(keyAndValues) < 1:
			raise Exception('length of keyAndValues must be at least 1')
		keyAndValues[0] = str(keyAndValues[0])
		for i in xrange(1, len(keyAndValues)):
			try:
				keyAndValues[i] = float(keyAndValues[i])
				intVal = int(keyAndValues[i])
				if intVal == keyAndValues[i]:
					keyAndValues[i] = intVal
			except ValueError:
				keyAndValues[i] = str(keyAndValues[i])
		self.__dataArr.insert(index, keyAndValues)
		
	def getIndices(self, key):
		result = []
		for i in xrange(0, len(self.__dataArr)):
			if self.__dataArr[i][0] == key:
				result.append(i)
		return result
		
	def getIndex(self, key):
		result = self.getIndices(key)
		if len(result) == 0:
			return -1
		if len(result) > 1:
			raise Exception('Key ' + key + ' has ambiguous values')
		return result[0]
		
	def getEntries(self, key):
		indices = self.getIndices(key)
		result = []
		for i in indices:
			result.append(self.__dataArr[i])
		return result
		
	def getEntry(self, key):
		index = self.getIndex(key)
		if index < 0:
			return None
		return self.__dataArr[index]
		
	def writeToFile(self, fileName):
		f = open(fileName, 'w')
 		for words in self.__dataArr:
 			for word in words:
 				f.write(str(word) + ' ')
 			f.write('\n')
 		f.close()
 		
 	def size(self):
 		return len(self.__dataArr)
 	
 	def getEntryByIndex(self, index):
 		return self.__dataArr[index]
 		
 	def removeEntryByIndex(self, index):
 		del self.__dataArr[index]
 	
 	def removeEntries(self, key):
 		indices = self.getIndices(key)
 		indices.reverse()
 		for i in indices:
 			self.removeEntryByIndex(i)

	def printAll(self):
		for words in self.__dataArr:
 			print words