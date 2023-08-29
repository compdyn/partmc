#!/usr/bin/env python

import numpy as np

class TSPlots(object):
    def __init__(self, data_setting = None, casesList = None, casesDict = None):

        assert not(data_setting) is None
        self.data_setting = data_setting
        assert not( (casesList is None) and (casesDict is None) )
        if casesDict is None:
            self.casesDict = {}
            for obj_case in casesList:
                self.casesDict[obj_case.caseName] = obj_case
        else:
            self.casesDict = casesDict
        
    def plot(self, ax, legend = False, grid = True):
        for dataName in self.data_setting:
            dataName_split = dataName.strip().split(":")
            if len(dataName_split) == 2:
                caseName, varName = dataName_split
                assert caseName in self.casesDict
                ensemble_name = self.casesDict[caseName].ensemble_nameList[0]
            else:
                caseName, ensemble_name, varName = dataName_split
                assert caseName in self.casesDict
            data, dimensions, info = self.casesDict[caseName].read_from_database(ensemble_name, varName)
            assert len(dimensions) == 1 and dimensions[0] == "time"
            if isinstance(data, list):
                data = np.array(data)
            timeList = self.casesDict[caseName].timeList
            plot_kwargs = self.data_setting[dataName]
            ax.plot(timeList, data, **plot_kwargs)
