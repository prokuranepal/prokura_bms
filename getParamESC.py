import scipy.io
import numpy as np
model = scipy.io.loadmat('./soc/readonly/CellModel.mat')
data = scipy.io.loadmat('./soc/readonly/CellData.mat')

time = data['time']
current = data['current']
voltage = data['voltage']
soc = data['soc']


def getParamESC(paramName, temp, model):
    theFields = model['model'].keys()
    print(theFields)


getParamESC(data, data, model)