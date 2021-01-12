import scipy.io
import SOCfromOCVtemp
import OCVfromSOCtemp

OCVmodel = scipy.io.loadmat('./soc/readonly/CellModel.mat')
data = scipy.io.loadmat('./soc/readonly/CellData.mat')

SOCmodel = scipy.io.loadmat('./soc/readonly/E2model.mat')
data = scipy.io.loadmat('./soc/readonly/CellData.mat')
soc = data['soc']

time = data['time']
current = data['current']
voltage = data['voltage']
soc = data['soc']

d = SOCfromOCVtemp.SOCfromOCVtemp(voltage, 25, OCVmodel)
e = OCVfromSOCtemp.OCVfromSOCtemp(soc, 25, SOCmodel)