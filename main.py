import scipy.io
import SOCfromOCVtemp

model = scipy.io.loadmat('./soc/readonly/CellModel.mat')
data = scipy.io.loadmat('./soc/readonly/CellData.mat')

time = data['time']
current = data['current']
voltage = data['voltage']
soc = data['soc']

d = SOCfromOCVtemp.SOCfromOCVtemp(voltage, 25, model)