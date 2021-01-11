# This function returns an estimate of soc from a fully rested open-circuit-voltage 
# of an LiPB cell

import scipy.io
import numpy as np
model = scipy.io.loadmat('./soc/readonly/CellModel.mat')
data = scipy.io.loadmat('./soc/readonly/CellData.mat')

time = data['time']
current = data['current']
voltage = data['voltage']
soc = data['soc']




def SOCfromOCVtemp(ocv, temp, model):
    # name = model['model'][0][0][0]
    ocvcol = ocv
    OCV = model['model'][0][0][4]
    SOCrel = model['model'][0][0][6]
    SOC0 = model['model'][0][0][5]
    tempcol = temp
    if isinstance(temp, list) == False:
        tempcol = temp * np.ones(ocv.size)
    soc=np.zeros(ocvcol.size)
    diffOCV=OCV[0][1]-OCV[0][0]

    I1=np.where(ocvcol <= OCV[0][0])
    # I1 = [value for value in ocvcol if (ocvcol[ocvcol <= OCV[0][0]])]
    I2 = np.where(ocvcol >= OCV[0][OCV.size - 1])
    # I2= [value for value in ocvcol if (ocvcol[ocvcol >= OCV[0][OCV.size - 1]])]
    I3 = np.where((ocvcol[0] > OCV[0][0]) & (ocvcol[0] < OCV[0][OCV.size - 1]))
    # I3 = [value for value in ocvcol if (ocvcol[(ocvcol > OCV[0][0]) & (ocvcol < OCV[0][OCV.size - 1])]).any()]
    I6 = np.isnan(ocvcol)
    # print(I6)
   
    
    dz = (SOC0[0][1] + np.multiply(tempcol,SOCrel[0][1]) - (SOC0[0][0] + np.multiply(tempcol,SOCrel[0][0])))
    # print(dz)
    
    soc[I1[0]] = (np.multiply((ocvcol[0][I1[0]] - OCV[0][0]),dz[I1[0]])/diffOCV) + SOC0[0][0] + np.multiply(tempcol[I1[0]],SOCrel[0][0])
    # while i < (len(I1[0])):
    # for j in I1[0]:
    #     soc[j] = (np.multiply((ocvcol[j] - OCV[0][0]),dz[j])/diffOCV) + SOC0[0][0] + np.multiply(tempcol[j],SOCrel[0][0])
        # i = i + 1
    # print(soc)
    #for socs higher than highest voltage
    #extrapolate off high end of table
    dz = (SOC0[0][SOC0.size - 1]+np.multiply(tempcol,SOCrel[0][SOCrel.size-1])) - (SOC0[0][SOC0.size-2]+np.multiply(tempcol,SOCrel[0][SOCrel.size - 2]))
    # print(dz.size, ocvcol[0].size,tempcol.size)
   
    soc[I2[0]] = np.multiply((ocvcol[0][I2[0]]-OCV[0][OCV.size - 1]),dz[I2[0]])/diffOCV + SOC0[0][SOC0.size -1] + np.multiply(tempcol[I2[0]],SOCrel[0][SOCrel.size - 1])
    # while i < (len(I2[0])):
    # for j in I2[0]:
    #     soc[j] = np.multiply((ocvcol[0][j]-OCV[0][OCV.size - 1]),dz[j])/diffOCV + SOC0[0][SOC0.size -1] + np.multiply(tempcol[j],SOCrel[0][SOCrel.size - 1])
        # i = i + 1

    I4 = np.zeros(len(I3[0]))
    # print(len(I3[0]), ocvcol[0].size, OCV[0].size)

    I4[I3[0]] = (ocvcol[0]- OCV[0][0])/diffOCV
    # while i < (len(I3[0])):
    # for j in I3[0]:
    #     I4.append((ocvcol[0][j]- OCV[0][0])/diffOCV)
        # i = i + 1
    I5 = np.floor(I4)
    I5 = I5.astype(int)
    print(I5[4567])

    # while i < (len(I3[0])):
    soc[I3[0]] = np.multiply(SOC0[0][I5],(1-(I4-I5))) + np.multiply(SOC0[0][I5+1],(I4-I5))
    # for j in I3[0]:
        # soc[I3[0][i]] = np.multiply(SOC0[0][I5[i] + 1],(1-(I4[i]-I5[i]))) + np.multiply(SOC0[0][I5[i]+2],(I4[i]-I5[i]))
        # (SOC0[0][I5[j]+1],(1-(I4[j]-I5[j]))) 
        # soc[j] =  np.multiply(SOC0[0][I5[j] + 1],(1-(I4[j]-I5[j]))) + np.multiply(SOC0[0][I5[j]+2],(I4[j]-I5[j]))
        # i = i + 1

    print(soc[0])

    soc[I3[0]] = soc[I3[0]] + np.multiply(tempcol,(np.multiply((1 - (I4 - I5)), SOCrel[0][I5]) + np.multiply(SOCrel[0][I5+1],(I4-I5))))
    # while i < (len(I3[0])):
    # for j in I3[0]:
        # soc[I3[0][i]] = soc[I3[0][i]] + np.multiply(np.multiply(tempcol[I3[0][i]],SOCrel[0][I5[i]+1]),(1 - (I4[i] - I5[i]))) + np.multiply(SOCrel[0][I5[i]+2],(I4[i]-I5[i]))
        # soc[j] = soc[j] + np.multiply(np.multiply(tempcol[j], SOCrel[0][I5[j]+1]),(1 - (I4[j] - I5[j]))) + np.multiply(SOCrel[0][I5[j]+2],(I4[j]-I5[j]))
        # i = i + 1
    
    soc[I6[0]] = 0
    # for j in I6[0]:
    #     soc[j] = 0
    
    # soc = np.reshape(soc,ocv.size);
    print(soc[0])
    

# for socs lower than lowest voltage
# extrapolate off low end of table

    # dz = SOC0[0] + np.multiply(tempcol[:], SOCrel[:])
    
        


SOCfromOCVtemp(voltage, 25, model)
