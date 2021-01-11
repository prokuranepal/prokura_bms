#This function returns the fully rested open-circuit-voltage of an LiPB
#cell given its soc.
#Syntax: ocv=OCVfromSOCtemp(soc,temp,model)
#where soc is cell state of charge between 0 and 1,
#temp is cell temperature in degrees celsius,
#and model is a cell model structure.

import scipy.io
import numpy as np
model = scipy.io.loadmat('./soc/readonly/E2model.mat')
data = scipy.io.loadmat('./soc/readonly/CellData.mat')
soc = data['soc']
def OCVfromSOCtemp(soc, temp, model):
    # name = model['model'][0][0][0]
    soccol = soc[0]
    SOC = model['model'][0][0][2]
    OCV0 = model['model'][0][0][0]
    OCVrel = model['model'][0][0][1]
    tempcol = temp
    if isinstance(temp, list) == False:
        tempcol = temp * np.ones(soccol.size) #replicate temperature for all soc
    ocv=np.zeros(soccol.size)
    # ocv = 0
    diffSOC=SOC[0][1]-SOC[0][0]

    I1=np.where(soccol <= SOC[0][0])
 
    I2 = np.where(soccol >= SOC[0][SOC.size - 1])

    I3 = np.where((soccol > SOC[0][0]) & (soccol < SOC[0][SOC.size - 1]))
    
    I6 = np.isnan(soccol) 
    
  
    dv = (OCV0[0][1] + np.multiply(soccol,OCVrel[0][1]) - (OCV0[0][0] + np.multiply(tempcol,OCVrel[0][0])))
    
    ocv[I1[0]] = (np.multiply((soccol[I1[0]] - SOC[0][0]),dv[I1[0]])/diffSOC) + OCV0[0][0] + np.multiply(tempcol[I1[0]],OCVrel[0][0])

   
    dv = (OCV0[0][OCV0.size - 1]+np.multiply(tempcol,OCVrel[0][OCVrel.size-1])) - (OCV0[0][OCV0.size-2]+np.multiply(tempcol,OCVrel[0][OCVrel.size - 2]))
   
    ocv[I2[0]] = np.multiply((soccol[I2[0]]-SOC[0][SOC.size - 1]),dv[I2[0]])/diffSOC + OCV0[0][OCV0.size -1] + np.multiply(tempcol[I2[0]],OCVrel[0][OCVrel.size - 1])

    I4 = np.zeros(len(I3[0]))

    
    I4[I3[0]] = (soccol- SOC[0][0])/diffSOC
    
    I5 = np.floor(I4)
    I5 = I5.astype(int)
 
    ocv[I3[0]] = np.multiply(OCV0[0][I5],(1-(I4-I5))) + np.multiply(OCV0[0][I5+1],(I4-I5))

    ocv[I3[0]] = ocv[I3[0]] + np.multiply(tempcol,(np.multiply((1 - (I4 - I5)), OCVrel[0][I5]) + np.multiply(OCVrel[0][I5+1],(I4-I5))))
  
    # soc[I6[0]] = 0
    print(ocv)
    # soc = np.reshape(soc,ocv.size);


OCVfromSOCtemp(soc, 25, model)


