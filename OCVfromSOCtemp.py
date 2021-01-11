#This function returns the fully rested open-circuit-voltage of an LiPB
#cell given its soc.
#Syntax: ocv=OCVfromSOCtemp(soc,temp,model)
#where soc is cell state of charge between 0 and 1,
#temp is cell temperature in degrees celsius,
#and model is a cell model structure.

import scipy.io
import numpy as np
model = scipy.io.loadmat('./soc/readonly/E2model.mat')


# current = data['current']
# voltage = data['voltage']
# soc = data['soc']




def OCVfromSOCtemp(soc, temp, model):
    # name = model['model'][0][0][0]
    soccol = soc
    SOC = model['model'][0][0][2]
    OCV0 = model['model'][0][0][0]
    OCVrel = model['model'][0][0][1]
    tempcol = temp
    if isinstance(temp, list) == False:
        tempcol = temp * np.ones(soccol.size) #replicate temperature for all ocvs
    # ocv=np.zeros(soccol.size)
    ocv = 0
    diffOCV=SOC[0][1]-SOC[0][0]

    I1=np.where(soccol <= SOC[0][0])
 
    I2 = np.where(soccol >= SOC[0][SOC.size - 1])

    I3 = np.where((soccol > SOC[0][0]) & (soccol < OCV[0][OCV.size - 1]))
    
    I6 = np.isnan(soccol) 
    
    # #for socs lower than lowest voltage
    # #extrapolate off low end of table
    dv = (OCV0[0][1] + np.multiply(soccol,OCVrel[0][1]) - (OCV0[0][0] + np.multiply(tempcol,OCVrel[0][0])))
    
    soc[I1[0]] = (np.multiply((soccol[0][I1[0]] - OCV[0][0]),dz[I1[0]])/diffOCV) + SOC0[0][0] + np.multiply(tempcol[I1[0]],SOCrel[0][0])

    # #for socs higher than highest voltage
    # #extrapolate off high end of table
    # dz = (SOC0[0][SOC0.size - 1]+np.multiply(tempcol,SOCrel[0][SOCrel.size-1])) - (SOC0[0][SOC0.size-2]+np.multiply(tempcol,SOCrel[0][SOCrel.size - 2]))
   
    # soc[I2[0]] = np.multiply((ocvcol[0][I2[0]]-OCV[0][OCV.size - 1]),dz[I2[0]])/diffOCV + SOC0[0][SOC0.size -1] + np.multiply(tempcol[I2[0]],SOCrel[0][SOCrel.size - 1])

    # I4 = np.zeros(len(I3[0]))

    # #for normal soc range...
    # #manually interpolate (10x faster than "interp1")
    # I4[I3[0]] = (ocvcol[0]- OCV[0][0])/diffOCV
    
    # I5 = np.floor(I4)
    # I5 = I5.astype(int)
 
    # soc[I3[0]] = np.multiply(SOC0[0][I5],(1-(I4-I5))) + np.multiply(SOC0[0][I5+1],(I4-I5))

    # soc[I3[0]] = soc[I3[0]] + np.multiply(tempcol,(np.multiply((1 - (I4 - I5)), SOCrel[0][I5]) + np.multiply(SOCrel[0][I5+1],(I4-I5))))
  
    # soc[I6[0]] = 0
    
    # soc = np.reshape(soc,ocv.size);


OCVfromSOCtemp(0.5, 25, model)


