#This function returns the fully rested open-circuit-voltage of an LiPB
#cell given its soc.
#Syntax: ocv=OCVfromSOCtemp(soc,temp,model)
#where soc is cell state of charge between 0 and 1,
#temp is cell temperature in degrees celsius,
#and model is a cell model structure.

import numpy as np

def OCVfromSOCtemp(soc, temp, model):
    # name = model['model'][0][0][0]
    soccol = np.array([])

    if isinstance(soc, np.ndarray) == False:
        soccol = np.append(soccol,soc)
        # print(soccol, "false")
        
    else:
        soccol = soc[0]
        # print(soccol, "true")

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

    # I4 = np.zeros(len(soccol))

    I4 = (soccol[I3[0]]- SOC[0][0])/diffSOC
    
    I5 = np.floor(I4)
    I5 = I5.astype(int)
    # print(I5)
    # print(ocv[I3[0]].size, I5.size, OCV0.size, I3[0].size, OCVrel[0].size)
    ocv[I3[0]] = np.multiply(OCV0[0][I5],(1-(I4-I5))) + np.multiply(OCV0[0][I5+1],(I4-I5))

    ocv[I3[0]] = ocv[I3[0]] + np.multiply(tempcol[I3[0]],(np.multiply((1 - (I4 - I5)), OCVrel[0][I5]) + np.multiply(OCVrel[0][I5+1],(I4-I5))))
  
    soc[I6[0]] = 0
    # print(ocv[9999])
    # soc = np.reshape(soc,ocv.size);


# OCVfromSOCtemp(soc, 25, model)


