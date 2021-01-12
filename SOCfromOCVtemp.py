# This function returns an estimate of soc from a fully rested open-circuit-voltage 
# of an LiPB cell

import numpy as np

def SOCfromOCVtemp(ocv, temp, model):
    # name = model['model'][0][0][0]
    ocvcol = np.array([])

    if isinstance(ocv, np.ndarray) == False:
        ocvcol = np.append(ocvcol,ocv)
        # print(ocvcol, "false")
        
    else:
        ocvcol = ocv[0]
        # print(ocvcol, "true")
    
    OCV = model['model'][0][0][4]
    SOCrel = model['model'][0][0][6]
    SOC0 = model['model'][0][0][5]
    tempcol = temp
    if isinstance(temp, list) == False:
        tempcol = temp * np.ones(ocvcol.size) #replicate temperature for all ocvs
    soc=np.zeros(ocvcol.size)
    diffOCV=OCV[0][1]-OCV[0][0]

    I1=np.where(ocvcol <= OCV[0][0])
    # print(ocvcol)
    I2 = np.where(ocvcol >= OCV[0][OCV.size - 1])

    I3 = np.where((ocvcol > OCV[0][0]) & (ocvcol < OCV[0][OCV.size - 1]))
    
    I6 = np.isnan(ocvcol) 
    
    #for socs lower than lowest voltage
    #extrapolate off low end of table
    dz = (SOC0[0][1] + np.multiply(tempcol,SOCrel[0][1]) - (SOC0[0][0] + np.multiply(tempcol,SOCrel[0][0])))
    
    soc[I1[0]] = (np.multiply((ocvcol[I1[0]] - OCV[0][0]),dz[I1[0]])/diffOCV) + SOC0[0][0] + np.multiply(tempcol[I1[0]],SOCrel[0][0])

    #for socs higher than highest voltage
    #extrapolate off high end of table
    dz = (SOC0[0][SOC0.size - 1]+np.multiply(tempcol,SOCrel[0][SOCrel.size-1])) - (SOC0[0][SOC0.size-2]+np.multiply(tempcol,SOCrel[0][SOCrel.size - 2]))
   
    soc[I2[0]] = np.multiply((ocvcol[I2[0]]-OCV[0][OCV.size - 1]),dz[I2[0]])/diffOCV + SOC0[0][SOC0.size -1] + np.multiply(tempcol[I2[0]],SOCrel[0][SOCrel.size - 1])

    # I4 = np.zeros(len(ocvcol))

    #for normal soc range...
    #manually interpolate (10x faster than "interp1")
    I4 = (ocvcol[I3[0]]- OCV[0][0])/diffOCV
    
    I5 = np.floor(I4)
    I5 = I5.astype(int)
 
    soc[I3[0]] = np.multiply(SOC0[0][I5],(1-(I4-I5))) + np.multiply(SOC0[0][I5+1],(I4-I5))

    soc[I3[0]] = soc[I3[0]] + np.multiply(tempcol[I3[0]],(np.multiply((1 - (I4 - I5)), SOCrel[0][I5]) + np.multiply(SOCrel[0][I5+1],(I4-I5))))
  
    soc[I6[0]] = 0
    soc = np.reshape(soc,ocvcol.size)
    return soc



