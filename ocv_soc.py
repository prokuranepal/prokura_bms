from scipy import interpolate
import numpy as np


def OCVfromSOCtemp(soc, temp, model):
    """Return the fully rested open-circuit-voltage(OCV)
    of a Li-ion cell given its SOC.

    Args:
        soc (float): cell state of charge between 0 and 1,
        temp (array): temperature in degree celcius
        model (ndarray): cell model structure
    """
    # name = model['model'][0][0][0]
    soccol = np.array([])

    if isinstance(soc, np.ndarray) == False:
        soccol = np.append(soccol, soc)
        # print(soccol, "false")

    else:
        soccol = soc[0]
        # print(soccol, "true")

#     SOC = model['model'][0][0][2]
#     OCV0 = model['model'][0][0][0]
#     OCVrel = model['model'][0][0][1]
    OCV0 = model['model']['OCV0'].item()
#     print(OCV.shape)
    OCVrel = model['model']['OCVrel'].item()
#     print(SOCrel.shape)
    SOC = model['model']['SOC'].item()
    tempcol = temp
    if isinstance(temp, list) == False:
        # replicate temperature for all soc
        tempcol = temp * np.ones(soccol.size)
    ocv = np.zeros(soccol.size)
    # ocv = 0
    diffSOC = SOC[0][1]-SOC[0][0]

    I1 = np.where(soccol <= SOC[0][0])

    I2 = np.where(soccol >= SOC[0][SOC.size - 1])

    I3 = np.where((soccol > SOC[0][0]) & (soccol < SOC[0][SOC.size - 1]))

    I6 = np.isnan(soccol)

    dv = (OCV0[0][1] + np.multiply(soccol, OCVrel[0][1]) -
          (OCV0[0][0] + np.multiply(tempcol, OCVrel[0][0])))

    ocv[I1[0]] = (np.multiply((soccol[I1[0]] - SOC[0][0]), dv[I1[0]]) /
                  diffSOC) + OCV0[0][0] + np.multiply(tempcol[I1[0]], OCVrel[0][0])

    dv = (OCV0[0][OCV0.size - 1]+np.multiply(tempcol, OCVrel[0][OCVrel.size-1])) - \
        (OCV0[0][OCV0.size-2]+np.multiply(tempcol, OCVrel[0][OCVrel.size - 2]))

    ocv[I2[0]] = np.multiply((soccol[I2[0]]-SOC[0][SOC.size - 1]), dv[I2[0]])/diffSOC + \
        OCV0[0][OCV0.size - 1] + \
        np.multiply(tempcol[I2[0]], OCVrel[0][OCVrel.size - 1])

    # I4 = np.zeros(len(soccol))

    I4 = (soccol[I3[0]] - SOC[0][0])/diffSOC

    I5 = np.floor(I4)
    I5 = I5.astype(int)
    # print(I5)
    # print(ocv[I3[0]].size, I5.size, OCV0.size, I3[0].size, OCVrel[0].size)
    ocv[I3[0]] = np.multiply(OCV0[0][I5], (1-(I4-I5))) + \
        np.multiply(OCV0[0][I5+1], (I4-I5))

    ocv[I3[0]] = ocv[I3[0]] + np.multiply(tempcol[I3[0]], (np.multiply(
        (1 - (I4 - I5)), OCVrel[0][I5]) + np.multiply(OCVrel[0][I5+1], (I4-I5))))

    ocv[I6[0]] = 0
    # print(ocv[9999])
    # soc = np.reshape(soc,ocv.size);
    return ocv


# This function returns an estimate of soc from a fully rested open-circuit-voltage
# of an LiPB cell

# import numpy as np


def SOCfromOCVtemp(ocv, temp, model):
    """Return an estimate of SOC from a fully rested OCV
    of a Liion cell.

    Args:
        ocv (float): cell voltage
        temp (array): temperature in degree celcius.
        model (ndarray): cell model structure.
    """
    # name = model['model'][0][0][0]
    ocvcol = np.array([])

    if isinstance(ocv, np.ndarray) == False:
        ocvcol = np.append(ocvcol, ocv)
        # print(ocvcol, "false")

    else:
        ocvcol = ocv[0]
        # print(ocvcol, "true")

    OCV = model['model']['OCV'].item()
#     print(OCV.shape)
    SOCrel = model['model']['SOCrel'].item()
#     print(SOCrel.shape)
    SOC0 = model['model']['SOC0'].item()
    tempcol = temp
#     print(SOC0.shape)
    if isinstance(temp, list) == False:
        # replicate temperature for all ocvs
        tempcol = temp * np.ones(ocvcol.size)
    soc = np.zeros(ocvcol.size)
    diffOCV = OCV[0][1]-OCV[0][0]

    I1 = np.where(ocvcol <= OCV[0][0])
    # print(ocvcol)
    I2 = np.where(ocvcol >= OCV[0][OCV.size - 1])

    I3 = np.where((ocvcol > OCV[0][0]) & (ocvcol < OCV[0][OCV.size - 1]))

    I6 = np.isnan(ocvcol)

    # for socs lower than lowest voltage
    # extrapolate off low end of table
    dz = (SOC0[0][1] + np.multiply(tempcol, SOCrel[0][1]) -
          (SOC0[0][0] + np.multiply(tempcol, SOCrel[0][0])))

    soc[I1[0]] = (np.multiply((ocvcol[I1[0]] - OCV[0][0]), dz[I1[0]]) /
                  diffOCV) + SOC0[0][0] + np.multiply(tempcol[I1[0]], SOCrel[0][0])

    # for socs higher than highest voltage
    # extrapolate off high end of table
    dz = (SOC0[0][SOC0.size - 1]+np.multiply(tempcol, SOCrel[0][SOCrel.size-1])) - \
        (SOC0[0][SOC0.size-2]+np.multiply(tempcol, SOCrel[0][SOCrel.size - 2]))

    soc[I2[0]] = np.multiply((ocvcol[I2[0]]-OCV[0][OCV.size - 1]), dz[I2[0]])/diffOCV + \
        SOC0[0][SOC0.size - 1] + \
        np.multiply(tempcol[I2[0]], SOCrel[0][SOCrel.size - 1])

    # I4 = np.zeros(len(ocvcol))

    # for normal soc range...
    # manually interpolate (10x faster than "interp1")
    I4 = (ocvcol[I3[0]] - OCV[0][0])/diffOCV

    I5 = np.floor(I4)
    I5 = I5.astype(int)

    soc[I3[0]] = np.multiply(SOC0[0][I5], (1-(I4-I5))) + \
        np.multiply(SOC0[0][I5+1], (I4-I5))

    soc[I3[0]] = soc[I3[0]] + np.multiply(tempcol[I3[0]], (np.multiply(
        (1 - (I4 - I5)), SOCrel[0][I5]) + np.multiply(SOCrel[0][I5+1], (I4-I5))))

    soc[I6[0]] = 0
    soc = np.reshape(soc, ocvcol.size)
    return soc


def getParamESC(param_name, temperature, model):
    """Retrieve the cell parameters from the model.

    Args:
        param_name (str): name of the parameter
        temperature (float): temperature where data is req
        model (ndarray): cell data describing all cell parameters

    Raises:
        ValueError: data not found
        ValueError: data not found

    Returns:
        float: the value of the parameter
    """
    theFields = model.dtype.names  # get list of fields stored in model
      # match = param_name in theFields; # see if any match desired data
    fields_arr = np.array(theFields)
    try:
        match = int(np.argwhere((param_name == fields_arr)))
    except:  # if not, throw an error
        raise ValueError(
            'Parameter {} does not exist in model'.format(paramName))

    fieldName = fields_arr[match]  # case-sensitive field name

    # if model contains data at only one temperature
    temps = model['temps'].item()[0]
    if temps.shape[0] == 1:
        if temperature not in temps:  # check whether requested data exists
            raise ValueError('Model does not contain requested\
        data at {} temperature'.format(temperature))
            theParam = model[fieldName]
            return theParam

    # Otherwise, model has multiple temperatures. Bound input "temp" between
    # mininum and maximum stored temperature to prohibit "NaN" in output
    theParamData = model[fieldName].item()
    temp = max(min(temperature, max(temps)), min(temps))
    # ind = find(model.temps == temp);
    ind = np.argwhere(temp == temps)  # see if there is an exact match to

    if ind.size:  # avoid call to (slow) interp1 whenever possible
        # if size(theParamData,1) == 1,
        if min(theParamData.shape) == 1:
            try:
                theParam = theParamData[:, ind]
            except:
                theParam = theParamData[ind, :]
        else:
            theParam = theParamData[ind, :]
            # TODO: Make this work for when shape is not 1
            # TODO: i.e there are multiple values for the same parms
    else:  # if there is not an exact match, we interpolate between parameter
        tck = interpolate.splrep(temps, theParamData.squeeze())  # values
        theParam = interpolate.splev(temp, tck)
    # stored at different temperatures

    return float(theParam)
