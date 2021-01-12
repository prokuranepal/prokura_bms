from scipy import interpolate


def SOCfromOCVtemp(ocv, temp, model):
    # TODO: implement the SOC - OCV relation
    return 0.0


def OCVfromSOCtemp():
    # TODO: implement the inverse ocv-soc relation
    return 0.0


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