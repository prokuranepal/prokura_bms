import numpy as np
from ocv_soc import SOCfromOCVtemp, getParamESC
# from scipy import linalg as LA


def init_spkf(v0, T0, SigmaX0, SigmaV, SigmaW, model):
    # initial state description
    ir0 = 0.0
    hk0 = 0.0
    SOC0 = SOCfromOCVtemp(v0, T0, model)

    # initialize data dictionary
    spkf_data = {}
    spkf_data['irInd'] = 1
    spkf_data['hkInd'] = 2
    spkf_data['zkInd'] = 3
    # initial state
    spkf_data['xhat'] = np.array((ir0, hk0, SOC0), ndmin=2).T
    spkf_data['SigmaX'] = SigmaX0
    spkf_data['SigmaV'] = SigmaV
    spkf_data['SigmaW'] = SigmaW
    spkf_data['Snoise'] = np.real(np.linalg.cholesky(
        np.diag([SigmaW, SigmaV])
    ))
    spkf_data['Qbump'] = 5

    # SPKF specific parameters
    Nx = len(spkf_data['xhat'])  # state-vector length
    spkf_data['Nx'] = Nx
    Ny = 1  # measurement-vector length
    spkf_data['Ny'] = Ny
    Nu = 1  # input-vector length
    spkf_data['Nu'] = Nu
    Nw = SigmaW.shape[0]  # process-noise-vector length
    spkf_data['Nw'] = Nw
    Nv = SigmaV.shape[0]  # sensor-noise-vector length
    spkf_data['Nv'] = Nv
    Na = Nx+Nw+Nv  # augmented-state-vector length
    spkf_data['Na'] = Na

    h = np.sqrt(3)
    h = 3
    # h = 3 why do this? I don't understand
    spkf_data['h'] = h  # spkf tuning factor
    weight1 = (h*h - Na)/(h*h)  # weighting factors when computing mean
    weight2 = 1/(2*h*h)  # and for covariance

    spkf_data['Wm'] = np.array([weight1, weight2*np.ones(2*Na, 1)])
    spkf_data['Wc'] = spkf_data['Wm']

    # previous value of current
    spkf_data['priorI'] = 0
    spkf_data['signIk'] = 0

    # store model data structure too
    spkf_data['model'] = model

    return spkf_data


def iter_spkf(vk, ik, Tk, deltat, spkf_data):
    model = spkf_data['model']

    # load all cell model parameters
    Q = getParamESC('QParam', Tk, model)
    G = getParamESC('GParam', Tk, model)
    M = getParamESC('MParam', Tk, model)
    M0 = getParamESC('M0Param', Tk, model)
    RC = np.exp(-deltat/abs(getParamESC()))
    R = getParamESC('RParam', Tk, model)
    R0 = getParamESC('R0Param', Tk, model)
    eta = getParamESC('etaParam', Tk, model)
    if ik < 0:
        ik = ik*eta

    # get data stored in spkf_data dict
    I = spkf_data['priorI']
    SigmaX = spkf_data['SigmaX']
    xhat = spkf_data['xhat']
    Nx = spkf_data['Nx']
    Nw = spkf_data['Nw']
    Nv = spkf_data['Nv']
    Na = spkf_data['Na']
    Snoise = spkf_data['Snoise']
    Wc = spkf_data['Wc']
    irInd = spkf_data['irInd']
    hkInd = spkf_data['hkInd']
    zkInd = spkf_data['zkInd']

    if np.abs(ik) > Q/100:
        spkf_data['signIk'] = np.sign(ik)

    signIk = spkf_data['signIk']

    # step 1a: State estimate time update
    # - create xhatminus augmented sigmax points
    # - extract xhatminus state sigmax points
    # - compute weighted average xhatminus(k)

    # step 1a-1: create augmented sigmax and xhat
    try:
        sigmaXa = np.linalg.chol(SigmaX)
    except np.linalg.LinAlgError as err:
        if 'Matrix is not positive definite' in str(err):
            print('Cholesky error: ' + str(err))
            theAbsDiag = np.abs(np.diag(SigmaX))
            sigmaXa = np.diag([np.max(SQRT(theAbsDiag)),
                               SQRT(spkf_data['SigmaW'])])
        else:
            raise

    sigmaXa = np.array([[np.real(sigmaXa), np.zeros((Nx, Nw+Nv))],
                        [np.zeros((Nw+Nv, Nx)), Snoise]])

    xhata = np.vstack([xhat, np.zeros((Nw+Nv, 1))])
    # NOTE: sigmaXa is lower-triangular

    # Step 1a-2: Calculate SigmaX points
    Xa = xhata[:, np.ones((1, 2*Na + 1), dtype='int')] + \
        spkf_data['h'] * np.vstack([np.zeros((Na, 1)),
                                    sigmaXa,
                                    -sigmaXa])
    # strange indexing to avoid repmat call
    # which is inefficient operation

    # Step 1a-3: Time update from last iteration until now
    # stateEqn(xold, current, xnoise)
    Xx = stateEqn(Xa[1:Nx, :], I, Xa[Nx+1:Nx+Nw, :])
    xhat = Xx*spkf_data['Wm']

    # Step 1b: Error covariance time update
    # Compute weihgted covariance sigmaminus(k)
    Xs = Xx - xhat[:, np.ones((1, 2*Na + 1), dtype='int')]
    SigmaX = Xs * np.diag(Wc) * Xs.T

    # Step 1c: Output Estimate
    # Compute weighted output estimate yhat(k)
    I = ik
    yk = vk
    Y = outputEqn(Xx, I, Xa[Nx+Nw+1:, :], Tk, model)
    yhat = Y*spkf_data['Wm']

    # Step 2a: Estimator gain matrix
    Ys = Y - yhat[:, np.ones((1, 2*Na+1), dtype='int')]
    SigmaXY = Xs*np.diag(Wc)*Ys.T
    SigmaY = Ys*np.diag(Wc)*Ys.T
    L = SigmaXY / SigmaY

    # Step 2b: State estimate measurement update
    r = yk - yhat  # residual. used to check sensor errors
    if r**2 > 100*SigmaY:
        L[:, 1] = 0.0

    xhat = xhat + L*r
    xhat[zkInd] = min([1.05, max([-0.05, xhat(zkInd)])])
    xhat[hkInd] = min([1, max([-1, xhat[hkInd]])])

    # Step 2c: Error covariance measurement update
    SigmaX = SigmaX - L*SigmaY*L
    _, S, V = np.linalg.svd(SigmaX)
    V = V.T
    HH = V*S*V.T
    SigmaX = (SigmaX + SigmaX.T + HH + HH.T) / 4  # help maintain
    # robustness

    # Q-bump code
    if r**2 > 4*SigmaY:
        print('Bumping Sigmax\n')
        SigmaX(zkInd, zkInd) = SigmaX(zkInd, zkInd) * spkf_data[
            Qbump]

    # Save data in spkf_data structure for next itern.
    spkf_data[priorI] = ik
    spkf_data[SigmaX] = SigmaX
    spkf_data[xhat] = xhata
    zk = xhat(zkInd)
    zkbnd = 3 * np.sqrt(SigmaX[zkInd, zkInd])

    # calculate new states for all of the old vectors
    def stateEqn(xold, current, xnoise):
        current = current + xnoise  # add noise to current
        xnew = 0*xold
        xnew[irInd, :] = (RC @ xold[irInd, :] +
                          (1-np.diag(RC)) @ current)
        Ah = exp(-abs(current*G*deltat/(3600*Q)))
        xnew[hkInd, :] = Ah*xold(hkInd,:) + (Ah-1)*np.sign(current)
        xnew[zkInd, :] = xold[zkInd,:] - current/3600/Q
        xnew[hkInd, :] = min([1, max([-1,xnew[hkInd,:]])])
        xnew[zkInd, :] = min([1.05, max([-0.05, xnew[zkInd, :]])])
        return xnew

    # calculate cell output voltage for all state vectors
    def outputEqn(xhat, current, ynoise, T, model):
        yhat = OCVfromSOCtemp(xhat[zkInd,:], T, model)
        yhat = yhat + M*xhat[hkInd,:] + M0*signIk
        yhat = yhat - R*xhat[irInd,:] - R0*current \
            + ynoise[1,:]
        return yhat

    # "Safe" square root
    def SQRT(x):
        X = np.sqrt(max(0,x))
        return X




