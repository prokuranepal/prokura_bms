from scipy.io import loadmat
import numpy as np
import matplotlib.pyplot as plt


from soc_estimator import init_spkf, iter_spkf


# load the data for the cell
pan_data = loadmat('readonly/PAN_CAPSTONE_DATA.mat')

T = 25  # set the current temperature

# retrieve data
DYNData = pan_data['DYNData']
script1 = DYNData['script1'].item()
time = script1['time'].item().T
deltat = time[1] - time[0]
time = time - time[0]
current = script1['current'].item().T
voltage = script1['voltage'].item().T
soc = script1['soc'].item().T

# load the model for the cell
PANmodel = loadmat('readonly/PANmodel.mat')
soc_hat = np.zeros(soc.shape)
soc_bound = np.zeros(soc.shape)


def tuneSPKF():
    """Set the values for SPKF parametersdatetime A combination of a date and a time. Attributes: ()

    Returns:
        array: Values of SPKF parameters.
    """
    SigmaW = np.array([0.06])
    SigmaV = np.array([0.00078])
    SigmaZ0 = np.array([0.0028])
    return SigmaW, SigmaV, SigmaZ0


SigmaW, SigmaV, SigmaZ0 = tuneSPKF()

SigmaX0 = np.diag(np.hstack([1e-6, 1e-6, SigmaZ0]))

spkfData = init_spkf(float(voltage[0]), T, SigmaX0, SigmaV, SigmaW, PANmodel)

# This simulation tests the SPKF when there is an inital SOC-estimation error
# The true initial SOC is 95%, but we will initialize the SOC estimate in the
# filter to 90% and see how quickly and well the filter converges toward the
# correct SOC.
spkfData['xhat'][spkfData['zkInd']] = 0.90


# Now, enter loop for remainder of time, where we update the SPKF
# once per sample interval
print('Please be patient. This code will take a minute or so to execute.\n')
for k in range(len(voltage)):
    vk = voltage[k]  # "measure" voltage
    ik = current[k]  # "measure" current
    Tk = T          # "measure" temperature
#     print('we are {}th iteration'.format(k))
    # Update SOC (and other model states)
    soc_hat[k], soc_bound[k], spkfData = iter_spkf(
        vk, ik, Tk, deltat, spkfData)
    # update waitbar periodically, but not too often (slow procedure)
    if (k+1) % 300 == 0:
        print('  Completed {} out of {} iterations...\n'.format(k+1, len(voltage)))


# Analyze the result
J1 = np.sqrt(np.mean((100*(soc-soc_hat))**2))
print('RMS soc estimation error: {}'.format(J1))

J2 = 100*soc_bound[-1]
print('Final value of SOC estimation error bounds: {}'.format(J2))

ind = np.argwhere(np.abs(soc-soc_hat) > soc_bound)
print('Percent of time error outside bounds = {}\n'.format(
    len(ind)/len(soc)*100))

# Plot the data
plt.figure(figsize=[16, 9])

plt.subplot(1,2,1) 
plt.plot(time/60,100*soc_hat,time/60,100*soc)
plt.plot(np.concatenate([time/60, [[np.nan]], time/60]),
         np.concatenate([100*(soc_hat+soc_bound), 
                         [[np.nan]], 100*(soc_hat-soc_bound)]), '--')
plt.title('SOC estimation using SPKF')
plt.grid()
plt.xlabel('Time (min)')
plt.ylabel('SOC (%)')
plt.legend(['Estimate','Truth','Bounds']);


plt.subplot(1, 2, 2)
plt.plot(time/60,100*(soc-soc_hat))
plt.plot(np.concatenate([time/60, [[np.nan]], time/60]),
         np.concatenate([100*(soc_bound), 
                         [[np.nan]], -100*(soc_bound)]), '--')
plt.xlabel('Time (min)')
plt.ylabel('SOC errors (%)')
plt.ylim([-4,4])
plt.grid()
plt.title('SOC estimation errors using EKF')
plt.legend(['Estimation error','Bounds'])
plt.show()
