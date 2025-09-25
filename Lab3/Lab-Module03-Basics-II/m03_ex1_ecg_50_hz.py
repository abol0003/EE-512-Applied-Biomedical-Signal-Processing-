"""
    The objective of this exercise is to study the influence of the
    parameterization of the Welch spectral estimator in order to highlight a 50
    Hz perturbation in an ECG signal.

"""

import numpy as np
import pylab as py
py.ion()
py.close('all')
import scipy.signal as sp

x = np.genfromtxt('ecg.dat')
fs = 500

"""
    Objective: Compare spectral estimation for different window lengths using
    welch estimation.

    Plot the log spectrum of the signal using windows of 100, 500, 2000.

    Q: Comment the results.
    Q: Which windows length is the most suitable for the observation of 50 Hz?
    Q: Why?
"""

f,X_100 = sp.welch(x, nperseg=100, nfft=4096, fs=fs)
f,X_500 = sp.welch(x, nperseg=500, nfft=4096, fs=fs)
f,X_2000 = sp.welch(x, nperseg=2000, nfft=4096, fs=fs)

py.figure(1, figsize=[5,8])
py.clf()
py.subplot(3,1,1)
py.plot(f, 10*np.log10(X_100), 'k')
py.grid()
py.xlabel('frequency (Hz)')
py.ylabel('power (dB)')
py.title('PSD with lengh = 100')

py.subplot(3,1,2)
py.plot(f, 10*np.log10(X_500), 'k')
py.grid()
py.xlabel('frequency (Hz)')
py.ylabel('power (dB)')
py.title('PSD with lengh = 500')

py.subplot(3,1,3)
py.plot(f, 10*np.log10(X_2000), 'k')
py.grid()
py.xlabel('frequency (Hz)')
py.ylabel('power (dB)')
py.title('PSD with lengh = 2000')




