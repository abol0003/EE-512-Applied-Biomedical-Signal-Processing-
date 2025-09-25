"""
    The objective of this exercice is to analyse the control of the autonomic
    nervous system at rest and after alcool consumption using breathing, mean
    blood pressure and interbeat signals.

"""

import numpy as np
import pylab as py
py.ion()
py.close('all')
import scipy.signal as sp
import m03_ex2_ext as my_plot

# Load signals of a subject at rest.
x = np.genfromtxt('heart_1.dat', delimiter='  ').T
x = {'rr':x[0], 'bp':x[1], 'resp':x[2]}
# Load signals of a subject after alcool consumption.
y = np.genfromtxt('heart_2.dat', delimiter='  ').T
y = {'rr':y[0], 'bp':y[1], 'resp':y[2]}
# Signals are sampled at 4 Hz.
fs = 4
# Generate the time for the recordings.
t = np.arange(len(x['rr']))/fs

"""
    Cardiac interbeats, mean blood pressure and respiration volume of a subject
    at rest.

    Q: Comment the different signals and their realtionships.
    Q: Which signals are related and how?

"""


my_plot.plot_time(x, t, 'normal state')

"""
    Cardiac interbeats, mean blood pressure and respiration volume of a subject
    after alcool consumption.

    Q: Comment the different signals and their realtionships.
    Q: Which signals are related and how?
    Q: What are the differences with rest recording of previous figure?

"""

my_plot.plot_time(y, t, 'after alcool consumption')

""" 
    Compute the intercorrelation of the signals of the subject at rest.

    Q: Comment the oscillation present in the different signals.
    Q: Which signals are related and how.

"""

def my_corr(x):
    rxx = np.correlate(x-np.mean(x), x-np.mean(x), mode='full')/len(x)
    return rxx

x['rxx_rr'] = my_corr(x['rr'])
x['rxx_bp'] = my_corr(x['bp'])
x['rxx_resp'] = my_corr(x['resp'])

my_plot.plot_rxx(x, 'Rxx for normal state')

""" 
    Compute the intercorrelation of the signals of the subject after alcool
    consumption.

    Q: Comment the oscillation present in the different signals.
    Q: Which signals are related and how.
    Q: What difference do you observe with the previous figure?

"""

y['rxx_rr'] = my_corr(y['rr'])
y['rxx_bp'] = my_corr(y['bp'])
y['rxx_resp'] = my_corr(y['resp'])

my_plot.plot_rxx(y, 'Rxx after alcool consumption')


""" 
    Compute the PSD of the signal for the subject at rest.

    Q: How the different peaks are related to the control of the autonomic
        nervous system?
    Q: Do the positions and amplitude of the peaks confirm you previous
        findings;

"""

def my_psd(x, half_win=250):
    interval = np.arange(-half_win,half_win+1)+len(x)//2
    x_sub = x[interval]
    psd = np.abs(np.fft.fft(sp.hann(len(interval))*x_sub, 2048))
    return psd

x['RR'] = my_psd(x['rxx_rr'])
x['BP'] = my_psd(x['rxx_bp'])
x['RESP'] = my_psd(x['rxx_resp'])

my_plot.plot_X(x, fs, 'PSD for normal state')

""" 
    Compute the PSD of the signal for the subject after alcool consumption.

    Q: How the different peaks are related to the control of the autonomic
        nervous system?
    Q: Do the positions and amplitude of the peaks confirm you previous
        findings;

"""

y['RR'] = my_psd(y['rxx_rr'])
y['BP'] = my_psd(y['rxx_bp'])
y['RESP'] = my_psd(y['rxx_resp'])

my_plot.plot_X(y, fs, 'PSD after alcool consumption')

"""
    Plot the PSDs of the signals for the two conditions.

    Q: Discuss the differences.

"""

my_plot.plot_XY(x, y, fs, 'Comparison of the PSD')
