"""
    The objective of this exercise is to study the signal of ECG during atrial
    fibrilation (AF). The signal analysed contains different type of AF with
    stable repolarisation loops and random AF.

"""

import numpy as np
import pylab as py
py.ion()
py.close('all')
import scipy.signal as sp

"""
    The first signal is an ECG with atrial fibrilation.

    Q: What are the differences of this ECG with a normal ECG?

"""

ecg = np.genfromtxt('ecg_af.dat')
ecg_fs = 300
t_ecg = np.arange(len(ecg))/ecg_fs

py.figure(1,figsize=[5,5])
py.plot(t_ecg, ecg, 'k')
py.xlabel('time (s)')
py.ylabel('amplitude (a.u.)')
py.title('ECG with atrial fibrilation')
py.xlim(210, 215)

"""
    We compute the autocorelation of the ECG signal.

    In order to dircard the modulation of the baseline we first apply a
    high-pass filter with a cut-off frequency of 0.5 Hz.

    Q: Do you see a specific pattern that permits to characterize the atrial
        fibrilation?

"""

b, a = sp.butter(2, 0.5/ecg_fs*2, btype='high')

ecg_hp = sp.filtfilt(b, a, ecg)

rxx_ecg = np.correlate(ecg_hp, ecg_hp, mode='full') 
k = np.arange(len(rxx_ecg))-len(rxx_ecg)//2

py.figure(2,figsize=[5,5])
py.plot(k, rxx_ecg, 'k')
py.xlabel('$k$')
py.ylabel('amplitude (a.u.)')
py.title('$R_{xx}$ of ECG with atrial fibrilation')
py.xlim(-500, 500)

"""
    Compute the PSD of the ECG signal.

    Q: What do you see?

"""

f, ECG = sp.welch(ecg_hp, nperseg=500, nfft=4096, noverlap=250, fs=ecg_fs)

py.figure(3, figsize=[5,5])
py.clf()
py.plot(f, ECG, 'k')
py.xlabel('frequency (Hz)')
py.ylabel('power (a.u.)')
py.xlim(0,60)

"""
    In order to higlight the signal related to the repolarisation of the
    atria and ECG signal with atrial fibrilation has been process, keeping only
    the P wave (repolarisation of the atria) and the QRST waves have been
    removed.

    During the measurement 4 time segments exhibit different behaviors.

    Q: What are the difference between the different segments ?

"""

p_wave = np.genfromtxt('AF_sync.dat')
p_wave_fs = 50
t_p_wave = np.arange(len(p_wave))/p_wave_fs

segments = [1500, 2000, 2500, 3000, 3500]


py.figure(4,figsize=[10,10])
for n in range(len(segments)-1):
    py.subplot(2, 2, int(n+1))
    idx = np.arange(segments[n], segments[n+1])
    py.plot(t_p_wave[idx], p_wave[idx], 'k')
    py.xlabel('time (s)')
    py.ylabel('amplitude (a.u.)')
    py.title('p_waves for segment '+str(n+1))

"""
    We compute the autocorelation of the p_wave signal.

    In order to dircard the modulation of the baseline we first apply a
    high-pass filter with a cut-off frequency of 0.5 Hz.
    

    Q: Do you see a specific pattern that permits to characterize the atrial
        fibrilation?
    Q: Discuss the organisation of the signals. Which one is the more organised,
        which one is closer to a noise?

"""



py.figure(5,figsize=[10,10])
for n in range(len(segments)-1):
    py.subplot(2, 2, int(n+1))
    idx = np.arange(segments[n], segments[n+1])
    rxx_p_wave = np.correlate(p_wave[idx], p_wave[idx], mode='full') 
    rxx_p_wave /= np.correlate(np.ones(len(idx)), np.ones(len(idx)), mode='full')
    k = np.arange(len(rxx_p_wave))-len(rxx_p_wave)//2
    py.plot(k, rxx_p_wave, 'k')
    py.xlabel('$k$')
    py.ylabel('amplitude (a.u.)')
    py.title('$R_{xx}$ for p_waves for segment '+str(n+1))

"""
    Compute the PSD of the p_wave signal.

    Q: What do you see?
    Q: Which one is the more organised?
    Q: Which ones looks like a noise?
    Q: Which ones exhibit a sustained repolarisation loop?

"""


py.figure(7, figsize=[10,10])
for n in range(len(segments)-1):
    idx = np.arange(segments[n], segments[n+1])
    f, P_WAVE = sp.welch(p_wave[idx], nperseg=250, nfft=4096, noverlap=100, fs=p_wave_fs)
    py.subplot(2, 2, int(n+1))
    py.plot(f, P_WAVE, 'k')
    py.xlabel('frequency (Hz)')
    py.ylabel('power (a.u.)')
    py.xlim(0,25)
    py.title('PSD for p_waves for segment '+str(n+1))


