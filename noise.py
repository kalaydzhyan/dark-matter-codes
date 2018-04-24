from __future__ import division
""" Importing main functions and classes"""
from main1 import *
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import pyplot as plt
import numpy as np
from scipy.fftpack import dct
import math, sys

""" Some global variables """
number_sources = 1000               # N=1000 gives better results than N=100
observation_time = 500             # Choose between 10, 100, 1000
sampling_time = 1
number_sessions = 10

""" Regime we are probing and angular resolution, in degrees"""
short_wavelength_regime = False
delta_zeta = 10

print "\n\n"    
print "--------------------------------------------------------"
print "     Testing signal with additional detector noise      "
print "--------------------------------------------------------"
print "\nNumber of sources: " + str(number_sources)

""" Creating sources. For reproducible results call Source(some_number) """
sources = Source(0)
for source_number in range(number_sources): 
    sources.add()

""" Getting data for the reference detector """
detector0 = Detector()
""" No noise for the reference detector, because it is correlated """
data = [(t, sources.signal(t, detector0.position) + 0*detector0.noise(t)) \
        for t in np.arange(0, number_sessions*observation_time, sampling_time)]
timeline, signal = zip(*data)

""" Let's add two more detectors """
detector1 = Detector([1, 0, 0])
data1 = [(t, sources.signal(t, detector1.position) + detector1.noise(t)) \
        for t in np.arange(0, number_sessions*observation_time, sampling_time)]

detector2 = Detector([np.sin(0.3), np.cos(0.3), 0])
data2 = [(t, sources.signal(t, detector2.position) + detector2.noise(t)) \
        for t in np.arange(0, number_sessions*observation_time, sampling_time)]

_, signal1 = zip(*data1)
_, signal2 = zip(*data2)
signal1 = np.array(signal1) - np.array(signal)
signal2 = np.array(signal2) - np.array(signal)

zeta = angle(detector1.position, detector2.position)
print 'Angle between two detectors = ' + str(np.rad2deg(zeta)) + ' degrees'

""" Signal for the first detector """
plt.plot(timeline, signal1, '-')
plt.xlabel('t (seconds)')
plt.ylabel('X(t)')
plt.title('Signal for the first detector')
#plt.savefig("signal_reference.png")
plt.show()

""" Spectrum of the signal """
frequencies, fourier_image = zip(*Fourier(data1, sampling_time))
plt.plot(frequencies, np.abs(fourier_image),'r-') 
plt.xlabel('f (Hz)')
plt.ylabel('$\sqrt{S_X(f)}$')
plt.title('Power spectrum of the signal for the first detector')
#plt.savefig("spectrum_reference.png")
plt.show()

""" Compute SIGNAL """

signals = []
printProgressBar(0, number_sessions, prefix = 'Progress:', suffix = 'Complete', length = 50)

for session_i in range(number_sessions):
    signal_array = [[signal1[t1]*signal2[t2]*np.sinc(np.pi*2*frequency_variation*(t1-t2)) \
                         for t1 in np.arange(session_i*observation_time, (session_i + 1)*observation_time)] \
                            for t2 in np.arange(session_i*observation_time, (session_i + 1)*observation_time)]
# 'Delta-function filter for debugging purposes:'
#    signal_array = [signal1[t1]*signal2[t1] \
#                         for t1 in np.arange(session_i*observation_time, (session_i + 1)*observation_time)]
    signals.append((sampling_time**2)*np.sum(signal_array)/observation_time)
    printProgressBar(session_i + 1, number_sessions, prefix = 'Progress:', suffix = 'Complete', length = 50)
    
print "Signal: ", np.mean(signals)
print "Noise: ", np.std(signals)
print "SNR: ", np.abs(np.mean(signals))/np.std(signals)
print "new SNR: ", np.abs(np.mean(signals))/np.sqrt(np.mean([ss**2 for ss in signals]))
