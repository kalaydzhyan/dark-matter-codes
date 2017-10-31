from __future__ import division
from main import *
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import pyplot as plt
import numpy as np
from scipy.fftpack import dct
import math, sys
from statsmodels.nonparametric.smoothers_lowess import lowess

""" Some global variables """
number_sources = 1000
observation_time = 1000     # Should be 1000 or more
sampling_time = 1

""" Regime we are probing and angular resolution, in degrees"""
short_wavelength_regime = False
delta_zeta = 10

print "\n\n"    
print "--------------------------------------------------------"
print " Simulation of stochastic backgrounds of various fields "
print "--------------------------------------------------------"
print "\nNumber of sources: " + str(number_sources)

""" Creating sources. For reproducible results call Source(some_number) """
sources = Source()
for source_number in range(number_sources): 
    sources.add()

""" Plotting the sources """
fig = plt.figure(figsize=(8, 8))
ax = fig.add_subplot(111, projection='3d', aspect=1.0)
ax.scatter(sources.positions(0), sources.positions(1), sources.positions(2), c='r', marker='o')
ax.set_xlabel('X')
ax.set_ylabel('Y')
ax.set_zlabel('Z')
plt.xticks(np.arange(-1, 1, .5))
plt.yticks(np.arange(-1, 1, .5))
plt.title('Angular positions of the sources')
plt.savefig("sources.png")
plt.show()

""" Getting data for the reference detector """
detector0 = Detector()
data = [(t, sources.signal(t, detector0.position)) \
        for t in np.arange(0, observation_time, sampling_time)]

""" Signal for the reference detector """
timeline, signal = zip(*data)
plt.plot(timeline, signal, '-')
plt.xlabel('t (seconds)')
plt.ylabel('X(t)')
plt.title('Signal for the reference detector')
plt.savefig("signal_reference.png")
plt.show()

""" Spectrum of the signal """
frequencies, fourier_image = zip(*Fourier(data, sampling_time))
plt.plot(frequencies, np.abs(fourier_image),'r-') 
plt.xlabel('f (Hz)')
plt.ylabel('$\sqrt{S_X(f)}$')
plt.title('Power spectrum of the signal for the reference detector')
plt.savefig("spectrum_reference.png")
plt.show()

""" Let's add two more detectors """
detector1 = Detector([1, 0, 0])
data1 = [(t, sources.signal(t, detector1.position)) \
        for t in np.arange(0, observation_time, sampling_time)]

detector2 = Detector([np.sin(0.3), np.cos(0.3), 0])
data2 = [(t, sources.signal(t, detector2.position)) \
        for t in np.arange(0, observation_time, sampling_time)]

zeta = angle(detector1.position, detector2.position)
print 'Angle between two detectors: ' + str(np.rad2deg(zeta)) + ' degrees'
print 'Distance between detectors: ', distance(detector1.position, detector2.position)


""" Compute cross-correlation """
_, signal1 = zip(*data1)
_, signal2 = zip(*data2)
signal1 = np.array(signal1) - np.array(signal)
signal2 = np.array(signal2) - np.array(signal)

xcorrelator_full = np.correlate(signal1, signal2, 'full')
xcorrelator = np.array([xcorrelator_full[t] for t in 
        range(np.int(len(xcorrelator_full)/2),len(xcorrelator_full))])
tau = range(0, len(xcorrelator), sampling_time)
averaging_factor = np.array([observation_time - t*sampling_time for t in tau])
xcorrelator /= averaging_factor

""" Cutting boundary effects """    
if (observation_time/sampling_time)>200:
    plt.plot(tau[:-200], xcorrelator[:-200], '-')
else:
    plt.plot(tau, xcorrelator, '-')
    
plt.xlabel('$\\tau$ (seconds)')
plt.ylabel('x-corr(t)')
plt.title('cross-correlation')
plt.savefig("cross_correlator.png")
plt.show()

frequencies, signal1image = zip(*Fourier(zip(timeline, signal1), sampling_time))
frequencies, signal2image = zip(*Fourier(zip(timeline, signal2), sampling_time))

x_fourier_image = signal1image * np.conj(signal2image)

#filtered = lowess(np.real(x_fourier_image), frequencies, is_sorted=True, frac=0.025, it=0)
#plt.plot(frequencies, frequencies, filtered, 'b-')
plt.plot(frequencies, np.real(x_fourier_image),'r-')
plt.xlabel('f (Hz)')
plt.ylabel('$S_c(f)$')
plt.title('Cross-spectrum')
plt.savefig("spectrum_cross.png")
plt.show()


""" Produce the angular curve """
ratios = []

if short_wavelength_regime:
    distance1 = 3000
    distance2 = 300
else:
    distance1 = 0.001
    distance2 = 0.0001

zetas = np.arange(0, math.pi, np.deg2rad(delta_zeta))
timeline = np.arange(0, observation_time, sampling_time)

for zeta in zetas:
    detector1 = Detector([distance1, 0, 0])
    signal1 = [sources.signal(t, detector1.position) \
               for t in timeline]

    detector2 = Detector([distance2 * np.cos(zeta), distance2 * np.sin(zeta), 0])
    signal2 = [sources.signal(t, detector2.position) \
               for t in timeline]

    """ Compare to the reference detector in the origin """
    signal1 = np.array(signal1) - np.array(signal)
    signal2 = np.array(signal2) - np.array(signal)
        
    frequencies, signal1image = zip(*Fourier(zip(timeline, signal1), sampling_time))
    frequencies, signal2image = zip(*Fourier(zip(timeline, signal2), sampling_time))
    
    x_fourier_image = signal1image * np.conj(signal2image)
    
    if frequency_variation*observation_time < 1:
        print 'Frequency variation must be larger than observation time, cannot resolve!'
    
    left_f = int((frequency_mean - frequency_variation) * observation_time)
    right_f = int((frequency_mean + frequency_variation) * observation_time)
    cross_spectrum_f = np.real(np.array(x_fourier_image))
    S_cross = np.mean(cross_spectrum_f[left_f:right_f])
    detector1_spectrum_f = np.abs(np.array(signal1image))**2
    S_d1 = np.mean(detector1_spectrum_f[left_f:right_f])

    ratios.append(S_cross/S_d1)
    
fig, ay = plt.subplots()
ay.plot(zetas, ratios,'go-')

if short_wavelength_regime:
    ratios_mean = np.mean(ratios)
    ratios_std = np.std(ratios)
    print 'Level position: ' + str(ratios_mean) + ' +/- ' + str(ratios_std)
    ay.axhspan(ratios_mean-ratios_std, ratios_mean+ratios_std, alpha=0.5, color='green')
    ay.axhline(y=ratios_mean, linestyle='dashed')

plt.xlabel('$\zeta$, rad')
plt.title('Angular curve')

if short_wavelength_regime:
    plt.savefig("angular_short.png")
else:
    plt.savefig("angular_long.png")
    
plt.show()

