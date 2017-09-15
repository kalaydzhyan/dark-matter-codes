from __future__ import division
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import pyplot as plt
import numpy as np
from scipy.fftpack import dct
import math, sys
from statsmodels.nonparametric.smoothers_lowess import lowess

""" Some global variables """
mass = 1 # mass of the scalar field, add it later
number_sources = 1000
observation_time = 1000
sampling_time = 1

""" Properties of the source frequencies """
frequency_mean = 0.1
frequency_variation = 0.01

""" Regime we are probing and angular resolution, in degrees"""
short_wavelength_regime = True
delta_zeta = 1


def normalize(vector):
    norm = math.sqrt(math.fsum([v**2 for v in vector]))
    if norm == 0:
        print "Null vector! Can't normalize."
        
    return [element/norm for element in vector]
     
# Source class 
class Source:
    
    def __init__(self, values=None):

        """list of all sources"""        
        self.list = []
        
        """ Random seed for testing purposes """
        if values is not None:
            np.random.seed(values)
        
    def __repr__(self):
        """string representation of the object"""
        representation = ""
        
        for one_source in self.list:
            representation += "\n Source properties." + \
            "\n\t name: \t\t" + one_source['name'] + \
            "\n\t amplitude: \t" + str(one_source['amplitude']) + \
            "\n\t frequency: \t" + str(one_source['frequency']) + \
            "\n\t wave vactor: \t" + str(one_source['wave vector']) + \
            "\n\t phase: \t" + str(one_source['phase']) + "\n"
            
        return representation
        
        
    def add(self, set_name = None):
        
        if set_name is None:
            name = str(len(self.list))
        else:
            name = set_name
            
        amplitude   = 10 * np.random.random()
        frequency   = np.random.normal(frequency_mean, frequency_variation)
        sin_delta   = (np.random.random() - 0.5) * 2
        alpha       = 2 * math.pi * np.random.random()
        phase       = 2 * math.pi * np.random.random()
        
        cos_delta   = math.sqrt(1 - sin_delta ** 2)
        wave_vector = [cos_delta * math.cos(alpha), \
                            cos_delta * math.sin(alpha), sin_delta]
        wave_vector = normalize(wave_vector)
        
        self.list.append({'name': name,
                          'amplitude': amplitude,
                          'frequency': frequency,
                          'wave vector': wave_vector,
                          'phase': phase
                          })

    def positions(self, i):
        """return coordinates of sources for, e.g., visualization"""
        try: 
            return [one_source['wave vector'][i] for one_source in self.list]
        except IndexError:
            sys.exit("Positions of sources: asked for a nonexistent coordinate")
            
        
    def signal(self, t, x):
        
        background = 0.0
        
        for one_source in self.list:
            background += one_source['amplitude'] * \
            math.sin(2*math.pi*one_source['frequency'] * t +\
                     np.inner(one_source['wave vector'], x) +\
                     one_source['phase'])
        
        return background

# Detector class     
class Detector:
    
    def __init__(self, values=None):
        
        """ Detector position """
        if values is None:
            self.position = [0 for _ in range(3)]
        else:
            self.position = values

def Fourier(data, sampling_interval):
    """ Fourier transform, conventions: X_k = Sum_m[x_m*Exp[-2 PI i m k / N]] """
    t, y = zip(*data)

    N = len(y)                  # number of data points
    k = np.arange(N)
    T = N * sampling_interval
    frq = k/T                   # total frequency range
    Y = np.fft.rfft(y)/N        # real FFT computing and normalization

    return zip(frq, Y)

def angle(X, Y):
    """ Angle between two vectors """
    return np.arccos(np.inner(normalize(X), normalize(Y)))



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
print 'Angle between two detectors = ' + str(np.rad2deg(zeta)) + ' degrees'

""" Compute cross-correlation """
_, signal1 = zip(*data1)
_, signal2 = zip(*data2)
signal1 = np.array(signal1) - np.array(signal)
signal2 = np.array(signal2) - np.array(signal)
xcorrelator = np.correlate(signal1, signal2, 'same')/observation_time
# figure out what is the appropriate regime here

tau = range(0, len(xcorrelator), sampling_time)

plt.plot(tau, xcorrelator, '-')
plt.xlabel('t (seconds)')
plt.ylabel('x-corr(t)')
plt.title('correlation')
plt.savefig("cross_correlator.png")
plt.show()

#frequencies, x_fourier_image = zip(*Fourier(zip(tau, xcorrelator), sampling_time))
#_, d1_fourier_image = zip(*Fourier(zip(timeline, signal1), sampling_time))
#filtered = lowess(x_fourier_image, frequencies, is_sorted=True, frac=0.025, it=0)
#plt.plot(frequencies, x_fourier_image,'r-', frequencies, filtered, 'b-')
frequencies, signal1image = zip(*Fourier(zip(timeline, signal1), sampling_time))
frequencies, signal2image = zip(*Fourier(zip(timeline, signal2), sampling_time))

x_fourier_image = signal1image * np.conj(signal2image)

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

