from __future__ import division
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import pyplot as plt
import numpy as np
import math, random, sys


""" Some global variables """
mass = 1 # mass of the scalar field, add it later
number_sources = 100
observation_time = 1000
sampling_time = 1

def normalize(vector):
    norm = math.sqrt(math.fsum([v**2 for v in vector]))
    return [element/norm for element in vector] 

# Source class 
class Source:
    
    def __init__(self, values=None):

        """list of all sources"""        
        self.list = []
        
        """ Random seed for testing purposes """
        if values is not None:
            random.seed(values)
        
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
            
        amplitude   = random.random()
        frequency   = random.random()/2
        sin_delta   = (random.random() - 0.5) * 2
        alpha       = 2 * math.pi * random.random()
        phase       = 2 * math.pi * random.random()
        
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
    
#    def plot sphere of sources
#    def remove source by name

def Fourier(data, sampling_interval):
    
    t, y = zip(*data)

    N = len(y)                  # number of data points
    k = np.arange(N)
    T = N * sampling_interval
    frq = k/T                   # total frequency range
    frq = frq[range(N//2)]      # cut by Kotel'nikov

    Y = np.fft.fft(y)/N         # FFT computing and normalization
    Y = Y[range(N//2)]
    
    return zip(frq, abs(Y))

print "\n\n"    
print "--------------------------------------------------------"
print " Simulation of stochastic backgrounds of various fields "
print "--------------------------------------------------------"
print "\nNumber of sources: " + str(number_sources)

""" Creating sources. For reproducible results call Source(0) """
sources = Source()
for source_number in range(number_sources): 
    sources.add()

""" Plotting the sources """
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.scatter(sources.positions(0), sources.positions(1), sources.positions(2), c='r', marker='o')
ax.set_xlabel('X')
ax.set_ylabel('Y')
ax.set_zlabel('Z')
plt.xticks(np.arange(-1, 1, .5))
plt.yticks(np.arange(-1, 1, .5))
plt.title('Positions of sources')
plt.show()

""" Getting data for one detector """
data = [(t, sources.signal(t, [0, 0, 1])) \
        for t in np.arange(0, observation_time, sampling_time)]

""" Signal for one detector """
timeline, signal = zip(*data)
plt.plot(timeline, signal, '-')
plt.xlabel('Time (s)')
plt.ylabel('signal')
plt.title('Signal for one detector')
plt.show()

""" Spectrum of the signal """
frequencies, fourier_image = zip(*Fourier(data, sampling_time))
plt.plot(frequencies, fourier_image,'r-') # plotting the spectrum
plt.xlabel('Freq (Hz)')
plt.ylabel('Spectrum')
plt.title('Spectrum of the signal')
plt.show()
