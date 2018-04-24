from __future__ import division
from matplotlib import pyplot as plt
import numpy as np
import math, sys

""" Some global variables """
mass = 1 # mass of the scalar field, add it later

""" Properties of the source frequencies """
frequency_mean = 0.1
frequency_variation = 1


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
            
        amplitude   = 10 + np.random.random()
        """ Flat spectrum: """
        frequency   = (np.random.random()-1/2)*frequency_variation
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
            
    def noise(self, t):
        
        amplitude=0.05
        mean = 0
        std = 1

        return amplitude*np.random.normal(mean, std)

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

# Print iterations progress
def printProgressBar (iteration, total, prefix = '', suffix = '', decimals = 1, length = 100, fill = 'X'):
    """
    Call in a loop to create terminal progress bar
    @params:
        iteration   - Required  : current iteration (Int)
        total       - Required  : total iterations (Int)
        prefix      - Optional  : prefix string (Str)
        suffix      - Optional  : suffix string (Str)
        decimals    - Optional  : positive number of decimals in percent complete (Int)
        length      - Optional  : character length of bar (Int)
        fill        - Optional  : bar fill character (Str)
    """
    percent = ("{0:." + str(decimals) + "f}").format(100 * (iteration / float(total)))
    filledLength = int(length * iteration // total)
    bar = fill * filledLength + '-' * (length - filledLength)
    print('\r%s |%s| %s%% %s' % (prefix, bar, percent, suffix))
    # Print New Line on Complete
    if iteration == total: 
        print '\n'

def distance(X, Y):
    XY = np.array(X) - np.array(Y)
    return np.sqrt(np.inner(XY, XY))

if __name__ == '__main__':
    print "\n Test run of the library.\n"

    number_sources = 1000
    observation_time = 200
    sampling_time = 1

    distance1 = 0.001
    distance2 = 0.0001
    delta_zeta = 10  

    timeline = [t for t in np.arange(0, observation_time, sampling_time)]
    
    sources = Source()
    for source_number in range(number_sources): 
        sources.add()

    detector0 = Detector()
    signal0 = [sources.signal(t, detector0.position) \
        for t in np.arange(0, observation_time, sampling_time)]

    ratios = []

    zetas = np.arange(0, math.pi, np.deg2rad(delta_zeta))
    timeline = np.arange(0, observation_time, sampling_time)

    for zeta in zetas:
        detector1 = Detector([distance1, 0, 0])
        signal1 = [sources.signal(t, detector1.position) + detector1.noise(t) \
                       for t in timeline]

        detector2 = Detector([distance2 * np.cos(zeta), distance2 * np.sin(zeta), 0])
        signal2 = [sources.signal(t, detector2.position) + detector2.noise(t) \
                       for t in timeline]

        signal1 = np.array(signal1) - np.array(signal0)
        signal2 = np.array(signal2) - np.array(signal0)
        
        frequencies, signal1image = zip(*Fourier(zip(timeline, signal1), sampling_time))
        frequencies, signal2image = zip(*Fourier(zip(timeline, signal2), sampling_time))
    
        x_fourier_image = signal1image * np.conj(signal2image)
    
        if frequency_variation*observation_time < 1:
            print 'Frequency variation must be larger than observation time, cannot resolve!'
    
        #left_f = int((frequency_mean - frequency_variation) * observation_time)
        #right_f = int((frequency_mean + frequency_variation) * observation_time)
        
        left_f = 1
        right_f = observation_time - 1
        
        cross_spectrum_f = np.real(np.array(x_fourier_image))
        S_cross = np.mean(cross_spectrum_f[left_f:right_f])
        detector1_spectrum_f = np.abs(np.array(signal1image))**2
        S_d1 = np.mean(detector1_spectrum_f[left_f:right_f])

        ratios.append(S_cross/S_d1)
    
    fig, ay = plt.subplots()
    ay.plot(zetas, ratios,'go', zetas, (distance2/distance1)*np.cos(zetas),'g-')

    plt.xlabel('$\zeta$, rad')
    plt.title('Angular curve')
    plt.savefig("angular_long_noise.png")

    plt.show()

    """
    noise_record = [detector0.noise(t) for t in timeline]
    frequencies, noise_image = zip(*Fourier(zip(timeline, noise_record), sampling_time))

    plt.plot(frequencies, np.abs(noise_image),'r-') 
    plt.xlabel('f (Hz)')
    plt.ylabel('$\sqrt{S_X(f)}$')
    plt.title('Power spectrum of the noise for the reference detector')
    plt.show()
    """