from __future__ import division
from matplotlib import pyplot as plt
import numpy as np
import pandas as pd
import sys, glob, re
import scipy.constants as phco


""" geocentric gravitational constant (m^3/s^2) """
GM = 3.986004418e14
""" Earth radius (m) """
Rearth = 6378000
""" Earth rotation (rad/s) """
wE = 7.292e-5

def position(dataframe, element=0):
    return np.array([dataframe['x'][element], dataframe['y'][element], dataframe['z'][element]])

def velocity(dataframe, element=0):
    return np.array([dataframe['Vx'][element], dataframe['Vy'][element], dataframe['Vz'][element]])

def unit_vector(vector):
    try:
        normalized = vector/np.linalg.norm(vector)
    except ZeroDivisionError:
        sys.exit("Can't normalize, zero vector!")
    return normalized

""" One-way relativistic effect without 1st order Doppler, precision 1/c^3 """

def relativistic_effect(i, lower_precision = False):
    assert gs_orbit.shape == iss_orbit.shape
    
    r_A = position(iss_orbit, i)
    r_B = position(gs_orbit, i)
    
    abs_r_A = np.linalg.norm(r_A)
    abs_r_B = np.linalg.norm(r_B)
    R_AB = r_B - r_A
    abs_R_AB = np.linalg.norm(R_AB)
    N_AB = R_AB/abs_R_AB
    
    v_A = velocity(iss_orbit, i)
    v_B = velocity(gs_orbit, i)
    abs_v_A = np.linalg.norm(v_A)
    abs_v_B = np.linalg.norm(v_B)
     
    common_factor = (4*GM/(phco.c**3))/((abs_r_A + abs_r_B)**2 - (abs_R_AB)**2)
    
    """ Here one can lower precision to 1/c^2, if needed """
    if lower_precision:
        q_A = q_B = 1
    else:
        q_A = 1 - common_factor*((abs_r_A + abs_r_B)*np.dot(N_AB, v_A) + abs_R_AB*np.dot(r_A, v_A)/abs_r_A)
        q_B = 1 - common_factor*((abs_r_A + abs_r_B)*np.dot(N_AB, v_B) - abs_R_AB*np.dot(r_B, v_B)/abs_r_B)

    result = 1 - (q_A/q_B)*(1 - (GM/abs_r_B + 0.5*abs_v_B**2)/(phco.c**2))/(1 - (GM/abs_r_A + 0.5*abs_v_A**2)/(phco.c**2))
    
    return result*(-1 if frequency_n=="1" else 1)

current_path = '/Users/tigrank/Tigran/research/ACES'
directories = glob.glob(current_path + '/v4.3.2_mb_53896_53907/gs999/*')

""" Collecting delays """

for dataset_number in range(len(directories)):

    print 'Calculating dataset number', dataset_number
    
    """ Extracting the m. Julian day for the dataset"""
    mjd = re.search('%s(.*)%s' % ('mjd','\.'), directories[dataset_number]).group(1)

    """ Collecting ISS positions and velocities """

    path = current_path + '/v4.3.2_mb_53896_53907/auxdata/iss/orbit/'
    filename = 'v4.3.2_orb_KU_' + mjd
    iss_orbit = pd.read_csv(path + filename, delim_whitespace = True, skiprows=1, \
                   names=['julian.day','coord.time', 'x', 'y','z','Vx','Vy','Vz'])

    """ Collecting ground station positions and velocities """

    path = current_path + '/v4.3.2_mb_53896_53907/auxdata/gs999/orbit/'
    filename = 'v4.3.2_orbGS_'  + mjd
    gs_orbit = pd.read_csv(path + filename, delim_whitespace = True, skiprows=1, \
                   names=['julian.day','coord.time', 'x', 'y','z','Vx','Vy','Vz'])

    path = directories[dataset_number] + '/theo/'

    """ Iterating over frequencies and carrier/code """
    
    for carrier_code in ["ca", "co"]:
        for frequency_n in ["1", "2", "3"]:
            filename = 'v4.3.2_theo_f' + frequency_n + '_' + carrier_code + '.dat'
            
            data = pd.read_csv(path + filename, delim_whitespace = True, skiprows=1, \
                   names=['time.tag','MJD', 'coord.time', 'prop.time','PToF','geometric.ToF','shapiro','tropo.delay',\
                          'iono.1/f^2', 'iono.1/f^3', 'STEC', 'desynchronisation', 'rangeX', 'rangeY', 'rangeZ', \
                          'sequence.number','Tm'])

            data['distance'] = np.sqrt(data['rangeX']**2 + data['rangeY']**2 + data['rangeZ']**2)

            data['test'] = data['geometric.ToF'] + data['tropo.delay'] + data['iono.1/f^2'] + data['iono.1/f^3'] \
               + data['PToF'] + 1e-3*(-1 if frequency_n=="1" else 1)
            """ This additional +/- 1ms doesn't really matter, because we'll take a derivative"""
    
            """ Taking derivative of that data """

            time_steps_experiment = np.gradient(np.array(data['coord.time']))
            data_derivatives = np.gradient(np.array(data['test']))/time_steps_experiment

            relativistic_array_c_cube = [relativistic_effect(i) for i in range(gs_orbit.shape[0])]
            relativistic_array_c_square = [relativistic_effect(i, lower_precision=True) for i in range(gs_orbit.shape[0])]

            """ Assuming 1 sec time steps """

            t_first = int(data['coord.time'][0])
            t_last = int(data['coord.time'][-1])

            plt.plot(np.array(gs_orbit['coord.time'][t_first:t_last]), relativistic_array_c_cube[t_first:t_last], c = 'k')
            plt.plot(np.array(gs_orbit['coord.time'][t_first:t_last]), relativistic_array_c_square[t_first:t_last], c = 'r')
            plt.scatter(data['coord.time'], data_derivatives, c = 'orange', marker = '+')
            plt.ylim(np.min(data_derivatives) - 1e-13, np.max(data_derivatives) + 1e-13)
            plt.legend(["Theory 1/c^3", "Theory 1/c^2", "Simulated data"])
            plt.title('One-way desynchronization, frequency #' + frequency_n + ',\n' + ("carrier" if carrier_code=="ca" else "code") + ', 1st order Doppler-subtracted')
            plt.xlabel('coordinate time (s)')
            plt.savefig('/Users/tigrank/Dropbox/output_for_makan/plots/desync_mjd' + str(mjd) + '_f' + frequency_n + '_' + carrier_code + '_dataset' + str(dataset_number)+'.png')
            plt.close()