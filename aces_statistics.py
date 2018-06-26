from __future__ import division
from matplotlib import pyplot as plt
import numpy as np
import pandas as pd
import math, sys, glob, re, os
import scipy.constants as phco


""" geocentric gravitational constant (m^3/s^2) """
GM = 3.986004418e14
""" Earth radius (m) """
Rearth = 6378000
""" Earth rotation (rad/s) """
wE = 7.292e-5

current_path = os.path.dirname(os.path.abspath(__file__))

""" Collecting delays """

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

""" Finding position of the object by the coordinate time """
def position_by_time(dataframe, time):
    
    if (dataframe['coord.time'].iloc[0] < time < dataframe['coord.time'].iloc[-1]):
        array = dataframe['coord.time']
        idx = np.searchsorted(array, time, side="left")[0]
        delta_t = time - dataframe['coord.time'][idx-1]
        """ Interpolating, assuming 1 sec sampling """ 
        object_position = position(dataframe, idx-1) + \
            (position(dataframe, idx)-position(dataframe, idx-1))*delta_t
    else:
        print "Cannot find the specified coordinate time!"
    
    return object_position

""" Finding velocity of the object by the coordinate time """
def velocity_by_time(dataframe, time):
    
    if (dataframe['coord.time'].iloc[0] < time < dataframe['coord.time'].iloc[-1]):
        array = dataframe['coord.time']
        idx = np.searchsorted(array, time, side="left")[0]
        delta_t = time - dataframe['coord.time'][idx-1]
        """ Interpolating, assuming 1 sec sampling """ 
        object_velocity = velocity(dataframe, idx-1) + \
            (velocity(dataframe, idx)-velocity(dataframe, idx-1))*delta_t
    else:
        print "Cannot find the specified coordinate time!"
    
    return object_velocity

def collect_data(dataset_number = 0, carrier_code = "ca", frequency_n = "1"):
    
    directories = glob.glob(current_path + '/v4.3.2_mb_53896_53907/gs999/*')

    if (dataset_number >= len(directories)):
        print "Dataset does not exist! Change dataset_number. \n"

    """ Extracting the m. Julian day for the dataset"""
    mjd = re.search('%s(.*)%s' % ('mjd','\.'), directories[dataset_number]).group(1)

    path = directories[dataset_number] + '/theo/'
    #path = current_path[0] + '/v4.3.2_mb_53896_53907/gs999/mjd53897.220682_ch1/theo/'
    
    filename = 'v4.3.2_theo_f' + frequency_n + '_' + carrier_code + '.dat'
    data = pd.read_csv(path + filename, delim_whitespace = True, skiprows=1, \
                   names=['time.tag','MJD', 'coord.time', 'prop.time','PToF','geometric.ToF','shapiro','tropo.delay',\
                          'iono.1/f^2', 'iono.1/f^3', 'STEC', 'desynchronisation', 'rangeX', 'rangeY', 'rangeZ', \
                          'sequence.number','Tm'])

    data['distance'] = np.sqrt(data['rangeX']**2 + data['rangeY']**2 + data['rangeZ']**2)

    
    """ Subtracting delays from PToF """
    data['all.delays'] = data['geometric.ToF'] + data['tropo.delay'] + data['iono.1/f^2'] + data['iono.1/f^3']
    data['test'] = data['all.delays'] + data['PToF'] + 1e-3*(-1 if frequency_n=="1" else 1)
    """ This additional +/- 1ms doesn't really matter, because we'll take a derivative"""
    
      
    return data

    #dataset_number = 0
    #carrier_code = "ca"      # ["ca", "co"]
    #frequency_n = "3"        # ["1", "2", "3"], where "1": GS -> ISS; "2" and "3": ISS -> GS

def collect_trajectories(dataset_number = 0):
    
    directories = glob.glob(current_path + '/v4.3.2_mb_53896_53907/gs999/*')

    if (dataset_number >= len(directories)):
        print "Dataset does not exist! Change dataset_number. \n"

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

    return gs_orbit, iss_orbit
    
    
""" one-way relativistic effect, theory and experiment."""
""" Theory, as defined by SYRTE, needs to be corrected later! """

def relativistic_effect(dataframe, frequency_n = '1'):
    
    r_A = dataframe["iss_positions"]
    r_B = dataframe["gs_positions"]
    
    abs_r_A = np.array(map(lambda x: np.linalg.norm(x), r_A))
    abs_r_B = np.array(map(lambda x: np.linalg.norm(x), r_B))
    R_AB = r_B - r_A
    abs_R_AB = np.array(map(lambda x: np.linalg.norm(x), R_AB))
    N_AB = R_AB/abs_R_AB
    
    v_A = dataframe["iss_velocities"]
    v_B = dataframe["gs_velocities"]
    abs_v_A = np.array(map(lambda x: np.linalg.norm(x), v_A))
    abs_v_B = np.array(map(lambda x: np.linalg.norm(x), v_B))
    
    potential_difference = GM*(1/abs_r_B - 1/abs_r_A)/(phco.c**2)
    second_doppler = .5*(abs_v_B**2 - abs_v_A**2)/(phco.c)**2
    
    dataframe["rel.effect_theory"] = (potential_difference + second_doppler)*(-1 if frequency_n=="1" else 1)
    
    time_steps_experiment = np.gradient(np.array(dataframe['coord.time']))
    data_derivatives = np.gradient(np.array(dataframe['test']))/time_steps_experiment
    
    dataframe["rel.effect_experiment"] = data_derivatives
    dataframe["difference"] = dataframe["rel.effect_experiment"] - dataframe["rel.effect_theory"]
    
    return dataframe

if __name__ == '__main__':
    print "\n Test run of the library.\n"
    
    """ Collecting 2-way data """
    data1 = collect_data(dataset_number = 0, carrier_code = "ca", frequency_n = "1")
    data2 = collect_data(dataset_number = 0, carrier_code = "ca", frequency_n = "2")
    gs_orbit, iss_orbit = collect_trajectories(dataset_number = 0)
    
    data1["gs_positions"] = map(lambda x: position_by_time(gs_orbit, x), data1['coord.time'])
    data1["gs_velocities"] = map(lambda x: velocity_by_time(gs_orbit, x), data1['coord.time'])
    data1["iss_positions"] = map(lambda x: position_by_time(iss_orbit, x), data1['coord.time'])
    data1["iss_velocities"] = map(lambda x: velocity_by_time(iss_orbit, x), data1['coord.time'])
    
    data2["gs_positions"] = map(lambda x: position_by_time(gs_orbit, x), data2['coord.time'])
    data2["gs_velocities"] = map(lambda x: velocity_by_time(gs_orbit, x), data2['coord.time'])
    data2["iss_positions"] = map(lambda x: position_by_time(iss_orbit, x), data2['coord.time'])
    data2["iss_velocities"] = map(lambda x: velocity_by_time(iss_orbit, x), data2['coord.time'])
    
    """ Check if these are correct datasets """
    assert np.absolute(data1["coord.time"][0] - data2["coord.time"][0]) < 2e-3
    

    """ theoretical and experimental residuals """
    
    data1 = relativistic_effect(data1, frequency_n = '1')
    data2 = relativistic_effect(data2, frequency_n = '2')
    
    data1.plot(x = 'coord.time', y = ["rel.effect_experiment", "rel.effect_theory"])
    data1.plot(x = 'coord.time', y = ["difference"])
    
    data2.plot(x = 'coord.time', y = ["rel.effect_experiment", "rel.effect_theory"])
    data2.plot(x = 'coord.time', y = ["difference"])
    
    print "uplink difference, theory vs. experiment: ", np.mean(data1['difference']), ' +/- ', np.std(data1['difference'])
    print "downlink difference, theory vs. experiment: ", np.mean(data2['difference']), ' +/- ', np.std(data2['difference'])