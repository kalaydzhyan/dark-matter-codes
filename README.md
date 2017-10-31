# Dark matter tests

The project deals with dark matter detection methods with space detectors. More codes will be added with time.

## curve.py + main.py

This code generates a stochastic background of scalar field waves and analyzes a response from two-detectors (two pairs of atomic clocks). Main functions and classes are defined in **main.py**.
The main result are the angular curves for cross-spectra, see preliminary notes in **stochastic.pdf**. 

Parameters controlling the output (with default values):
```
short_wavelength_regime = True
```
This parameter switches between the short- and long-wavelength regimes (see the notes), i.e., the scalar wave wavelength is shorter or longer that the distance between detector in a pair.

```
delta_zeta = 10
number_sources = 1000
observation_time = 1000
sampling_time = 1
```
Angle step (in degrees) for the calculation of the angular function, number of wave sources and observation/sampling time (in seconds). These parameters change the precision of the simulation as well as the performance of the code.

## noise.py + main1.py

Similar code, to test the signal-to-noise ratio dependence on the observation time. File **main1.py** added temporarily, to use a flat power spectrum of scalar waves, instead of quasi-monochromatic in **main.py**.