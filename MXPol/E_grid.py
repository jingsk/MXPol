import numpy as np
from scipy import signal
from io import config

dyn_config = config['dynamics']
E_waveform = dyn_config['E_waveform']
E_wavelength = dyn_config['E_wavelength_fs']
E_wave_amp = dyn_config['E_wave_amp_mV_per_ang']
tot_step = dyn_config['tot_step'] 

if dyn_config['E_waveform'] == 'uniform':
    E_grid = E_wave_amp * np.ones(tot_step+ 1)

if dyn_config['E_waveform'] == 'square':
    step_grid = np.arange(tot_step + 1)
    E_grid = E_wave_amp * signal.square(2 * np.pi * step_grid / E_wavelength)

if dyn_config['E_waveform'] == 'sine':
    step_grid = np.arange(tot_step + 1)
    E_grid = E_wave_amp * np.sin(2 * np.pi * step_grid / E_wavelength)

