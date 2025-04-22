import numpy as np
from scipy import signal
from io import config

dyn_config = config['dynamics']
E_waveform = dyn_config['E_waveform']
E_wavelength = dyn_config['E_wavelength_fs']
tot_step = dyn_config['tot_step'] 

def E_grid(E_wave_amp):
    if dyn_config['E_waveform'] == 'uniform':
        grid = E_wave_amp * np.ones(tot_step+ 1)

    if dyn_config['E_waveform'] == 'square':
        step_grid = np.arange(tot_step + 1)
        grid = E_wave_amp * signal.square(2 * np.pi * step_grid / E_wavelength)

    if dyn_config['E_waveform'] == 'sine':
        step_grid = np.arange(tot_step + 1)
        grid = E_wave_amp * np.sin(2 * np.pi * step_grid / E_wavelength)
    else:
        raise NotImplementedError('waveform not yet implemented')
    return grid

