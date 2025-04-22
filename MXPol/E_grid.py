import numpy as np
from scipy import signal
from MXPol.io import config
import warnings

#dyn_config = config['dynamics']
# E_waveform = dyn_config['E_waveform']
# E_wavelength = dyn_config['E_wavelength_fs']
# tot_step = dyn_config['tot_step'] 

def E_grid(
        E_wave_amp,
        E_waveform,
        E_period,
        time
        ):
    if E_waveform == 'uniform':
        Ex_grid = E_wave_amp * np.ones_like(time)

    elif E_waveform == 'square':
        if time[time<E_period].size % 2 == 1:
            warnings.warn("Grid size in each period is an odd number. Expect more positive field than negative.")
        Ex_grid = E_wave_amp * signal.square(2 * np.pi * time / E_period)

    elif E_waveform == 'sine':
        Ex_grid = E_wave_amp * np.sin(2 * np.pi * time / E_period)
        
    else:
        raise NotImplementedError('waveform not yet implemented')
    E_field = np.column_stack([
    Ex_grid, 
    np.zeros_like(Ex_grid), 
    np.zeros_like(Ex_grid)
    ])
    return E_field



