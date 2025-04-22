import numpy as np
from io import config

def wave_grid(n, A, waveform = 'smooth square'):
    wl = n
    x = np.arange(n)
    if waveform =='triangular':    
        grid = A *(8 / np.pi**2) * (
            np.sin(2  * np.pi * x/wl) / 1**2 - 
            np.sin(6  * np.pi * x/wl) / 3**2 + 
            np.sin(10 * np.pi * x/wl) / 5**2
        )
    if waveform =='smooth square':
        d = wl /200
        grid = A * np.sin(2  * np.pi * x/wl) / (np.sqrt(d**2 + np.sin(2  * np.pi * x/wl)**2))
    #atoms_flexed.translate(dis_vec)
    if waveform == 'uniform':
        grid = np.zeros(n)
    return grid

build_config = config['build']
a_grid = build_config['b'] + \
    np.abs(wave_grid(
        build_config['n_grid'], 
        build_config['wave_amp'], 
        build_config['waveform']
        )
    )