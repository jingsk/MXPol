from ase.calculators.calculator import Calculator, all_changes
from becqsdr.calculator import calculate_bec
import numpy as np

class bec_under_field(Calculator):
    """
    calculate forces and due to applied E-field
    """
    
    implemented_properties = ['energy','charges', 'forces']
    #implemented_properties += ['stress', 'stresses']  # bulk properties
    default_parameters = {
        'E_field': 0 , # common unit of V/m, q* eV/A in ase unit
    }
    nolabel = True

    def __init__(self, **kwargs):
      """
      Parameters
      ----------
      sigma: E_field
        A vector of E field along x,y,z 
        Applied electric field in V/A (ase units)
      """
    
      Calculator.__init__(self, **kwargs)
    
    def set_E_field(self, E_field):
        self.parameters.E_field = E_field

    def calculate(
      self,
      atoms=None,
      properties=None,
      system_changes=all_changes,
    ):
        if properties is None:
            properties = self.implemented_properties
        
        Calculator.calculate(self, atoms, properties, system_changes)
        
        natoms = len(self.atoms)
        if isinstance(self.parameters.E_field, int):
            E_field = self.parameters.E_field * np.ones(3)
        elif isinstance(self.parameters.E_field, list):
            assert len(self.parameters.E_field) == 3
            E_field = np.array(self.parameters.E_field)
        elif isinstance(self.parameters.E_field, np.ndarray):
            E_field = self.parameters.E_field
        else:
            raise ValueError('"E_field" must be a three-component vector.')
        enn = self.parameters.enn
        
        energies = np.zeros(natoms)
        forces = np.zeros((natoms, 3))
        #charges = np.zeros(natoms)
        stresses = np.zeros((natoms, 3, 3))
        atoms.info['bec'] = calculate_bec(atoms, enn)
        bec = atoms.info['bec'].numpy()
        
        for i in range(natoms):
            forces[i] += bec[i] @ E_field
        
        #energy = energies.sum()
        self.results['energy'] = 0
        #self.results['charges'] = charges
        #self.results['energies'] = energies
        self.results['free_energy'] = 0
        self.results['forces'] = forces