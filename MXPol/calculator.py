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

def calc_Ps_from_ref(atoms, ref_atoms, calc, nsteps = 5):
    '''calculates spontaneous polarization according to ref_atoms
    which is assumed to have zero polarization which is accurate 
    to the first order. Prepare ref atoms using 
    MXPol.build.supercell_from_ref.
    Equations are from Modern Theory of Polarization 
    by Resta, R., & Vanderbilt, D. (eq 13) and 
    10.1088/0953-8984/22/16/165902 (eqs 4,5)'''
    from ase.mep.neb import interpolate
    from becqsdr.model import E3NN
    from ase.calculators.vasp import Vasp
    #ref_atoms = supercell_from_ref(atoms)
    #write('ref2.vasp',ref_atoms, format='vasp',sort=True)
    lin_path = [ref_atoms.copy()]*(nsteps-1) + [atoms.copy()]
    interpolate(lin_path, mic=True)
    if isinstance(calc, E3NN):
        bec_path = [calculate_bec(a, calc) for a in lin_path]
    elif isinstance(calc, Vasp):
        from becqsdr.io import read_vasp_bec
        bec_path = []
        for a in lin_path:
            a.calc = calc
            a.get_potential_energy()
            xml_file = f'{calc.directory}/vasprun.xml'
            arr = np.loadtxt(f'{calc.directory}/ase-sort.dat',dtype=int)
            mapping = {k:i for i, (k,v) in enumerate(arr)}
            bec_path.append(read_vasp_bec(xml_file)[[mapping[i] for i in range(len(atoms))]])
    #Ps = np.average(bec_path, axis=0) @ (atoms_unit_cell.get_positions() - ref_atoms.get_positions())
    Ps = []
    bec_avg = np.average(bec_path, axis=0)
    d_pos = atoms.get_positions() - ref_atoms.get_positions()
    for p, b in zip(d_pos, bec_avg):
        #b = b.numpy()
        #p-= ref_pos
        Ps.append(b @ p)
    Ps = np.array(Ps)
    print(Ps.shape)
    Ps = Ps.reshape(-1,4,3)
    return np.average(Ps,axis=1) #,d_pos

def calc_Ps(atoms, calc, a_length, b_length, nsteps = 5):
    '''wrapper for calc_Ps_from_ref.
    assumes that the atoms are sorted according to 
    MXPol stride (a then b in a group of 4)'''
    from MXPol.build import supercell_from_ref
    ref_atoms = supercell_from_ref(atoms, a_length, b_length, sorted=False)
    #write('ref2.vasp',ref_atoms, format='vasp',sort=True)
    return calc_Ps_from_ref(atoms, ref_atoms, calc, nsteps)