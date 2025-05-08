import numpy as np
from monochalcogenpy.build import unit_cell

def reflect_Ge(atoms):
    '''
    Reflect the Ge atoms by Se anchor points
    Currently hardcoded atomic species Ge to be 0,1 and Se to be 2,3.
    TODO: Needs automatic detection
    '''

    pos = atoms.get_positions()
    pos[0,0] = 2 * pos[2,0] - pos[0,0]
    pos[1,0] = 2 * pos[3,0] - pos[1,0]
    atoms.set_positions(pos)
    return atoms

def supercell(a_grid, b, c, reflect):
    '''
    make 1d supercell with specified a_grid, b, c.
    Reflection reflects Se at midway point.
    '''
    n_grid = len(a_grid)
    atoms_list = []
    for i, a in enumerate(a_grid):
        atoms = unit_cell(a,b,c, orientation = 'ac', use_symm=False)
        #this is to make sure for an even grid size, the number of 
        # opposing polarization are equal. 
        if isinstance(reflect, bool):
            if i > n_grid//2-1:
                #atoms = reflect(atoms.copy(), [1.0,0,0], center = [a/2, b/2, c/2])
                atoms = reflect_Ge(atoms.copy()) if reflect else atoms.copy()
        elif isinstance(reflect, list):
            atoms = reflect_Ge(atoms.copy()) if reflect[i] else atoms.copy()
        atoms_list.append(atoms)
    
    nanosheet = atoms_list[0].copy()
    for atoms in atoms_list[1:]:
        this_atoms = atoms.copy()
        offset = nanosheet.get_cell()[0,0]
        this_atoms.translate([offset,0,0])
        cell = this_atoms.get_cell()[:]
        cell[0,0] += offset
        nanosheet.set_cell(cell,scale_atoms = False)
        nanosheet += this_atoms
    return nanosheet
