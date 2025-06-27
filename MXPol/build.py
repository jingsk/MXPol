import numpy as np
from monochalcogenpy.build import unit_cell

def reflect_Ge(atoms, axis):
    '''
    Reflect the Ge atoms by Se anchor points
    Currently hardcoded atomic species Ge to be 0,1 and Se to be 2,3.
    TODO: Needs automatic detection
    '''

    pos = atoms.get_positions()
    #overall_shift = pos[2,axis] - pos[0,axis]
    pos[0,axis] = 2 * pos[2,axis] - pos[0,axis]
    pos[1,axis] = 2 * pos[3,axis] - pos[1,axis]
    #pos += overall_shift
    atoms.set_positions(pos)
    return atoms

def init_unit_cell(a,b,c, relaxed = False, **db_kwargs):
    '''
    return unit cell of specified lattice vectors from either 
    1. monochalcogenpy.build (relaxed=False)
    2. specified db (relaxed=True)
    '''
    if relaxed:
        from ase.db import connect
        path = db_kwargs['path'] 
        db_read = connect(path + db_kwargs['db_name'])
        atoms = db_read.get(a=a,b=b,relaxed=True, model='macecalculator2').toatoms()
        atoms.set_constraint()
        return atoms
    else:
        orientation = 'ac' if a>b else 'bc' 
        return unit_cell(a,b,c, orientation = orientation, use_symm=False)

def atoms_grid(a_grid, b_grid, c, reflect, relaxed = False, **db_kwargs):
    '''
    return a grid of atoms with specified a_grid, b_grid, c.
    Reflection reflects Se at midway point.
    '''
    atoms_array = []
    #n_grid = len(a_grid)
    for j, b in enumerate(b_grid):
        atoms_list = []
        for i, a in enumerate(a_grid):
            #atoms = db_read.get(a=a,b=b,relaxed=True, model='macecalculator2').toatoms()
            atoms = init_unit_cell(a,b,c, relaxed = relaxed, **db_kwargs) #.set_constraint()
            #this is to make sure for an even grid size, the number of 
            # opposing polarization are equal. 
            if (reflect[i,j]) and (a != b):
                flip_dir = 0 if a>b else 1
                atoms = reflect_Ge(atoms.copy(), flip_dir)
            atoms_list.append(atoms)
        atoms_array.append(atoms_list)
        
    wire_list = []
    for atoms_list in atoms_array:
        nanowire = atoms_list[0].copy()
        for atoms in atoms_list[1:]:
            this_atoms = atoms.copy()
            offset = nanowire.get_cell()[0,0]
            this_atoms.translate([offset,0,0])
            cell = this_atoms.get_cell()[:]
            cell[0,0] += offset
            nanowire.set_cell(cell,scale_atoms = False)
            nanowire += this_atoms
        wire_list.append(nanowire)
    return wire_list

def construct_atoms_from_grid(atoms_grid):
    '''
    given an atoms grid make a supercell
    '''
    nanosheet = atoms_grid[0].copy()
    for nanowire in atoms_grid[1:]:
        this_wire = nanowire.copy()
        offset = nanosheet.get_cell()[1,1]
        this_wire.translate([0,offset,0])
        cell = this_wire.get_cell()[:]
        cell[1,1] += offset
        nanosheet.set_cell(cell,scale_atoms = False)
        nanosheet += this_wire
    return nanosheet 

def supercell(a_grid, b_grid, c, reflect, relaxed = False, **db_kwargs):
    '''
    user will call this function. Give a_grid, b_grid, c to make 
    supercell.
    '''
    atoms_grid = atoms_grid(a_grid, b_grid, c, reflect, relaxed = relaxed, **db_kwargs)
    return construct_atoms_from_grid(atoms_grid)

# def supercell_from_ref(ref_atoms, relaxed = False, **db_kwargs):
#     #atoms_grid = atoms_grid(a_grid, b_grid, c, reflect, relaxed = relaxed, **db_kwargs)
#     return construct_atoms_from_grid(atoms_grid)