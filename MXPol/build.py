import numpy as np
from monochalcogenpy.build import unit_cell
from itertools import product

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

def init_unit_cell(a,b,c, mode = 'monochalcogenpy', **db_kwargs):
    '''
    return unit cell of specified lattice vectors from either 
    1. monochalcogenpy.build (mode = monochalcogenpy)
    2. specified db (mode = relaxed)
    3. monochalcogenpy.build but zero polarization (mode = reference)
    '''
    if mode == 'monochalcogenpy':
        orientation = 'ac' if a>b else 'bc' 
        return unit_cell(a,b,c, orientation = orientation, use_symm=False)
    elif mode == 'relaxed':
        from ase.db import connect
        path = db_kwargs['path'] 
        db_read = connect(path + db_kwargs['db_name'])
        atoms = db_read.get(a=a,b=b,relaxed=True, model='macecalculator2').toatoms()
        atoms.set_constraint()
        return atoms
    elif mode == 'reference':
        ref_atoms = unit_cell(b,b,c)
        ref_atoms.set_cell([a, b, c],scale_atoms=True)
        return ref_atoms

# def make_atoms_grid(a_grid, b_grid, c, reflect, mode, **db_kwargs):
#     '''
#     return a grid of atoms with specified a_grid, b_grid, c.
#     Reflection reflects Se at midway point.
#     '''
#     atoms_array = []
#     if isinstance(reflect, bool):
#         if reflect:
#             reflect = np.ones([len(a_grid),len(b_grid)], dtype=bool)
#         else:
#             reflect = np.zeros([len(a_grid),len(b_grid)], dtype=bool)
#     #n_grid = len(a_grid)
#     for j, b in enumerate(b_grid):
#         atoms_list = []
#         for i, a in enumerate(a_grid):
#             #atoms = db_read.get(a=a,b=b,relaxed=True, model='macecalculator2').toatoms()
#             atoms = init_unit_cell(a,b,c, mode = mode, **db_kwargs) #.set_constraint()
#             #this is to make sure for an even grid size, the number of 
#             # opposing polarization are equal. 
#             if (reflect[i,j]) and (a != b):
#                 flip_dir = 0 if a>b else 1
#                 atoms = reflect_Ge(atoms.copy(), flip_dir)
#             atoms_list.append(atoms)
#         atoms_array.append(atoms_list)
        
#     wire_list = []
#     for atoms_list in atoms_array:
#         nanowire = atoms_list[0].copy()
#         for atoms in atoms_list[1:]:
#             this_atoms = atoms.copy()
#             offset = nanowire.get_cell()[0,0]
#             this_atoms.translate([offset,0,0])
#             cell = this_atoms.get_cell()[:]
#             cell[0,0] += offset
#             nanowire.set_cell(cell,scale_atoms = False)
#             nanowire += this_atoms
#         wire_list.append(nanowire)
#     return wire_list


def make_atoms_grid(ab_grid, c, reflect, mode, **db_kwargs):
    '''
    return a grid of atoms with specified a_grid, b_grid, c.
    Reflection reflects Se at midway point.
    '''
    atoms_array = np.empty([ab_grid.shape[0], ab_grid.shape[1]], dtype=object)
    if isinstance(reflect, bool):
        if reflect:
            reflect = np.ones([ab_grid.shape[0],ab_grid.shape[1]], dtype=bool)
        else:
            reflect = np.zeros([ab_grid.shape[0],ab_grid.shape[1]], dtype=bool)
    #n_grid = len(a_grid)
    for i,j in product(range(ab_grid.shape[0]), range(ab_grid.shape[1])):
        a,b = ab_grid[i,j]
        
        #atoms = db_read.get(a=a,b=b,relaxed=True, model='macecalculator2').toatoms()
        atoms = init_unit_cell(a,b,c, mode = mode, **db_kwargs) #.set_constraint()
        #this is to make sure for an even grid size, the number of 
        # opposing polarization are equal. 
        if (reflect[i,j]) and (a != b):
            flip_dir = 0 if a>b else 1
            atoms = reflect_Ge(atoms.copy(), flip_dir)
        atoms_array[i,j] = atoms.copy()
    return atoms_array.tolist()
        
    # wire_list = []
    # for atoms_list in atoms_array:
    #     nanowire = atoms_list[0].copy()
    #     for atoms in atoms_list[1:]:
    #         this_atoms = atoms.copy()
    #         offset = nanowire.get_cell()[0,0]
    #         this_atoms.translate([offset,0,0])
    #         cell = this_atoms.get_cell()[:]
    #         cell[0,0] += offset
    #         nanowire.set_cell(cell,scale_atoms = False)
    #         nanowire += this_atoms
    #     wire_list.append(nanowire)
    return wire_list

def construct_atoms_from_grid(atoms_grid):
    '''
    given an atoms grid make a supercell
    '''
    wire_list = []
    for atoms_list in atoms_grid:
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
    nanosheet = wire_list[0].copy()
    for nanowire in wire_list[1:]:
        this_wire = nanowire.copy()
        offset = nanosheet.get_cell()[1,1]
        this_wire.translate([0,offset,0])
        cell = this_wire.get_cell()[:]
        cell[1,1] += offset
        nanosheet.set_cell(cell,scale_atoms = False)
        nanosheet += this_wire
    return nanosheet 

def supercell(a_grid, b_grid, c, reflect, mode = 'monochalcogenpy', **db_kwargs):
    '''
    user will call this function. Give a_grid, b_grid, c to make 
    supercell.
    '''
    ab_grid = np.array(np.meshgrid(a_grid, b_grid)).transpose([1,2,0])
    atoms_grid = make_atoms_grid(ab_grid, c, reflect, mode = mode, **db_kwargs)
    return construct_atoms_from_grid(atoms_grid)

def supercell_from_ref(ref_atoms, a_length, b_length, sorted=False, **sort_kwargs):
    '''
    Make a zero polarization structure for using as a reference 
    to calculate sponteneous polarization. 
    Assuming that stretching a=b structure makes a zero polarization 
    unit cell, which is True up to first order approximation.
    '''
    from MXPol.utils import atoms_mapping, element_index_grid
    if sorted == True:
        mapping = atoms_mapping(atoms_key = sort_kwargs['unsorted_atoms'], 
                                atoms_val = ref_atoms)
    else:
        mapping = {}

    grid = element_index_grid(a_length, b_length, 
                              offset = 2, mapping=mapping)
    
    ab_grid = np.zeros([a_length, b_length, 2], dtype=float)
    for i,k in product(np.arange(a_length), np.arange(b_length)):
        j = (i + 1) % a_length
        l = (k + 1) % b_length
        #print(grid[i,k],grid[j,k])
        a = ref_atoms.get_distance(grid[i,k],grid[i,l],mic=True,vector=True)[0]
        b = ref_atoms.get_distance(grid[i,k],grid[j,k],mic=True,vector=True)[1]
        if a <=0:
            a += ref_atoms.cell[0,0]
        if b <=0:
            b = ref_atoms.cell[1,1]
        ab_grid[i,k] = [a,b]
    c = ref_atoms.cell[2,2]
    atoms_grid = make_atoms_grid(ab_grid, c, reflect=False,mode = 'reference')
    return construct_atoms_from_grid(atoms_grid)

# def supercell_from_ref(ref_atoms, relaxed = False, **db_kwargs):
#     #atoms_grid = atoms_grid(a_grid, b_grid, c, reflect, relaxed = relaxed, **db_kwargs)
#     return construct_atoms_from_grid(atoms_grid)