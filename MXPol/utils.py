import numpy as np
from itertools import product

def atoms_mapping(atoms_key, atoms_val):
    '''
    assuming two atoms are the same. 
    Find mapping from one to another
    '''
    _, ind_key, ind_val = intersect2D(
        atoms_key.get_positions(),atoms_val.get_positions()
        )
    return {k:v for k,v in zip(ind_key, ind_val)}

# def sorted_to_unsorted_mapping(atoms_sorted, atoms_unsorted):
#     pos_sorted =atoms_sorted.get_positions()
#     pos_unsorted =atoms_unsorted.get_positions()
    
#     _, sort_ind, unsort_ind = intersect2D(pos_sorted,pos_unsorted)
#     return {k:v for k,v in zip(unsort_ind,sort_ind)}

# def read_mapping(arr_element):
#     return mapping[arr_element]
    
def intersect2D(Array_A, Array_B):
  """
  Find row intersection between 2D numpy arrays, a and b.
  Returns another numpy array with shared rows and index of items in A & B arrays
  obtained from https://gist.github.com/RashidLadj/bac71f3d3380064de2f9abe0ae43c19e
  """
  # [[IDX], [IDY], [value]] where Equal
  # ''' Using Tuple ''' #
  #IndexEqual = np.asarray([(i, j, x) for i,x in enumerate(Array_A) for j, y in enumerate (Array_B)  if(tuple(x) == tuple(y))],dtype='object').T
  
  # ''' Using Numpy array_equal ''' #
  IndexEqual = np.asarray([(i, j, x) for i,x in enumerate(Array_A) for j, y in enumerate (Array_B)  if(np.array_equal(x, y))],dtype='object').T
  
  idx, idy, intersectionList = (IndexEqual[0], IndexEqual[1], IndexEqual[2]) if len(IndexEqual) != 0 else ([], [], [])

  return intersectionList, idx, idy

def element_index_grid(a_length, b_length, offset = 2, mapping = {}):
    natoms_unit_cell = 4
    grid = natoms_unit_cell* np.arange(a_length*b_length).reshape(a_length, b_length)+offset
        
    if len(mapping) != 0:
        for i,j in product(range(a_length), range(b_length)):
            grid[i,j] = mapping[grid[i,j]]
    return grid