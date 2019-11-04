import numpy as np
import re
from pymatgen import MPRester, Composition, Structure

def pdos(dos):
    """Extract projected density of state arrays from the doc file and save them into library file. 
    Args:
       dos (doc): doc
       size (int): number of energy levels for the specific material electronic structure 
    Returns:
       dict: library of pdos arrays obtained by summing all the atomic pdos 
    """
    if dos is None:
        raise TypeError("The dos object for this particular structure is not provided by Materials Project") 
    else:
        # split the dos text file into list of pdos elements
        pdos = str(dos.__dict__['pdos'].values()).split('<')
        # detect the key index of specific orbital and extract the corresponding numeric arrays   
        orbital_list = ['Orbital.s', 'Orbital.py', 'Orbital.pz', 'Orbital.px', 'Orbital.dz2', 'Orbital.dx2',\
                            'Orbital.dxy', 'Orbital.dxz', 'Orbital.dyz']
        # initialize the pdos_dict with 0 for each orbital
        pdos_dict = {}
        for item in orbital_list:
            pdos_dict[item] = np.zeros(len(dos.energies))
        # search all the orbital-resolved partial density of states and calculate the sum    
        # for density of states with spin unpolarized calculations (1 spin chanel)
        if 'Spin.up' in pdos[-1]:
            print('spin unpolarized calculation')
            for orbital in orbital_list:
                orbital_index = [inx+1 for inx, item in enumerate(pdos) if orbital in item]
                orbital_str = [re.search('>:(.*)}', pdos[x]).group(1) for x in orbital_index]
                orbital_array = np.array([[float(x) for x in item[item.index('[')+1:item.index(']')].split(',')] for item in orbital_str])
                orbital_sum = orbital_array.sum(axis=0)
                pdos_dict[orbital] = orbital_sum
        # for density of states with spin polarized calculations (2 spin channels)
        if 'Spin.down' in pdos[-1]:
            print('spin polarized calculation')
            for orbital in orbital_list:
                orbital_index_up = [inx+1 for inx, item in enumerate(pdos) if orbital in item]
                orbital_index_down = [inx+2 for inx, item in enumerate(pdos) if orbital in item]
                orbital_str_up = [re.search('>:(.*),', pdos[x]).group(1) for x in orbital_index_up]
                orbital_str_down = [re.search('>:(.*)]}', pdos[x]).group(1) for x in orbital_index_down]
                orbital_array_up = np.array([[float(x) for x in item[item.index('[')+1:item.index(']')].split(',')] for item in orbital_str_up])
                orbital_array_down = np.array([[float(x) for x in item[item.index('[')+1:].split(',')] for item in orbital_str_down])
                orbital_array = np.concatenate((orbital_array_up, orbital_array_down), axis=0)
                orbital_sum = orbital_array.sum(axis=0)
                pdos_dict[orbital] = orbital_sum
    
        return pdos_dict


m = MPRester('## API key ###')
oxide = 'CaTiO3'
properties = ["material_id", "spacegroup", "pretty_formula"]
data = m.query(oxide, properties)
for entry in data:
    spacegroup = entry['spacegroup']['crystal_system']
    material_id = entry['material_id']
    if spacegroup == 'cubic':
        formula = entry['pretty_formula']
        dos = m.get_dos_by_material_id(material_id)
        partial_dos = pdos(dos)
        print('partial_dos', partial_dos)