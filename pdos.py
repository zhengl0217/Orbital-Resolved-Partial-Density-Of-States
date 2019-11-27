"""
This module defines a function to extract orbital resolved partial density of states(pdos) from Materials Project.
"""
import numpy as np
import re
from pymatgen import MPRester

__author__ = "Zheng Li"
__email__ = "zhengl@vt.edu"
__date__ = "Nov. 4, 2019"

def pdos(dos):
    """Extract projected density of state arrays from the doc file and save them into library file. 
    Args:
       dos (doc): doc
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
        # for density of states with spin unpolarized calculations (1 spin chanel), sum all the pdos.
        if 'Spin.up' in pdos[-1]:
            print('spin unpolarized calculation')
            for orbital in orbital_list:
                orbital_index = [inx+1 for inx, item in enumerate(pdos) if orbital in item]
                orbital_str = [re.search('>:(.*)}', pdos[x]).group(1) for x in orbital_index]
                orbital_array = np.array([[float(x) for x in item[item.index('[')+1:item.index(']')].split(',')] for item in orbital_str])
                orbital_sum = orbital_array.sum(axis=0)
                pdos_dict[orbital] = orbital_sum
        # for density of states with spin polarized calculations (2 spin channels), sum all the pdos (up and down channels).
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

if __name__ == "__main__":
    m = MPRester('## API key ###')
    data = m.query('CaTiO3', ["material_id", "spacegroup", "pretty_formula"])
    for entry in data:
        material_id = entry['material_id']
        spacegroup = entry['spacegroup']['crystal_system']
        if spacegroup == 'cubic':
            dos = m.get_dos_by_material_id(material_id)
            partial_dos = pdos(dos)
            print('partial_dos', partial_dos)
