#!/usr/bin/env python
import numpy as np
import re
from pymatgen import MPRester, Composition, Structure
import matplotlib.pyplot as plt

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

properties = ['task_id', "energy", "energy_per_atom", "volume",
                            "formation_energy_per_atom", "nsites",
                            "unit_cell_formula", "pretty_formula",
                            "is_hubbard", "elements", "nelements",
                            "e_above_hull", "hubbards", "is_compatible",
                            "spacegroup", "task_ids", "band_gap", "density",
                            "icsd_id", "icsd_ids", "cif", "total_magnetization",
                            "material_id", "oxide_type", "tags", "elasticity"]

#oxide = 'CaTiO3'
#oxide = 'SrNiO3'
#oxide = 'SrAgO3'
#oxide = 'LaCoO3'
#oxide = 'SrMnO3'
oxide = 'BaCoO3'
#oxide = 'LaFeO3'
#oxide = 'LaCuO3'
m = MPRester('FfPyT6VLIvhlRE9c')
data = m.query('BaRhO3', properties)
for entry in data:
    spacegroup = entry['spacegroup']['crystal_system']
    print('spacegroup', spacegroup)
    material_id = entry['material_id']
    #if material_id == mp-2723:
    if spacegroup == 'cubic':#'hexagonal':#'cubic':
        try:
            formula = entry['pretty_formula']
            dos = m.get_dos_by_material_id(material_id)
            lib = pdos(dos)
            lib_dos = sum(lib.values())
            efermi = dos.efermi                                                                                                               
            print('efermi', efermi)                                                                                                        
            energies = dos.energies                                                                                                        
        except:
            continue

densities = sum(dos.densities.values())                
densities_d = (lib['Orbital.dz2']+ lib['Orbital.dx2']+ lib['Orbital.dxy']+ lib['Orbital.dxz']+ lib['Orbital.dyz'])
densities_p = (lib['Orbital.px'] + lib['Orbital.py'] + lib['Orbital.pz'])
densities_s = lib['Orbital.s']
densities_tot = densities_d + densities_p + densities_s
plt.plot(densities, energies-efermi, c='grey', label = 'total')
plt.plot(densities_d, energies-efermi, c='b', label= 'd')
plt.plot(densities_s, energies-efermi, c='g', label= 's')
plt.plot(densities_p, energies-efermi, c='r', label= 'p')
#plt.plot(densities, energies-efermi, c='k', label = 'total_model')                                                                
plt.xlim(0, 0.5)
plt.ylim(-10, 20)
plt.legend()
plt.title(formula + ',' +spacegroup)
plt.show()
