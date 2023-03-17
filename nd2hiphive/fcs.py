"""
Construct a ForceConstantPotential from training data generated previously.

"""

from ase.io import read
from hiphive import ClusterSpace, StructureContainer, ForceConstantPotential
from hiphive.utilities import prepare_structures
from trainstation import Optimizer
from hiphive.force_constants import ForceConstants
from hiphive.cutoffs import estimate_maximum_cutoff


# read structures containing displacements and forces
prim = read('ndhiPhive_prim.xyz')
atoms_ideal = read('ndhiPhive_superc.xyz')

max_cutoff = estimate_maximum_cutoff(atoms_ideal)
print('maximum cutoff is', max_cutoff)

# set up cluster space
cutoffs = [6.35,6.35]
cs = ClusterSpace(prim, cutoffs)
print(cs)
cs.print_orbits()

# ... and structure container
sc = StructureContainer(cs)

with open('ndhiPhive_dispfor.xyz', 'r') as f:
    line = f.readline()
    n_snapshots_to_read = int(line.split()[0])
    n_atoms = int(line.split()[1])
    print("n_atoms %i n_snapshots_to_read %i" % (n_atoms, n_snapshots_to_read))
    for it in range(n_snapshots_to_read):
        print('Reading step', it)
        atoms_tmp = atoms_ideal.copy()
        disps = []
        forces = []
        for it2 in range(n_atoms):
            line = f.readline()
            disp_tmp = [float(fld) for fld in line.split()[0:3]]
            disps.append(disp_tmp)
            force_tmp = [float(fld) for fld in line.split()[3:6]]
            forces.append(force_tmp)
        atoms_tmp.new_array('displacements', disps, dtype=float)            
        atoms_tmp.new_array('forces', forces, dtype=float)
        sc.add_structure(atoms_tmp)

print(sc)

# train model
print('Fitting the parameters...')
opt = Optimizer(sc.get_fit_data())
opt.train()
print(opt)

# construct force constant potential
#fcp = ForceConstantPotential(cs, opt.parameters)
#print(fcp)

# enforce rotational sum rules
from hiphive import enforce_rotational_sum_rules
parameters_rot = enforce_rotational_sum_rules(cs, opt.parameters, ['Huang', 'Born-Huang'])
fcp = ForceConstantPotential(cs, parameters_rot)
print(fcp)

import numpy as np
import matplotlib.pyplot as plt

from ase import Atoms
from ase.build import bulk
from hiphive import ForceConstantPotential

from phonopy import Phonopy
from phonopy.structure.atoms import PhonopyAtoms

# add here the content of ndhiPhive_phonopy.py created by nd2hiphive
atoms_phonopy = PhonopyAtoms(cell=[[  3.2430944463499642,  0.0000000000000000,  0.0000000000000000],
                                   [ -1.6215472231749821,  2.8086021774112981,  0.0000000000000000],
                                   [  0.0000000000000000,  0.0000000000000000, 39.1818803292834161]],
                scaled_positions=[[  0.3333333333333333,  0.6666666666666666,  0.0842169776690000],
                                   [  0.6666666666666669,  0.3333333333333333,  0.2422572444300000],
                                   [  0.3333333333333333,  0.6666666666666666,  0.1995238704700000],
                                   [  0.3333333333333333,  0.6666666666666666,  0.2849805655400000],
                                   [  0.6666666666666669,  0.3333333333333333,  0.1234833703200000],
                                   [  0.6666666666666669,  0.3333333333333333,  0.0448833825980000]],
                symbols=['Mo']*2+['Se']*2+['S']*2)

phonopy = Phonopy(atoms_phonopy, supercell_matrix=[[4,0,0],
                                                   [0,4,0],
                                                   [0,0,1]],
                             primitive_matrix=None)

##

supercell = phonopy.get_supercell()
supercell = Atoms(cell=supercell.cell, numbers=supercell.numbers, pbc=True,
                  scaled_positions=supercell.get_scaled_positions())
fcs = fcp.get_force_constants(supercell)

print(fcs)
fcs.write_to_phonopy('FORCE_CONSTANTS','text')
fcs.write_to_phono3py('fc3.hdf5')

