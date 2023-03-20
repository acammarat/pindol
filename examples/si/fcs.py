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
cutoffs = [11.0,9.5]
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
fcp = ForceConstantPotential(cs, opt.parameters)
print(fcp)

import numpy as np
import matplotlib.pyplot as plt

from ase import Atoms
from ase.build import bulk
from hiphive import ForceConstantPotential

from phonopy import Phonopy
from phonopy.structure.atoms import PhonopyAtoms

# add here the content of ndhiPhive_phonopy.py created by nd2hiphive
atoms_phonopy = PhonopyAtoms(cell=[[  0.0000000000000000,  2.7241199674880803,  2.7241199674880803],
                                   [  2.7241199674880803,  0.0000000000000000,  2.7241199674880803],
                                   [  2.7241199674880803,  2.7241199674880803,  0.0000000000000000]],
                scaled_positions=[[  0.2500000000000000,  0.2500000000000000,  0.2500000000000000],
                                   [  0.0000000000000000,  0.0000000000000000,  0.0000000000000000]],
                symbols=['Si']*2)

phonopy = Phonopy(atoms_phonopy, supercell_matrix=[[6,0,0],
                                                   [0,6,0],
                                                   [0,0,6]],
                             primitive_matrix=None)


##

supercell = phonopy.get_supercell()
supercell = Atoms(cell=supercell.cell, numbers=supercell.numbers, pbc=True,
                  scaled_positions=supercell.get_scaled_positions())
fcs = fcp.get_force_constants(supercell)

print(fcs)
fcs.write_to_phonopy('FORCE_CONSTANTS','text')
fcs.write_to_phono3py('fc3.hdf5')

