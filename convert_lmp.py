import numpy as np
from copy import deepcopy
from src.data import Data
import sys

data = Data.read_file(sys.argv[1])
data.add_section([[0,0,0,0]], "Bonds")
data.bonds=0
original_data = deepcopy(data)

# O_type_id = np.where(data.atom_type_labels=="O")[0][0]+1
# H_type_id = np.where(data.atom_type_labels=="H")[0][0]+1
Au_type_id = np.where(data.atom_type_labels=="Au")[0][0]+1

gold_atoms = list(filter(lambda x: x.type_id == Au_type_id, data.section_atoms.lines))
# gold_velocities = list(filter(lambda x: x[0] in np.array(gold_atoms).T[0]))

dummy_atom_id = len(data.atom_type_labels) + 1
dummy_bond_id = data.bond_types + 1
data.atom_type_labels = np.append(data.atom_type_labels, "DE")
data.atoms += len(gold_atoms)
data.bonds += len(gold_atoms)
data.bond_types += 1
data.atom_types += 1

data.section_masses.lines = np.append(data.section_masses.lines, [[dummy_atom_id, 1]], axis=0)
data.section_masses.lines[data.section_masses.lines[:,0] == Au_type_id] = [Au_type_id, 195.96657]
# data.section_pair_coeffs.lines = np.append(data.section_pair_coeffs.lines, [[dummy_atom_id, 0.2, 2.62459]], axis=0)
# data.section_bond_coeffs.lines = np.append(data.section_bond_coeffs.lines, [[len(data.section_bond_coeffs.lines), 50]], axis=0)

dummy_atoms = []
bonds = []
for i, atom in enumerate(gold_atoms):
    atom.q = 1.0
    dummy_atom = deepcopy(atom)
    dummy_atom.type_id = dummy_atom_id
    # thwack the dummy atoms on at the end, easier than inserting them next to the gold atoms and shifting other IDs
    # add one to account for i starting at 1
    # IMPORTANT: Assumes that the existing atoms are numbered sequentially with no gaps, achieved via reset_ids in lammps
    dummy_atom.atom_id = int(original_data.atoms) + i + 1
    dummy_atom.q = -1.0
    
    # initialising them with 0 velocity, talk to lore about this
    #  data.section_velocities.lines = np.append(data.section_velocities.lines, [[dummy_atom.atom_id, 0.0, 0.0, 0.0]], axis=0)
    
    # Shift them to be 1 angstrom away from the gold atom centre, stops the simulation from exploding
    # according to eduard
    # print(data.section_velocities.lines.shape)
    dummy_atom.coords = [1/np.sqrt(3) + coord for coord in dummy_atom.coords]
    # Assumes that the dummy bond id is the final one, which should be a safe assumption
    data.section_bonds.lines = np.append(data.section_bonds.lines, [[original_data.bonds+i+1, dummy_bond_id, atom.atom_id, dummy_atom.atom_id]], axis=0)
    print(data.section_bonds.lines.shape)
    
    dummy_atoms.append(dummy_atom)

data.section_atoms.lines = np.append(data.section_atoms.lines, dummy_atoms)

data.write_file(file="/home/esmith/AuNP/resolv_runs/template/Au_111_18_22_6_pol.data")


