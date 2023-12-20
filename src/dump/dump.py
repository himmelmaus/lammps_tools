import numpy as np
from copy import copy, deepcopy

from ..atom import Atom
from ..data import Data
from ..helpers import readlines
# from .frame import Frame


class Dump:
    
    def __init__(self, trajectory_file, timestep=250, data_file = None):
        
        """
        Don't need to do that much here given can't actually eager process the trajectory file
        """
        
        self.filename = trajectory_file
        self.timestep = timestep
        
        file = open(self.filename, "r")
        header_lines = readlines(file, 9) # file.readlines(hint=9)
        
        self.init_timestep = int(header_lines[1])
        self.current_step = self.init_timestep
        self.n_atoms = int(header_lines[3])
        self.xbounds = list(map(float, header_lines[5].split()))
        self.ybounds = list(map(float, header_lines[6].split()))
        self.zbounds = list(map(float, header_lines[7].split()))
            
        if data_file:
            self.data_file = data_file
            self.init_data = Data.read_file(data_file)
            self.atoms = np.array(self.init_data.section_atoms.lines)
            self.atom_types = self.init_data.atom_type_labels[:,0]
            # Convert to textual atom type labels
            self.atoms[:, 2] = [self.init_data.atom_type_labels[self.init_data.atom_type_labels[:,1] == atom[2]][0] for atom in self.atoms]
        else:
            body_lines = readlines(file, self.n_atoms) # file.readlines(hint=self.n_atoms)
            # Defaults to atomic which matches the dump file encoding
            self.atoms = np.array(list(map(Atom.parse_line, body_lines)))
            self.atom_types = list(set(map(lambda x: x.type_id, self.atoms)))
        atom_groups = {}
        for atom_type in self.atom_types:
            atom_groups[atom_type] = [atom.atom_id for atom in self.atoms if atom.type_id == atom_type]
        file.close()
        
    def ingest(self, start=0, frame_funcs = []):
        """Main function for iterating through the frames of a file and processing based off of that"""
        file = open(self.filename, "r")
        func_outputs = [[] for _ in range(len(frame_funcs))]
        
        if start > 0:
            assert start/self.timestep == int(start/self.timestep)
            skip = self.lines_per_frame*(start/self.timestep)
            # Just skips us forward this many lines
            readlines(file, skip)
          
        while True:  
            old_atoms = deepcopy(self.atoms)
            atom_lines = readlines(file, self.lines_per_frame)
            
            assert int(atom_lines[3]) == self.n_atoms
            self.current_step = int(atom_lines[1])
            self.xbounds = list(map(float, atom_lines[5]))
            self.ybounds = list(map(float, atom_lines[6]))
            self.zbounds = list(map(float, atom_lines[7]))
            
            if len(atom_lines[0]) == 0:
                break
            
            self.atoms = np.array(list(map(lambda a, l: a.update_position(l), zip(self.atoms, atom_lines[9:]))))
            new_atoms = deepcopy(self.atoms)
            
            for i, func in enumerate(frame_funcs):
                func_outputs[i].append(func(old_atoms, new_atoms, dump=self))
            
        return 
            
        
        
        
    
    
    @property
    def lines_per_frame(self):
        return 9 + self.n_atoms
    
    @property
    def molecules(self):
        """
        Returns a list of the atoms which share molecule IDs, the first element (e.g. molecule ID 0) is always the atoms
        not in molecules
        """
        
        mol_ids = set(map(lambda x: x.molecule_id, self.init_data.section_atoms.lines))
        
        molecules = {}
        for id in mol_ids:
            molecules[id] = [atom.atom_id for atom in self.init_data.section_atoms.lines if atom.molecule_id == id]
            
        return molecules