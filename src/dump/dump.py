from copy import copy, deepcopy
import os


import numpy as np

from ..atom import Atom
from ..data import Data
from ..helpers import readlines
# from .frame import Frame


class Dump:

    init_data = None
    skip = None
    header_lines = 9
    
    has_timestep = True
    current_step = None
    
    data_class = Data
    dump_style = "lammps"

    def __init__(self, trajectory_file, timestep=250, data_file = None, atom_type_labels=True):
    
        """
        Don't need to do that much here given can't actually eager process the trajectory file
        """
        
        self.filename = trajectory_file
        self.timestep = timestep
        self.atom_type_labels = atom_type_labels
        
        file = open(self.filename, "r")
        header_lines, _ = readlines(file, self.header_lines) # file.readlines(hint=9)
        
        self.process_header(header_lines)
        
        self.init_timestep = self.current_step
            
        if data_file:
            self.data_file = data_file
            self.init_data = self.data_class.read_file(data_file)
            self.atoms = np.array(self.init_data.section_atoms.lines)
            self.atoms = sorted(self.atoms, key=lambda a: a.atom_id)
            self.atom_types = self.init_data.atom_type_labels
            if atom_type_labels:
                try:
                    # Convert to textual atom type labels
                    for atom in self.atoms:
                        atom.type_id = self.init_data.atom_type_labels[atom.type_id - 1]
                except Exception as exc:
                    raise AttributeError("Atom type labels not found in input data file") from exc
        else:
            body_lines, _ = readlines(file, self.n_atoms) # file.readlines(hint=self.n_atoms)
            self.process_body(body_lines)
            
        atom_groups = {}
        for atom_type in self.atom_types:
            atom_groups[atom_type] = [atom.atom_id for atom in self.atoms if atom.type_id == atom_type]
        file.close()
        
    def ingest(self, start=0, frame_funcs = [], save_output = True, max_steps = None, skip = None, show_step = False):
        """Main function for iterating through the frames of a file and processing based off of that"""
        file = open(self.filename, "r")
        func_outputs = [[] for _ in range(len(frame_funcs))]
        
        
        
        if start > 0:
            assert start/self.timestep == int(start/self.timestep)
            jump = self.lines_per_frame*(start/self.timestep)
            # Just skips us forward this many lines
            readlines(file, jump)
            # print("hello")

        if skip is not None:
            self.skip = skip
        
          
        while True:  

            old_atoms = deepcopy(self.atoms)

            if self.skip:
                atom_lines, _ = readlines(file, self.lines_per_frame*int(self.skip/self.timestep))
                atom_lines = atom_lines[-self.lines_per_frame:]
            else:
                atom_lines, _ = readlines(file, self.lines_per_frame)

            if file.tell() == os.fstat(file.fileno()).st_size:
                print("Reached EOF")
                break

            if max_steps and self.current_step >= max_steps:
                print("Reached step limit")
                break
            
            try:
                if self.dump_style == "lammps":
                    assert int(atom_lines[3]) == self.n_atoms
                elif self.dump_style == "cp2k":
                    assert int(atom_lines[0]) == self.n_atoms
                assert self.n_atoms == len(atom_lines[self.header_lines:])
            except AssertionError as e:
                raise ValueError("Different number of atoms from beginning") from e
            
            self.process_header(atom_lines[:self.header_lines])

            if show_step:
                print("Timestep: ", self.current_step)
            
            if len(atom_lines[0]) == 0:
                break
            
            self.atoms = np.array(list(map(lambda a: a[0].update_position(a[1]), zip(self.atoms, atom_lines[self.header_lines:]))))
            new_atoms = deepcopy(self.atoms)
            
            if save_output:
                for i, func in enumerate(frame_funcs):
                    func_outputs[i].append(func(old_atoms, new_atoms, dump=self))
                    continue
            for func in frame_funcs:
                func(old_atoms, new_atoms, dump=self)
        
        # print(func_outputs)
        return 
            
        
    def process_header(self, header_lines):
        
        self.current_step = int(header_lines[1].strip("\n"))
        self.n_atoms = int(header_lines[3])
        self.xbounds = list(map(float, header_lines[5].split()))
        self.ybounds = list(map(float, header_lines[6].split()))
        self.zbounds = list(map(float, header_lines[7].split()))
        
        return
    
    
    def process_body(self, body_lines):
        
        # Defaults to atomic which matches the dump file encoding
        self.atoms = np.array(list(map(self.atom_parser, body_lines)))
        self.atom_types = list(set(map(lambda x: x.type_id, self.atoms)))
        
        return
    
    def atom_parser(self, body_line):
        # to allow for DRY code in the CP2K case
        return Atom.parse_line(body_line)
    
    @property
    def lines_per_frame(self):
        return self.header_lines + self.n_atoms
    
    @property
    def mol_ids(self):
        return list(set(map(lambda x: x.molecule_id, self.init_data.section_atoms.lines)))
    
    @property
    def molecules(self):
        """
        Returns a list of the atoms which share molecule IDs, the first element (e.g. molecule ID 0) is always the atoms
        not in molecules
        """
        
        assert self.init_data is not None, "A lammps data file must be provided to self.init_data in order to access molecule information"
        
        molecules = {}
        for id in self.mol_ids:
            molecules[id] = [atom.atom_id for atom in self.init_data.section_atoms.lines if atom.molecule_id == id]
            
        return molecules
    

class CP2KDump(Dump):
    
    header_lines = 2
    has_timestep = False 
    dump_style = "cp2k"
    
    def process_header(self, header_lines):
        
        if not self.has_timestep and self.current_step is None:
            self.current_step = 0
        else:
            self.current_step += 1
            
        self.n_atoms = int(header_lines[0])
        
        # temporary hard coding while I fix it up
        self.xbounds = np.array([0, 15.0])
        self.ybounds = np.array([0, 15.0])
        self.zbounds = np.array([0, 15.0])
        
    def atom_parser(self, body_line):
        # to allow for DRY code in the CP2K case
        return Atom.parse_line(body_line, atom_style="CP2K")