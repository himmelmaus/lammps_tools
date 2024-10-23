import numpy as np
from collections.abc import Iterable

class Atom:
    
    
    def __init__(
        self, atom_id, type_id, x, y, z, *, q=0, molecule_id=None, atom_style=None, image_flags=[0, 0,0 ]
    ):
        self.atom_id = atom_id
        self.type_id = type_id
        self.coords = np.array([x, y, z])
        self.q = q
        self.image_flags = image_flags
        if molecule_id is not None:
            self.molecule_id = molecule_id
        if atom_style is not None:
            self.__atom_style = atom_style
        else:
            self.__atom_style = "atomic"

    def __str__(self):
        if self.__atom_style == "atomic":
            # atom-ID atom-type x y z
            return f"{self.atom_id} {self.type_id} {' '.join(map(str,self.coords))} {' '.join(self.image_flags)}\n"

        if self.__atom_style == "full":
            # atom-ID molecule-ID atom-type q x y z
            return f"{self.atom_id} {self.molecule_id} {self.type_id} {self.q} {' '.join(map(str,self.coords))} {' '.join(map(str,self.image_flags))}\n"
        
    def update_position(self, line):
        
        if isinstance(line, str):
            tokens = line.strip().split()
        elif isinstance(line, Iterable):
            tokens = np.copy(line)
        else:
            raise ValueError("Atom.line must be either a string or an iterable.")
        
        if self.__atom_style.lower() == "cp2k":
            assert len(tokens) == 4
            
            assert tokens[0] == self.type_id
            self.coords = np.array(list(map(float, tokens[1:])))
            
        elif True or self.__atom_style == "atomic":
            assert len(tokens) in [5, 8]
            if len(tokens) == 8:
                images = tokens[-3:]
                tokens = tokens[:-3]
            elif len(tokens) == 5:
                images = [0,0,0]
            
            assert int(tokens[0]) == self.atom_id
            self.coords = np.array(list(map(float, tokens[2:])))
            self.image_flags = images
            
        return self
            

    @classmethod
    def parse_line(cls, line, atom_style="atomic"):
        
        # Assumes we aren't caring about image flags
        
        if isinstance(line, str):
            tokens = line.strip().split()
        elif isinstance(line, Iterable):
            tokens = np.copy(line)
        else:
            raise ValueError("Atom.line must be either a string or an iterable.")

        if atom_style == "atomic":
            
            assert len(tokens) in [5, 8]
            if len(tokens) == 8:
                images = tokens[-3:]
                tokens = tokens[:-3]
            elif len(tokens) == 5:
                images = [0,0,0]
                
            return cls(
                int(tokens[0]),
                int(tokens[1]),
                *map(float, tokens[2:]),
                atom_style=atom_style,
                image_flags = images
            )
        
        elif atom_style.lower() == "cp2k":
            
            assert len(tokens) == 4
            
            return cls(
                # TODO fixmeeeeeeeeeeee, don't know what to do when we don't get given atom IDs
                0,
                tokens[0],
                *map(float, tokens[1:]),
                atom_style=atom_style
            )

        elif atom_style == "full":
            assert len(tokens) in [7, 10]
            if len(tokens) == 10:
                images = tokens[-3:]
                tokens = tokens[:-3]
                return cls(
                int(tokens[0]),
                int(tokens[2]),
                *map(float, tokens[4:]),
                q=float(tokens[3]),
                #TODO FIXME FIX THIS!!!!
                molecule_id=0,
                atom_style=atom_style,
                image_flags = images
            )
            # else:
            #     images = [0,0,0]
            return cls(
                int(tokens[0]),
                int(tokens[2]),
                *map(float, tokens[4:]),
                q=float(tokens[3]),
                # TODO FIXME FIX THIS!!!!
                molecule_id=0,
                atom_style=atom_style,
                # image_flags = images
            )