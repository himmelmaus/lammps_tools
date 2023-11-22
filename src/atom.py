import numpy as np
from collections.abc import Iterable

class Atom:
    
    
    def __init__(
        self, atom_id, type_id, x, y, z, *, q=0, molecule_id=None, atom_style=None
    ):
        self.atom_id = atom_id
        self.type_id = type_id
        self.coords = [x, y, z]
        self.q = q
        if molecule_id is not None:
            self.molecule_id = molecule_id
        if atom_style is not None:
            self.__atom_style = atom_style
        else:
            self.__atom_style = "atomic"

    def __str__(self):
        if self.__atom_style == "atomic":
            # atom-ID atom-type x y z
            return f"{self.atom_id} {self.type_id} {' '.join(map(str,self.coords))}\n"

        if self.__atom_style == "full":
            # atom-ID molecule-ID atom-type q x y z
            return f"{self.atom_id} {self.molecule_id} {self.type_id} {self.q} {' '.join(map(str,self.coords))}\n"

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
                tokens = tokens[:-3]
            return cls(
                int(tokens[0]),
                int(tokens[1]),
                *map(float, tokens[2:]),
                atom_style=atom_style,
            )

        if atom_style == "full":
            assert len(tokens) in [7, 10]
            if len(tokens) == 10:
                tokens = tokens[:-3]
            return cls(
                int(tokens[0]),
                int(tokens[2]),
                *map(float, tokens[4:]),
                q=float(tokens[3]),
                molecule_id=int(tokens[1]),
                atom_style=atom_style,
            )
