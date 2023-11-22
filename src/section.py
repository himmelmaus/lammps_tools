import numpy as np
from collections.abc import Iterable

from .atom import Atom
from .helpers import num_str, snake_case

class Section:
    
    def __init__(self, lines, title, *, style):
        self.title = title.strip()
        if style is not None:
            self.style = style.strip()

        if title != "Atoms":
            # Might come back to this and make each entry its own class instance, and make Atom inherit from this
            # but I'm not sure how useful that'd be. May also just scrap Atom, but maintaining an inbuilt functionality
            # for changing atom style and such would be handy.
            if isinstance(lines, str):
                self.lines = np.array(
                    list(map(lambda x: list(map(float, x.strip().split(" "))), lines))
                )
            elif isinstance(lines, Iterable):
                self.lines = np.array(
                    lines
                )
            else:
                raise ValueError("Section.lines must be either a string or an iterable.")

        else:
            self.lines = np.array(
                list(map(lambda x: Atom.parse_line(x, atom_style=style), lines))
            )
            
    def print_lines(self):
        
        if getattr(self, "style", None) is not None:
            outlines = [f"{self.title} # {self.style}\n\n"]
        else:
            # I think it's fine but this is giving me passing by reference flashbacks
            outlines = [self.title, "\n\n"]
            
        if self.title != "Atoms":
            outlines.extend(list(map(lambda x: " ".join(map(num_str, x)) + "\n", self.lines)))
            
        else:
            outlines.extend(list(map(str, self.lines)))

        return outlines
    
    @property
    def section_id(self):
        # Converts to snake case
        return snake_case(self.title)