from typing import Any
import numpy as np


class Analysis:
    
    unit_z = [0,0,1]
    
    def __init__(bin_height=0.5):
        self.bin_height = bin_height
    
    def __call__(self, old_atoms, new_atoms, *, dump):
        bins = np.linspace(dump.zbounds[0], dump.zbounds[0], self.bin_height)
        
        
    def 