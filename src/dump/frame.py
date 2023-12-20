import numpy as np
from copy import copy


class Frame:
    
    def __init__(self, frame_text):
        self.text = copy(frame_text)