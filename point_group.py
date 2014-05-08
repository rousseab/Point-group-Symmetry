import parse as P
import numpy as N
import os

class PointGroup(object):
    """
    Object which will contain the symmetries and representation 
    of a point group.
    """

    def __init__(self,name, symmetries,representations,characters):

        self.name            = name
        self.symmetries      = symmetries
        self.representations = representations
        self.characters      = characters
