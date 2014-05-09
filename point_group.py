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


    def get_irreducible_representations(self,reducible_characters):

        NG = 1.*len(self.symmetries)

        irreducible_representations = {}

        for rep in self.characters:
            o = N.dot( self.characters[rep] , reducible_characters)/NG

            irreducible_representations[rep] = N.round(o,8)

        return irreducible_representations
