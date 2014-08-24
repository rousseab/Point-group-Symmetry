import numpy as N

from point_group import PointGroup
from fortran_parser import load_pickle_file

from orbital_representation import get_representation_matrix


class SymmetricSites(object):
    """
    This object will represent the orbit of some generating positions, along with the orbitals decorating
    these sites.

    The genetors should be a list of (position, angular momentum) couples.
    """
    def __init__(self,Group, generators):

        self.Group      = Group
        self.generators = generators

    def extract_irreducible_representations(self):

        self.setup_representation()
        self.build_representation()

        self.characters = N.trace(self.representation,axis1=1,axis2=2)

        list_irr = self.Group.get_irreducible_representations(self.characters)

        return list_irr 


    def setup_representation(self):

        self.representation_dimension = 0

        self.list_positions = []
        self.list_l         = []

        for generator in self.generators: 

            position = generator[0]
            l        = generator[1]

            for sym in self.Group.symmetries:

                new_p = N.dot(sym,position)

                if self.new_vector(new_p):
                    self.list_positions.append(new_p)
                    self.list_l.append(l)
                    self.representation_dimension += 2*l+1

        self.list_positions = N.array(self.list_positions)


    def build_representation(self):

        self.representation = complex(0.,0.)*N.zeros([len(self.Group.symmetries), self.representation_dimension, self.representation_dimension])


        for isym, sym in enumerate(self.Group.symmetries):

            for I, (p, l) in enumerate(zip(self.list_positions,self.list_l)):

                new_p = N.dot(sym,p)

                J = self.find_position_index(new_p)

                D = get_representation_matrix(l, sym)

                
                i1, i2 = self.find_dimension_range(I)
                j1, j2 = self.find_dimension_range(J)

                self.representation[isym,i1:i2,j1:j2] = D


    def new_vector(self,vector):
        tol = 1e-8

        for old_vector in self.list_positions:
            if N.linalg.norm(old_vector-vector) < tol:
                return False

        return True

    def find_position_index(self,vector):
        tol = 1e-8

        for J, v in enumerate(self.list_positions):
            if N.linalg.norm(vector-v) < tol:
                return J

    def find_dimension_range(self,I):

        i1 = 0
        i2 = 0

        for J, l in enumerate(self.list_l):
            i2 += 2*l+1
            if J == I:
                break
            i1 += 2*l+1


        return i1, i2



