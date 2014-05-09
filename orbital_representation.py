
from math import factorial as fac
import numpy as N

from rotation_matrix import R_to_euler
from point_group import PointGroup
from fortran_parser import load_pickle_file


class Row(object):
    """
    This simple object will represent a row in a representation matrix, where the
    index is given by Row[m], m = -l,...,l.
    """
    def __init__(self,l):
        self.l    = l
        self._rep = complex(0.,0.)*N.zeros([2*l+1])

    def __getitem__(self,index):
        if index >= -self.l and index <= self.l:
            return self._rep[index+self.l]
        else:
            print 'ERROR! Trying to access element which does not exits!'

    def __setitem__(self,index,value):
        if index >= -self.l and index <= self.l:
            self._rep[index+self.l] = value
        else:
            print 'ERROR! Trying to access element which does not exits!'

class D_matrix(object):
    """
    This object will represent the d matrices as given in the paper

    "A fast and stable method for rotating spherical harmonic expansions",
    Journal of Computational Physics 228 (2009) 5621-5627

    This is just a convenience object to make the recurrence easier to write down and understand.
    """
    def __init__(self,l):
        self.l    = l
        self._rep = []
        for n in N.arange(2*l+1):
            self._rep.append(Row(l))


    def __getitem__(self,index):
        if index >= -self.l and index <= self.l:
            return self._rep[index+self.l]
        else:
            print 'ERROR! Trying to access element which does not exits!'


    def get_matrix(self):

        Matrix = complex(0.,0.)*N.zeros([2.*self.l+1,2*self.l+1])

        m_range = N.arange(-self.l,self.l+1)

        for im1, m1 in enumerate(m_range):
            for im2, m2 in enumerate(m_range):

                Matrix[im1,im2] = self._rep[im1][m2]
        return Matrix

def get_representation_matrix(l, Group_Operation):
    """
    this function computes the d matrix, as presented in 

    "A fast and stable method for rotating spherical harmonic expansions", 
    Journal of Computational Physics 228 (2009) 5621-5627

    Note that the algorithm uses the explicit sum quoted as being from Wigner. This is known to be unstable, but we will
    hardly need to consider cases beyond l=2, and this algorithm seems the easiest to code.
    """

    # The sign indicates if the operation is a pure rotation, or a rotation combined with inversion
    sign = N.linalg.det(Group_Operation)

    R    = sign*Group_Operation

    # Only a pure rotation has Euler angles!
    alpha, beta, gamma = R_to_euler(R, 'zyz')

    co = N.cos(beta/2.)
    si = N.sin(beta/2.)

    D = D_matrix(l)

    m_range = N.arange(-l,l+1)

    for mp in m_range:
        for m in m_range:

            exp = N.exp(-1j*m*(alpha+gamma))

            prefactor = (-1.)**(mp-m)*N.sqrt( fac(l+mp)*fac(l-mp)*fac(l+m)*fac(l-m) )

            s_min = N.max( [0,m-mp])
            s_max = N.min( [l+m,l-mp])

            #print 'mp = %i,  m = %i'%(mp,m)
            #print '     s_min = %i,  s_max = %i'%(s_min,s_max)
                
            for s in N.arange(s_min,s_max+1):

                den = fac(l+m-s)*fac(s)*fac(mp-m+s)*fac(l-mp-s)

                e1  = 2*(l-s)+m-mp
                e2  = 2*s-m+mp

                D[mp][m] += prefactor*(-1.)**s*co**e1*si**e2/den*exp


    Representation_Matrix = sign**l * D.get_matrix()

    return Representation_Matrix 




