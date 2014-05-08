import unittest
import numpy as N

import sys
sys.path.append('../')

import os
from fortran_parser import PointGroupFortranParser

class TestFortranParser(unittest.TestCase):

    def setUp(self):
        # This may not be the best way to test for Fortran 
        self.fortran_path     = '/Users/Bruno/work/depository/abinit-7.6.3/src/43_ptgroups/'
        self.fortran_filename = self.fortran_path+'ptg_Oh.F90'

        self.G     = PointGroupFortranParser(self.fortran_filename)

    def test_directory_exists(self):
        # make sure the fortran files directory exists
        self.assertTrue( os.path.isdir(self.fortran_path) )

    def test_file_exists(self):
        # make sure the fortran file exists
        self.assertTrue( os.path.isfile(self.fortran_filename) )

    def test_check_symmetry_matrix(self):
        self.sym37 = N.array([[0, -1, 0],[-1, 0, 0], [0, 0, 1]])

        isym = 36 # Careful! python indices start at zero

        error = N.linalg.norm(self.G.symmetries[isym,:,:]-self.sym37)

        self.assertTrue( error < 1e-8)

    def test_number_of_irr_reps(self):
        self.assertEqual( len(self.G.representations), 10)

    def test_irr_reps(self):

        irr_T2g_9 = N.array([ [0.0, 0.0, 1.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]])
        isym = 8 # Careful! python indices start at zero

        error = N.linalg.norm( self.G.representations["T2g"][isym,:,:]-irr_T2g_9)
        self.assertTrue( error < 1e-8)


if __name__ == '__main__':
    unittest.main()
