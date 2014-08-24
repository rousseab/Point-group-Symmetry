import numpy as N
from orbital_representation import get_representation_matrix

from rotation_matrix import R_to_euler
from point_group import PointGroup
from fortran_parser import load_pickle_file
from sites import SymmetricSites


groups = load_pickle_file()

Group = groups['Oh']

#l = 1
#generators = [ ([0.0,0.0,1.0], l) ]

#l = 2
#generators = [ ([0.0,0.0,0.0], l) ]

l = 2
generators = [([0.0,0.0,1.0], 1), ([0.0,0.0,0.0], 2) ]




sites      = SymmetricSites( Group, generators )


list_irr = sites.extract_irreducible_representations()

print '####################'
print '# l = %i'%l
print '####################'
for irr in list_irr:
    weight = list_irr[irr]

    x = N.real(weight)
    y = N.imag(weight)

    str = '          '
    if N.abs(weight) > 1e-8:
        if N.abs(x) > 1e-8:
            str += '%4.3f '%x
        if N.abs(y) > 1e-8:
            str += '%+4.3f j'%y
    
        str += '  %s'%irr
        print str

