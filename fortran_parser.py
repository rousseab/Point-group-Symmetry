import parse as P
import numpy as N
import os

class PointGroupFortranParser(object):
    """
    Object which will contain the symmetries and representation 
    of a point group.
    """

    def __init__(self,fortran_filename):

        self.filename   = fortran_filename

        self.group_name = os.path.basename(fortran_filename).replace('ptg_','').replace('.F90','')

        self.extract_group_information_from_Fortran()

    def parsing_string_to_matrix(self, string):

        list_str = string.strip().split(',')
        list_complex = []
        for str in list_str:
            list_complex.append(complex(str.replace('*','')))

        dim = int(N.sqrt(len(list_complex)))

        matrix = N.array(list_complex).reshape(dim,dim)

        return matrix

    def extract_group_information_from_Fortran(self):
        """
        Extract useful parameters. Only parse through the file once!
        """

        file  = open(self.filename,'r')
        lines = file.readlines()
        file.close()

        for line in lines:
            # extract number of symmetries
            result = P.search(' nsym = {:d}',line)
            if result != None:
                nsym = result[0]
                sym  = N.zeros([nsym,3,3])

            # extract number of classes
            result = P.search(' nclass = {:d}',line)
            if result != None:
                nclass = result[0]
                dic_class = {}

            # extract symmetry matrices. By the time we get a match, the 
            # the sym array will have been generated!
            result = P.search('sym(:,:,{:d}) = RESHAPE( (/{}/) ,(/3,3/) )',line)
            if result != None:
                isym   = result[0]-1
                matrix = self.parsing_string_to_matrix(result[1])
                sym[isym,:,:] = N.real(matrix)


            # read in the representation names
            result = P.search('Irr({:d})%name = "{}"',line)
            if result != None:
                key = result[0] # number of irr. rep.
                dic_class[key] = [result[1]] # name of irr. rep.

            # read in the representation dimension
            result = P.search('Irr({:d})%dim = {:d}',line)
            if result != None:
                irep= result[0] # number of irr. rep.
                dim = result[1]
                dic_class[irep].append(complex(0.,0.)*N.zeros([nsym,dim,dim])) # add zero-array for representation

            # read in the representation matrix
            result = P.search('Irr({:d})%mat(:,:,{:d}) =  RESHAPE( (/{}/),{}',line)
            if result != None:
                irep = result[0]
                isym = result[1]-1
                str  = result[2]
                matrix = self.parsing_string_to_matrix(str)

                dic_class[irep][1][isym,:,:] = matrix

        # feed it to the object

        self.symmetries = sym

        self.representations = {}
        self.characters      = {}

        for key in dic_class.keys():
            name, rep = dic_class[key]
            self.representations[name] = rep 
            self.characters[name] = rep.trace(axis1=1,axis2=2)

