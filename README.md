Point-group-Symmetry
====================

Python code which will determine the symmetry representations of combinations of atomic orbitals on different sites arranged according to a given point group.

   The goal of this project is to have a simple python code
which will determine the symmetry representations of
combinations of orbitals on different sites arranged according
to a given point group.

For example, a question this code seeks to address is
"what are the irreducible representations describing the
symmetry of p-orbitals on the vertices of an octahedron (Oh symmetry)?"

This will be useful, for example, in understanding density-of-states
peaks of abinitio computations pertaining to transition metal oxides,
which often have a transition metal surrounded by oxygens on the
vertices of a polyhedron. Simple crystal field theory describes how the
atomic state should hybridize.

This is neither a hard or a new problem. However, I 
cannot find a simple key-in-hand solution. I thus take this opportunity
to start this little project to gain experience with Python, unittest, git, etc...

This first commit pertains to extracting point group information
from fortran code, freely available from the ABINIT distribution,
and making it available to my python module.
