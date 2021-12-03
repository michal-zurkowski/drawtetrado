# distutils: language = c++

from optimizer cimport SolveFailsafe

def solve(edges, rotations, alignments):
    return SolveFailsafe(edges, rotations, alignments).positions
