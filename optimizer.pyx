# distutils: language = c++

from optimizer cimport Solution, SolveFailsafe

def solve(edges, rotations, alignments):
    return SolveFailsafe(edges, rotations, alignments).positions
