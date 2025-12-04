# distutils: language = c++
# distutils: extra_compile_args = -std=c++20

from optimizer cimport SolveFailsafe

def solve(edges, rotations, alignments):
    return SolveFailsafe(edges, rotations, alignments).positions
