# cython: language_level=3
from libcpp.pair cimport pair
from libcpp.vector cimport vector

cdef extern from "svg_optimizer.cpp":
    cdef cppclass Solution:
        vector[int] positions
        int score
        vector[pair[int,int]] dangling_edges

    Solution SolveFailsafe(vector[int]& edges, vector[int]& rotations, vector[int]& alignments)
