####################################################################
# Helper Functions for Machine Learning
#
# author: Gregor Sturm (gregor.sturm@cs.tum.edu)
###################################################################

import math
import numpy as np

def seq2binvec(seq):
    """creates binary features such as A = [1,0,0,0] for a list of nucleotides"""
    nt_map = {"A": [1,0,0,0], "C": [0,1,0,0], "G": [0,0,1,0], "T": [0,0,0,1]}
    nt_map["U"] = nt_map["T"]
    nt_vec = []
    for nt in seq:
        nt_vec.extend(nt_map[nt])
    return nt_vec

def seq2nvec(seq):
    """creates numeric features such as A = 0, C = 1, ... for a list of nucleotides"""
    nt_map = {"A": 0, "C": 1, "G": 2, "T": 3}
    nt_map["U"] = nt_map["T"]
    nt_vec = []
    for nt in seq:
        nt_vec.append(nt_map[nt])
    return nt_vec

def normalize(array):
    """normalize a vector to values between 0 and 1"""
    amax = float(np.nanmax(array))
    amin = float(np.nanmin(array))
    if amax - amin == 0:
        return [0] * len(array)
    array = [(x-amin)/(amax-amin) for x in array]
    return array

def pseudocount(array):
    """ add pseudo-counts to array """
    array = [x+1 for x in array]
    return array

def log_scale(array):
    """ make log scale with pseudo-counts """
    array = pseudocount(array) #apply pseudo-counts
    array = [math.log(x) for x in array]
    return array

def normalize_between(array, lower, upper):
    """ normalize an array between a lower and upper bound"""
    scale = abs(lower-upper)
    normalized = normalize(array)
    scaled = [x * scale for x in normalized]
    scaled = [x + lower for x in scaled]
    return scaled

def smooth_forward(array, size):
    """ take the average of the next *size* positions
    for each value of the array
    """
    tmp_array = array + [array[-1]] * (size -1)
    array_sm = []
    for i in range(len(array)):
        array_sm.append(np.mean(tmp_array[i:i+size]))
    return array_sm
