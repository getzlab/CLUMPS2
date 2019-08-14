import random
import scipy as sp
import traceback
import sys

class AcetylSampler(object):
    """
    Samples only Lysine residues.
    --------------------------
    Takes in pdb residue id (pr) and a mapping default dictionary.
    Samples only from Lysines.
    """
    def __init__(self, res_id, mapping, res='LYS'):
        self.res_id = res_id
        self.res_id_idx = range(len(self.res_id))
        self.mapping = mapping
        self.lysines = list({k for k,v in self.mapping.items() if list(v)[0] ==res})
        self.lysines_idx = [idx for idx,x in enumerate(self.res_id) if x in self.lysines]

    def sample(self, mi):
        return (sorted(random.sample(self.lysines_idx, len(mi))),  sp.random.permutation(len(mi)))
