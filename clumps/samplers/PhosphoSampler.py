import random
import scipy as sp
import traceback
import sys

class PhosphoSampler(object):
    """
    Samples only Phosphorylation residues.
    --------------------------
    Takes in pdb residue id (pr) and a mapping default dictionary.
    Samples only from Serines, Tyrosines, Threonines
    """
    def __init__(self, res_id, mapping):
        res = {"SER", "TYR", "THR"}

        self.res_id = res_id
        self.res_id_idx = range(len(self.res_id))
        self.mapping = mapping
        self.phospho_res = list({k for k,v in self.mapping.items() if list(v)[0] in res})
        self.phospho_res_idx = [idx for idx,x in enumerate(self.res_id) if x in self.phospho_res]

    def sample(self, mi):
        return (sorted(random.sample(self.phospho_res_idx, len(mi))),  sp.random.permutation(len(mi)))
