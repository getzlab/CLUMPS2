import ftplib
import os
import random
import re
import scipy as sp
import numpy as np
from gzip import GzipFile
from lxml import etree
from prody import parsePDBStream, parsePDBHeader
from scipy.spatial.distance import euclidean
from collections import defaultdict

from prody import confProDy

import contextlib
import tempfile
import subprocess

AMINO_ACID_MAP = {'ALA':'A','ARG':'R','ASN':'N','ASP':'D','CYS':'C','GLN':'Q','GLU':'E', \
         'GLY':'G','HIS':'H','ILE':'I','LEU':'L','LYS':'K','MET':'M','PHE':'F',\
         'PRO':'P','SER':'S','THR':'T','TRP':'W','TYR':'Y','VAL':'V','PYX':'C',\
         'SEP':'S','TPO':'T','TYS':'Y','MK8':'L','M3L':'K','DPR':'P','DSN':'S',\
         'ALN':'A','DLY':'K','MVA':'V','MLE':'L','DLE':'L','CSO':'C','PTR':'Y',\
         'BMT':'T','DAL':'A','FAK':'K','MSE':'M','SMF':'A','HYP':'P'}

BASE_COMPLEMENT = {'a':'t', 'c':'g', 'g':'c', 't':'a'}

MUTATION_INDICES = {}
i = 0
for b1 in 'acgt':
    for b3 in 'acgt':
        for b2 in 'ac':
            for n in 'acgt':
                if n == b2:
                    continue
                MUTATION_INDICES[(b1+b2+b3, n)] = i
                i += 1

CODING_MUTATIONS = set(["Translation_Start_Site", "Frame_Shift_Del", \
                   "Frame_Shift_Ins", "In_Frame_Del", "In_frame_Del", \
                   "In_Frame_Ins", "Missense_Mutation", "Missense", \
                   "Nonsense_Mutation", "Nonsense", "Nonstop_Mutation", \
                   "Silent", "Synonymous", "Splice_Site", "Start_Codon_DNP", \
                   "Start_Codon_Del", "Start_Codon_Ins", "Start_Codon_ONP", \
                   "Stop_Codon_DNP", "Stop_Codon_Del", "Stop_Codon_Ins", "RNA", \
                   "lincRNA", "De_novo_Start_InFrame", "Start_Codon_SNP", \
                   "Read-through", "Splice_Region", "Splice_Site_DNP", \
                   "Splice_Site_Del", "Splice_Site_Ins", "Splice_Site_ONP", \
                   "Splice_Site_SNP", "Splice_site", "Splice_site_SNP"])

NONCODING_MUTATIONS = set(["3'UTR", "5'Flank", "3'Flank", "5'UTR", \
                      "De_novo_Start_OutOfFrame", "IGR", "Intron", \
                      "Non-coding_Transcript"])

def reverse_complement(abc):
    """
    Reverse codon.
    """
    return BASE_COMPLEMENT[abc[2]] + BASE_COMPLEMENT[abc[1]] + BASE_COMPLEMENT[abc[0]]

def encode(abc, n):
    """
    Index mutation.
    """
    if abc[1] == 'g' or abc[1] == 't':
        abc = reverse_complement(abc)
        n = BASE_COMPLEMENT[n]
    return MUTATION_INDICES[(abc,n)]

def mkdir(path):
    try:
        os.makedirs(path)
    except FileExistsError:
        pass

def hill(t,n,x):
    """
    Hill function.
    """
    return x**n/(t**n + x**n)

def parse_resmap(resmap):
    """
    Parse resmap.
    --------------------------
    Parse residue map.
    """
    ur = []  ## uniprot (u1) residues
    pr = []  ## pdb residues
    prd = {} ## pdb residues as dict (for containment checks)

    for i in resmap.split():
        x,y = i.split(':')
        ur.append(int(x))
        pr.append(int(y))
        prd[int(y)] = True

    return np.r_[ur],np.r_[pr],prd

@contextlib.contextmanager
def gunzipper(gz_file):
    """
    Gunzipper
    --------------------------
    Takes in a file name. Returns a valid file object in an unzipped form. Uses
    context manager to create a temporary file and save into that format. This
    was written for programmatic extraction of PDB files using Prody.

    Inputs:
        gz_file: gzipped file

    Outputs:
        temp_file: temporary file
    """
    with tempfile.NamedTemporaryFile('r', suffix=os.path.splitext(gz_file)[1]) as temp_file:
        subprocess.check_call("gzip -dc {} >> {}".format(gz_file, temp_file.name), executable='/bin/bash', shell=True)
        yield temp_file

confProDy(verbosity='none')

def get_distance_matrix(pdbch, pdb_structures_dir, point='centroid', pdb_resids=None, return_centroid_coordinates=False):
    """
    Get Distance Matrix
    --------------------------
    Returns a distance matrix summarizing the Euclidean distances
    between residues in structure pdbch.
    """
    if return_centroid_coordinates and point != 'centroid':
        raise Exception('return_centroid_coordinates is True but point argument is not set to centroid')

    pdb_file = '%s/%s/pdb%s.ent.gz' % (pdb_structures_dir, pdbch[0][1:3], pdbch[0])

    # Open PDB File Using Prody
    with gunzipper(pdb_file) as pfile:
        if pdbch[1]:
            aa = parsePDBStream(pfile, chain=pdbch[1])  ## gets ALL atoms
        else:
            aa = parsePDBStream(pfile)  ## gets ALL atoms

    xx = aa.getResnums()
    yy = aa.getCoords()
    zz = aa.getResnames()

    # if list of PDB residues is not provided, look them up
    if pdb_resids is None:
        pdb_resids = {}
        for i in range(len(xx)):
            if zz[i] in AMINO_ACID_MAP:
                pdb_resids[xx[i]] = True
        pdb_resids = sorted(pdb_resids.keys())

    # otherwise, perform sanity check that provided residue list comprises valid amino acids
    else:
        if len(set(pdb_resids) - set(xx[np.r_[[z in AMINO_ACID_MAP for z in zz]]])):
            raise ValueError("Invalid PDB residues specified!")

    mapped_pdb_to_aa = defaultdict(set)
    for idx,resnum in enumerate(xx):
        if resnum in pdb_resids:
            mapped_pdb_to_aa[resnum].add(zz[idx])

    coords = {}
    for i in range(len(xx)):
        if xx[i] not in pdb_resids:
            continue
        if xx[i] not in coords:
            coords[xx[i]] = []
        coords[xx[i]].append(yy[i])  ## add coordinates of an atom belonging to this residue

    ## Euclidean distance matrix
    D = np.zeros(len(pdb_resids)*np.r_[1, 1])

    if point == 'centroid':
        ## distance between centroids
        ## calculate residue centroid positions
        centroids = {}
        for k in coords:
            centroids[k] = np.mean(np.array(coords[k]), 0)

        co = [centroids[i] for i in pdb_resids]  ## pdb residue coordinates

        for i in range(len(pdb_resids)):
            for j in range(i):
                D[i][j] = euclidean(co[i], co[j])

    elif point == 'min':
        ## min-distance (atom pairs)
        co = [coords[i] for i in pdb_resids]  ## pdb atom coordinates
        for i in range(len(pdb_resids)):
            for j in range(i):
                m = 10000000
                for x in co[i]:
                    for y in co[j]:
                        e = euclidean(x, y)
                        if e < m:
                            m = e
                D[i][j] = m
    else:
        raise Exception('Unknown setting for point: %s' % point)

    if return_centroid_coordinates:
        return (D, pdb_resids, mapped_pdb_to_aa, co)
    else:
        return (D, pdb_resids, mapped_pdb_to_aa)

def transform_distance_matrix(D, ur, XPO):
    """
    Transform distance matrix.
    --------------------------
    Transforms distance matrix.
    """
    DDt = []  ## array of transformed distance matrices
    for soft_thresh_idx in range(len(XPO)):
        den = 2.0 * XPO[soft_thresh_idx]**2
        m = []
        for i in range(len(ur)):
            mrow = sp.zeros(i, dtype=sp.float32)
            for j in range(i):
                mrow[j] = sp.exp(-(D[i][j]**2)/den)
            m.append(mrow)
        DDt.append(m)

    return DDt

def transform_distance_matrix2(D, XPO):
    """
    Transform distance matrix.
    --------------------------
    Transforms distance matrix.
    """
    DDt = []  ## array of transformed distance matrices
    for soft_thresh_idx in range(len(XPO)):
        den = 2.0 * XPO[soft_thresh_idx]**2
        DDt.append(np.exp(-D**2/den))

    return DDt

def load_prot_file(protein_dir, uniprot):
    """
    Load Protein File
    --------------------------
    Loads file from the split_proteins directory that has information
    about the protein from the input .maf.
    """
    try:
        with open(os.path.join(protein_dir,uniprot)) as f:
            return f.read()
    except:
        print("{} NOT FOUND.".format(os.path.join(protein_dir,uniprot)))

def get_pdb_muts_overlap(residues, protein_mutations, hill_exp, use_provided_values):
    """
    Get PDB & Mutations Overlap
    --------------------------
    Get indices of mutated residuces, the normalized mutation count
    of each residue, and the cancer types contributing to the mutations.
    """
    mi = []  ## index of mutated residue
    mv = []  ## normalized mutation count at each residue
    mt = []  ## cancer types contributing mutations

    for i in range(len(residues)):
        if residues[i] in protein_mutations:
            mi.append(i)
            mt.append(protein_mutations[residues[i]][1])

            if use_provided_values:
                mv.append(protein_mutations[residues[i]][0])
            else:
                mv.append(hill(2.0, hill_exp, protein_mutations[residues[i]][0]))

    return mi,mv,mt

def map_pos_with_weights(protein_dir, uniprot_id, mut_freqs, ttype, mutation_types, use_provided_values, sample_freq_weight):
    """
    Map position with mutational weights.
    --------------------------
    Returns protein_mutations:
        Residue --> [mutational_frequency, {tumor_type},{(sample,tumor_type)}]
    """
    protein_mutations = defaultdict(lambda: [0,set(),set()])

    gm = load_prot_file(protein_dir,uniprot_id)

    for mut_data in map(lambda x:x.split('\t'), gm.split('\n')):
        if mut_data == [''] or not mut_data[5]:
            continue

        if mut_data[6][0] in mutation_types:
            if ttype and not mut_data[0].startswith(ttype):
                continue

            pos = int(mut_data[5])

            if use_provided_values:
                print("USEPROVIDEVALUES: if there are several lines \
                      for residue %d, taking avg. value." % p)

            if use_provided_values:
                # Use provided data for mutation weights
                protein_mutations[pos][0] += float(mut_data[7])
            elif sample_freq_weight:
                # Sample-mutation frequency based on weighting of mutations
                protein_mutations[pos][0] += mut_freqs[(mut_data[0], mut_data[1])]
            else:
                # Weight all mutations equally (CLUMPS1)
                protein_mutations[pos][0] += 1

            protein_mutations[pos][1].add(mut_data[0])
            protein_mutations[pos][2].add((mut_data[1],mut_data[0]))

    if use_provided_values:
        # average value
        for p in protein_mutations:
            if len(protein_mutations[p][2]) > 1:
                protein_mutations[p][0] /= len(protein_mutations[p][2])

    return protein_mutations

def wap(mut_indices, mvcorr, Mmv, DDt):
    """
    WAP Score
    --------------------------
    Compute WAP score to summarize pairwise distances between
    mutated residues in 3D protein structure.
    """
    s = sp.zeros(len(DDt), sp.float64)
    for mat in range(len(DDt)):
        d = DDt[mat]
        for i in range(len(mut_indices)):
            dcol = d[mut_indices[i]]
            for j in range(i):
                s[mat] += Mmv[mvcorr[i]][mvcorr[j]] * dcol[mut_indices[j]]
    
    return s

def fwap(mi, mv, DDt):
    scores = np.zeros(len(DDt))
    for xpo_idx in range(len(DDt)):
        scores[xpo_idx] = mv.T@DDt[xpo_idx][mi, :][:, mi]@mv
    return scores

def get_fragment_annot(pdb, ch, pdb_dir):
    """
    Get pdb-fragment annotation.
    --------------------------
    Returns the PDB annotation for a pdb-chain.
    """
    poly = parsePDBHeader(os.path.join(pdb_dir, pdb[1:3], "pdb{}.ent.gz".format(pdb)), 'polymers')
    for p in poly:
        if p.chid != ch:
            continue
        return p.fragment
    return 'NA'
