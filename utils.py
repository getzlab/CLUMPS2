import ftplib
import os
import random
import re
import scipy as sp
import urllib
import numpy as np

from gzip import GzipFile
from lxml import etree
from prody import parsePDBStream, parsePDBHeader
from scipy.spatial.distance import euclidean
from collections import defaultdict

import contextlib
import tempfile
import subprocess

PDB_STRUCTURES = './dat/pdbs/ftp.wwpdb.org/pub/pdb/data/structures/divided/pdb'

aamap = {'ALA':'A','ARG':'R','ASN':'N','ASP':'D','CYS':'C','GLN':'Q','GLU':'E', \
         'GLY':'G','HIS':'H','ILE':'I','LEU':'L','LYS':'K','MET':'M','PHE':'F',\
         'PRO':'P','SER':'S','THR':'T','TRP':'W','TYR':'Y','VAL':'V','PYX':'C',\
         'SEP':'S','TPO':'T','TYS':'Y','MK8':'L','M3L':'K','DPR':'P','DSN':'S',\
         'ALN':'A','DLY':'K','MVA':'V','MLE':'L','DLE':'L','CSO':'C','PTR':'Y',\
         'BMT':'T','DAL':'A','FAK':'K','MSE':'M','SMF':'A','HYP':'P'}

def hill(t,n,x):
    return x**n/(t**n + x**n)

def parse_resmap(resmap):
    """
    Parse resmap.
    """
    ur = []  ## uniprot (u1) residues
    pr = []  ## pdb residues
    prd = {} ## pdb residues as dict (for containment checks)

    for i in resmap.split():
        x,y = i.split(':')
        ur.append(int(x))
        pr.append(int(y))
        prd[int(y)] = True

    return ur,pr,prd

@contextlib.contextmanager
def gunzipper(gz_file):
    """
    gunzipper.

    Takes in a file name.
    Returns a valid file object in an unzipped form.
    """
    with tempfile.NamedTemporaryFile('r', suffix=os.path.splitext(gz_file)[1]) as temp_file:
        subprocess.check_call("gzip -dc {} >> {}".format(gz_file, temp_file.name), executable='/bin/bash', shell=True)
        yield temp_file

def get_distance_matrix(pdbch, point='centroid', pdb_resids=None, return_centroid_coordinates=False):
    """
    Returns a distance matrix summarizing the Euclidean distances
    between residues in structure pdbch.
    """
    if return_centroid_coordinates and point != 'centroid':
        raise Exception('return_centroid_coordinates is True but point argument is not set to centroid')

    pdb_file = '%s/%s/pdb%s.ent.gz' % (PDB_STRUCTURES, pdbch[0][1:3], pdbch[0])

    # Open PDB File Using Prody
    with gunzipper(pdb_file) as pfile:
        if pdbch[1]:
            aa = parsePDBStream(pfile, chain=pdbch[1])  ## gets ALL atoms
        else:
            aa = parsePDBStream(pfile)  ## gets ALL atoms

    xx = aa.getResnums()
    yy = aa.getCoords()
    zz = aa.getResnames()

    if pdb_resids is None:
        pdb_resids = {}
        for i in range(len(xx)):
            if zz[i] in aamap:
                pdb_resids[xx[i]] = True
        pdb_resids = sorted(pdb_resids.keys())

    coords = {}
    for i in range(len(xx)):
        if xx[i] not in pdb_resids:
            continue
        if xx[i] not in coords:
            coords[xx[i]] = []
        coords[xx[i]].append(yy[i])  ## add coordinates of an atom belonging to this residue

    ## Euclidean distance matrix
    D = []
    for i in range(len(pdb_resids)):
        D.append(sp.zeros(i, dtype=sp.float32))

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
        return (D, pdb_resids, co)
    else:
        return (D, pdb_resids)

def transform_distance_matrix(D, ur, XPO):
    """
    Transform distance matrix.
    """
    ## transform distance matrix
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

def load_mut_freqs(in_file):
    """
    Load Mutational Frequencies

    Load mutational frequences for samples for allelic weighting.
    """
    mfreq = {}
    with open(in_file,'r') as f:
        for idx,line in enumerate(f):
            if idx > 0:
                line = line.strip().split('\t')
                mfreq[(line[0],line[1])] = float(line[4])
    return mfreq

def load_prot_file(protein_dir, uniprot):
    """
    Load Protein File

    Loads file from the splitProteinDir that has information about the protein
    from the input .maf.
    """
    if not os.path.isfile(os.path.join(protein_dir,uniprot)):
        print("file not found.")
    else:
        with open(os.path.join(protein_dir, uniprot)) as f:
            gm = f.read()

    return gm


def get_pdb_muts_overlap(residues, protein_mutations, hill_exp, USEPROVIDEDVALUES):
    """
    Get PDB / Mutations Overlap

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

            if USEPROVIDEDVALUES:
                mv.append(protein_mutations[residues[i]][0])
            else:
                mv.append(hill(2.0, hill_exp, protein_mutations[residues[i]][0]))

    return mi,mv,mt

def map_pos_with_weights(protein_dir, uniprot_id, mut_freqs, TTYPE, MUTTYPES, USEPROVIDEDVALUES, SAMPLEMUTFREQWEIGHT):
    """
    Map position with mutational weights.

    Returns protein_mutations:
        Residue --> [mutational_frequency, {tumor_type},{(sample,tumor_type)}]

    """
    protein_mutations = defaultdict(lambda: [0,set(),set()])

    gm = load_prot_file(protein_dir,uniprot_id)

    for mut_data in map(lambda x:x.split('\t'), gm.split('\n')):
        if mut_data == [''] or not mut_data[5]:
            continue

        if mut_data[6][0] in MUTTYPES:
            if TTYPE and not mut_data[0].startswith(TTYPE):
                continue

            pos = int(mut_data[5])

            if USEPROVIDEDVALUES:
                print("USEPROVIDEVALUES: if there are several lines \
                      for residue %d, taking avg. value." % p)

            if USEPROVIDEDVALUES:
                # Use provided data for mutation weights
                protein_mutations[pos][0] += float(mut_data[7])
            elif SAMPLEMUTFREQWEIGHT:
                # Sample-mutation frequency based on weighting of mutations
                protein_mutations[pos][0] += mut_freqs[(mut_data[0], mut_data[1])]
            else:
                # Weight all mutations equally (CLUMPS1)
                protein_mutations[pos][0] += 1

            protein_mutations[pos][1].add(mut_data[0])
            protein_mutations[pos][2].add((mut_data[1],mut_data[0]))

    if USEPROVIDEDVALUES:
        # average value
        for p in protein_mutations:
            if len(protein_mutations[p][2]) > 1:
                protein_mutations[p][0] /= len(protein_mutations[p][2])

    return protein_mutations
