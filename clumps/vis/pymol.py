import __main__
__main__.pymol_argv = [ 'pymol', '-qc']   ## Quiet and no GUI

import cgi
import pickle
import pymol
import sys
import pandas as pd

input_file = '/home/sanand/getzlab-CLUMPS2/outputs/exac/clumps_output_exac.tsv'

COL = {
    'bgcol':'white',
    'focusmappedcol':'gray70',
    'focusnonmappedcol':'paleyellow',
    'nonfocuscol':'aquamarine'
}

'''
python /src/createPymolSession.py 3j82
B:P60709 B all 0 0
/res/huniprot2pdb.run18.split/ /out splitByProtein/P60709 ACTB
'''

def render(gene, summary_file, pdb=None):
    """
    Render.

    Inputs:
        - gene: ENSEMBLE Gene Name
        - summary: CLUMPS output summary file
        - pdb: specific pdb-chain
    pdb provided if multiple structures picked per gene entry. If pdb is not
    provided, selects the most statistically signifcant gene.
    """
    df = pd.read_csv(summary_file, sep='\t')
    vals = df[df['GENE_NAMES']==gene]

    if pdb is not None:
        entry = vals[vals['PDBID-CHAIN']==pdb]
    else:
        entry = vals.iloc[0]

    pid,chain = entry['PDBID-CHAIN'].split('-')
    residues = entry['RESIDUES']

render('BPIFA2',input_file)
