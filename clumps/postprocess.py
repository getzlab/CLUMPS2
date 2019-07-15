import os
import sys
import pandas as pd
from collections import defaultdict
from tqdm import tqdm

from statsmodels.sandbox.stats.multicomp import multipletests
import urllib.request as urllib2
import argparse

from utils import get_fragment_annot

def main():
    parser = argparse.ArgumentParser(description='Post-process CLUMPS.')

    parser.add_argument('-i', '--input_dir', required=True, help='<Required> Results directory from CLUMPS.')
    parser.add_argument('-d', '--proteins_dir', required=True, type=str, help='<Required> Split protein directory with mutation mapping information.')
    parser.add_argument('-p', '--uniprot_map', required=False, default='./dat/huniprot/huniprot2pdb.run18.filt.txt', help='Aggregate file of all uniprot mappings to PDB IDs with residue-level information.')
    parser.add_argument('-o', '--output_file', required=False, type=str, default='clumps_output.tsv', help='Output file from CLUMPS.')
    parser.add_argument('-c', '--cancer_genes', required=False, type=str, default='./dat/allCancerGenes.txt', help='List of cancer genes, tab-delimited.')
    parser.add_argument('--pdb_dir', required=False, default='./dat/pdbs/ftp.wwpdb.org/pub/pdb/data/structures/divided/pdb/', help='Directory of PDB files for parsing headers.')

    args = parser.parse_args()
    cancer_genes_df = pd.read_csv('./dat/allCancerGenes.txt',sep='\t')

    #----------------------------------------
    # Uniprot to Gene Mapping
    #----------------------------------------
    mapdi = {}

    with urllib2.urlopen('http://www.uniprot.org/uniprot/?query=reviewed:yes+AND+organism:9606&columns=id,entry%20name,reviewed,protein%20names,genes,organism,length,database(geneid)&format=tab') as f:
        for idx,line in tqdm(enumerate(f), desc='Mapping Uniprots'):
            if idx == 0:
                hdr = line.decode('utf-8').strip().split('\t')
                iEntry = hdr.index('Entry')
                iLength = hdr.index('Length')
                iGn = hdr.index('Gene names')
                iEntrez = hdr.index('Cross-reference (geneid)')
            else:
                line = line.decode('utf-8').strip('\n').split('\t')
                mapdi[line[iEntry]] = [map(lambda x:x.strip(';').strip('.'), line[iEntrez].split()), line[iGn], line[iLength]]

    #----------------------------------------
    # Cancer Proteins
    #----------------------------------------
    canproteins = dict()

    for u1 in mapdi:
        for e in mapdi[u1][0]:
            try:
                canproteins[u1] = cancer_genes_df[cancer_genes_df['EntrezGeneID']==int(e)].loc[:,'GeneSymbol'].values[0]
            except:
                pass

    #----------------------------------------
    # Parse Through Uniprot mappings to find proteins we have mutations for
    #----------------------------------------
    # Uniprot --> set({residue_id: (sample, ttype)})
    prot2muts = defaultdict(lambda: defaultdict(set))
    import traceback

    with open(args.uniprot_map, 'r') as f:
        for line in tqdm(f, desc='Parsing through UNIPROT mappings'):
            line = line.strip('\n').split('\t',4)
            # UNIPROT1, UNIPROT2, PDB-CH, STRAND, AA-MAP
            try:
                with open(os.path.join(args.proteins_dir, line[0]), 'r') as mut_file:
                    for mut_line in mut_file:
                        mut_line = mut_line.strip().split('\t')
                        # TUMOR_TYPE, SAMPLE, na, UNIPROT, MUT_SITE, SITE, MUT_TYPE

                        sample = (mut_line[1], mut_line[0])
                        prot2muts[line[0]][int(mut_line[5])].add(sample)
            except:
                pass

    #----------------------------------------
    # Parse Through UNIPROT Mappings to find coverage
    #----------------------------------------
    # (Uniprot1, Uniprot2, PDB-chain, residue_id) --> [start, end, pdb-identity, num_covered_mutated_residues]
    struct2startend = {}

    # (Uniprot1, Uniprot2, PDB-chain, residue_id) --> set([sample, ttype])
    struct2covsamples = {}

    with open(args.uniprot_map, 'r') as f:
        for line in tqdm(f, desc='Parsing through UNIPROT mappings'):
            line = line.strip('\n').split('\t',4)
            # UNIPROT1, UNIPROT2, PDB-CH, STRAND, AA-MAP

            if prot2muts[line[0]]:
                muts = prot2muts[line[0]]
                covmuts = set([])       ## covered mutation sites
                covsamples = set([])    ## samples contributing covered mutations

                for res in [int(i.split(':')[0]) for i in line[4].split()]:
                    # For each residue in a PDB structure
                    if res in list(muts.keys()):
                        # Add residues if mutations are found
                        covmuts.add(res)
                        covsamples.update(muts[res])

                        start_res = line[4].split(' ',1)[0].split(':')[0]
                        end_res = line[4].rsplit(' ',1)[-1].split(':')[0]
                        iden = 100.0

                        # Get percent identity
                        if line[3] != '-':
                            for i in line[3].split(' '):
                                if i.startswith('pdb_identity'):
                                    iden = float(i.split(':')[1])

                        struct2startend[(line[0],line[1],line[2],start_res)] = [int(start_res), int(end_res), iden, len(covmuts), covmuts]
                        struct2covsamples[(line[0],line[1],line[2],start_res)] = covsamples

    #----------------------------------------
    # Parse Through CLUMPS output
    #----------------------------------------
    # Uniprot1 --> [[Uniprot2, PDB-chain, start_res, end_res, pdb_ident, [pvals]], ...]
    # Mapping of uniprot ids to
    u1structs = {}

    for result_file in tqdm(os.listdir(args.input_dir), desc='Parsing CLUMPS output'):
        fileline,u1,u2,pdbch,rs = result_file.split('_')

        if u1 not in u1structs:
            u1structs[u1] = list()

        try:
            with open(os.path.join(args.input_dir, result_file)) as f:
                dat = f.readlines()

                if not dat or dat == ['#\n']:
                    # empty file: not tested dur to some filter
                    continue

                if dat[-1] != '#0\n':
                    print('WW: potentially insufficient number of simulations for %s' % fn)

                pval_dim = len(dat[0].split('\t'))
                P = [[0,0] for x in range(pval_dim)]

                for i in range(0,len(dat),2):
                    l = dat[i].strip().split('\t')
                    for j in range(pval_dim):
                        e,d = map(int, l[j].split('/'))  ## enumerator, denominator
                        P[j][0] += e
                        P[j][1] += d

                # Array of p-values
                P = [x[0]/float(x[1]) for x in P]
        except:
            print('EE', result_file)
            continue

        if (u1,u2,pdbch,rs) not in struct2startend:
            print("Not found in mapping: ", result_file)
            continue

        start_res,end_res,iden,covmuts,covmuts_sites = struct2startend[(u1,u2,pdbch,rs)]

        if covmuts < 3:
            # filter based on number of covered mutated residues
            continue

        covsamples = len(struct2covsamples[(u1,u2,pdbch,rs)])
        u1structs[u1].append([u2,pdbch,start_res,end_res,iden,P,covmuts,covsamples,covmuts_sites])

    #----------------------------------------
    # Filter output by definitino of cluster
    #----------------------------------------
    u1structs_filt = {}

    for u1 in tqdm(u1structs, desc='Filtering structures'):
        ss = u1structs[u1]
        ss.sort(key=lambda x:x[3]-x[2], reverse=True)  ## key is coverage length

        u1structs_filt[u1] = []
        clusters = []  ## clusters of homologous structures

        for i in ss:
            newcl = True  ## does this found a new cluster
            for cl in clusters:
                num_overlapping = 0
                # how many members of the cluster cl does this struct overlap with (at least 90% jaccard overlap)
                for j in cl:
                    ol = max(0, min(i[3], j[3]) - max(i[2], j[2]))
                    if float(ol)/((i[3]-i[2]) + (j[3]-j[2]) - ol) > 0.9:
                        num_overlapping += 1
                if num_overlapping/float(len(cl)) >= 0.5:
                    # overlaps with at least half of the members of the cluster by at least 90%
                    cl.append(i)
                    newcl = False
                    break
            if newcl:
                clusters.append([i])

        for cl in clusters:
            #cl.sort(key=lambda x:x[4]*(x[3]-x[2]), reverse=True)  ## sort by identity * mapped length
            cl.sort(key=lambda x:x[5][2])  ## sort by p-value (6A)
            #j = cl[int(round(len(cl)/2.0))-1] ## median p
            j = cl[0] ## min p
            ok = True
            for i in u1structs_filt[u1]:
                ol = max(0, min(i[3], j[3]) - max(i[2], j[2]))
                if float(ol)/min(i[3]-i[2], j[3]-j[2]) > 0.1:
                    ok = False
            if ok:
                u1structs_filt[u1].append(j)

    #----------------------------------------
    # Format output
    #----------------------------------------
    outdata = []
    min_mut_samples = 0

    for u1 in u1structs_filt:
        for l in u1structs_filt[u1]:
            cancerannot = ''

            if u1 in canproteins:
               cancerannot = canproteins[u1]

            pdb,ch = l[1].split('-')
            if u1 in mapdi:
                if mapdi[u1][1]:
                    gename = mapdi[u1][1].split()[0]
                else:
                    gename = ''
                protlen = mapdi[u1][2]
                present = 1
            else:
                gename = ''
                protlen = ''
                present = 0
            if present and not gename: # or genename.startswith('HLA-'):
                print("filtering out %s due to missing gene name (Ig?)" % u1)
                continue
            if l[7] < min_mut_samples:  ## filter on the number of samples
                continue

            sites = list(l[8])
            sites.sort()
            sites = [str(x) for x in sites]

            l = {'UNIPROT_ID':u1,
                 'GENE_NAMES':gename,
                 'IN_CANCER_GENE_LISTS':cancerannot,
                 'MAPPED_UNIPROT_ID':l[0],
                 'PDBID-CHAIN':l[1],
                 'PDB_FRAGMENT': get_fragment_annot(pdb, ch, args.pdb_dir),
                 'UNIPROT_SEQ_LENGTH':str(protlen),
                 'MAP_START':str(l[2]),
                 'MAP_END':str(l[3]),
                 'PERCENT_IDENTITY':'%.1f' % l[4],
                 'NSITES':str(l[6]),
                 'SITES':'+'.join(sites),
                 'NSAMPLES':str(l[7]),
                 'CLUMPS_P':l[5][2],
                 'CLUMPS_Q_FULL':-1,
                 'CLUMPS_Q_RESTRICTED':-1}
            outdata.append(l)

    # Multiple Hypothesis Testing
    pvals = [i['CLUMPS_P'] for i in outdata]
    qvals = multipletests(pvals, method='fdr_bh')[1]

    for i in range(len(outdata)):
        outdata[i]['CLUMPS_Q_FULL'] = qvals[i]

    ## just for cancer proteins
    pvals = [i['CLUMPS_P'] for i in outdata if i['IN_CANCER_GENE_LISTS']]
    idx = [i for i in range(len(outdata)) if outdata[i]['IN_CANCER_GENE_LISTS']]
    qvals = multipletests(pvals, method='fdr_bh')[1]

    for i in range(len(idx)):
        outdata[idx[i]]['CLUMPS_Q_RESTRICTED'] = qvals[i]

    # Sort
    outdata.sort(key=lambda x:x['CLUMPS_P'])

    output_df = pd.DataFrame.from_dict(outdata)
    output_df = output_df.loc[:,['GENE_NAMES','CLUMPS_P','CLUMPS_Q_FULL','NSAMPLES', 'NSITES','SITES','UNIPROT_ID','MAPPED_UNIPROT_ID','PDBID-CHAIN','PERCENT_IDENTITY','PDB_FRAGMENT','UNIPROT_SEQ_LENGTH','MAP_START','MAP_END','IN_CANCER_GENE_LISTS','CLUMPS_Q_RESTRICTED']]

    print("Saved output to {}".format(args.output_file))
    output_df.to_csv(args.output_file, sep='\t')

if __name__ == "__main__":
    main()
