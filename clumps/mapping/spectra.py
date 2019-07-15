import os
from twobitreader import TwoBitFile
import pandas as pd
from collections import defaultdict
import numpy as np

from ..utils import reverse_complement, encode
from ..utils import CODING_MUTATIONS, NONCODING_MUTATIONS, MUTATION_INDICES

def calc_muts_spectra(input_maf, hgfile='../../dat/hg19.2bit', out_file='sampleMutSpectra.txt'):
    """
    Calculate mutational context.

    Finds each type of mutation and creates a categorical variable. Saves this output as
    the mutational spectra file.
    """
    hg = TwoBitFile(hgfile)
    sample_muts_context = defaultdict(lambda: [0]*96)

    df = pd.read_csv(input_maf, sep='\t').loc[:,['Tumor_Sample_Barcode','Variant_Type','Chromosome','Start_position','Reference_Allele','Tumor_Seq_Allele2','ttype']]
    df = df[df['Variant_Type']=='SNP']

    for idx,row in df.iterrows():
        ref_base = row['Reference_Allele'].lower()
        new_base = row['Tumor_Seq_Allele2'].lower()

        if ref_base == '-' or new_base == '-':
            continue
        if len(ref_base) > 1 or len(new_base) > 1:
            continue

        pos = int(row['Start_position'])
        chromosome = str(row['Chromosome'])

        if chromosome == '23':
            chromosome = 'X'
        elif chromosome == '24':
            chromosome = 'Y'
        elif chromosome == 'MT':
            chromosome = 'M'

        abc = hg['chr'+chromosome][pos-2:pos+1].lower()

        if abc[1] != ref_base and ref_base != '--':
            print(abc, ref_base, line)
            print('non-matching reference.')
            continue

        pat = (row['Tumor_Sample_Barcode'], row['ttype'])

        try:
            sample_muts_context[pat][encode(abc, new_base)] += 1
        except:
            ## because of Ns
            print("Because of Ns")
            continue

    hdr = list(MUTATION_INDICES.items())
    hdr.sort(key=lambda x:x[1])
    index_col = ['patient', 'ttype'] + [i[0][0]+'-'+i[0][1] for i in hdr]

    df = pd.DataFrame.from_dict(sample_muts_context).T
    df = df.reset_index()
    df.columns = index_col
    df.to_csv(out_file, sep='\t',index=None)

def calc_muts_freq(input_maf, out_file='sampleMutFreq.txt'):
    """
    Calculation mutation frequency of samples.

    For each sample, find the number of mutations, tumor type, and compute
    the rank score and zlog score.
    """
    print("SHOULD BE A MAF WITH ALL CODING MUTATIONS (NOT JUST MISSENSES)")

    # Get coding mutations
    di = defaultdict(lambda: defaultdict(lambda: 0))
    df = pd.read_csv(input_maf, sep='\t').loc[:,['Tumor_Sample_Barcode','Variant_Classification','ttype']]

    # Compute frequency of tumor mutation
    for idx,row in df.iterrows():
        if row['Variant_Classification'] in CODING_MUTATIONS:
            pass
        elif row['Variant_Classification'] in NONCODING_MUTATIONS:
            continue
        else:
            raise Exception("Mutation type {} unknown.".format(row['Variant_Classification']))
        di[row['ttype']][row['Tumor_Sample_Barcode']] += 1

    # Write out frequencies
    with open(out_file, 'w') as f:
        f.write('TTYPE\tSAMPLE\tMUT_COUNT\tTTYPE_RANK_SCORE\tZLOG_SCORE\n')
        print('Tumor Types: {}'.format(len(di)))

        for tumor_type in di:
            samples = list(di[tumor_type].items())
            print('\t{} | {} samples'.format(tumor_type,len(samples)))

            samples.sort(key=lambda x:x[1])

            logmf = [np.log10(x[1]) for x in samples]
            meanmf = np.mean(logmf)
            sdmf = np.std(logmf)

            for idx,sample in enumerate(samples):
                zlog = 1.0/max(1, (np.log10(sample[1]-meanmf) / sdmf))
                f.write('\t'.join([tumor_type, sample[0],  ## ttype and sample
                        '%d' % sample[1],  ## raw count
                        '%g' % (float(idx)/float(len(samples))),  ## rank score
                        '%g' % zlog]) + '\n')
