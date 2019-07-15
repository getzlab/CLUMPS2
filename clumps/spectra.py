import os
import sys
from twobitreader import TwoBitFile
import pandas as pd
from collections import defaultdict
import numpy as np

compl = {'a':'t', 'c':'g', 'g':'c', 't':'a'}

mutation_indices = {}
i = 0
for b1 in 'acgt':
    for b3 in 'acgt':
        for b2 in 'ac':
            for n in 'acgt':
                if n == b2:
                    continue
                mutation_indices[(b1+b2+b3, n)] = i
                i += 1

coding_muts = set(["Translation_Start_Site", "Frame_Shift_Del", \
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

noncoding_muts = set(["3'UTR", "5'Flank", "3'Flank", "5'UTR", \
                      "De_novo_Start_OutOfFrame", "IGR", "Intron", \
                      "Non-coding_Transcript"])

def reverse_complement(abc):
    """
    Reverse codon.
    """
    compl = {'a':'t', 'c':'g', 'g':'c', 't':'a'}
    return compl[abc[2]] + compl[abc[1]] + compl[abc[0]]

def encode(abc, n):
    """
    Index mutation.
    """
    compl = {'a':'t', 'c':'g', 'g':'c', 't':'a'}
    if abc[1] == 'g' or abc[1] == 't':
        abc = reverse_complement(abc)
        n = compl[n]
    return mutation_indices[(abc,n)]

def calc_muts_spectra(input_maf, hgfile='./dat/hg19.2bit', out_file='sampleMutSpectra.txt'):
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

    hdr = list(mutation_indices.items())
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

    # Get Column Names
    df = pd.read_csv(input_maf, sep='\t').loc[:,['Tumor_Sample_Barcode','Variant_Classification','ttype']]

    # Compute frequency of tumor mutation
    for idx,row in df.iterrows():
        if row['Variant_Classification'] in coding_muts:
            pass
        elif row['Variant_Classification'] in noncoding_muts:
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
