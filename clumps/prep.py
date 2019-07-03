import os
import sys

from tqdm import tqdm
import argparse
from Bio.Seq import Seq
import gzip
from twobitreader import TwoBitFile
import pandas as pd
from collections import defaultdict
import numpy as np

from mapper import GPmapper, make_gp_file, split_muts_file
from spectra import calc_muts_spectra, calc_muts_freq

def makedir(dirname):
    try:
        os.mkdir(dirname)
    except FileExistsError:
        print("{} already exists.".format(dirname))

def main():
    parser = argparse.ArgumentParser(description='Prepare inputs for CLUMPS.')
    parser.add_argument('-i','--input', required=True, type=str, help='<Required> Input file for CLUMPS. Default is to expect an input .maf file, but a processed acetylomics file, already in protein coordinates, is allowed using -a option.')
    parser.add_argument('-o', '--output_dir', type=str, required=False, help='Output directory for CLUMPS input files.')
    parser.add_argument('-a','--acetylomics', dest='acetyl_flag', action='store_true', help='Flag to indicate input is already in protein coordinates.')

    parser.add_argument('--hgfile', type=str, required=False, default='./dat/hg19.2bit', help='2bit human genome build file.')
    parser.add_argument('--fasta', type=str, required=False, default='./dat/UP000005640_9606.fasta.gz', help='Protein primary sequence fasta file.')
    parser.add_argument('--gpmaps', type=str, required=False, default='./dat/genomeProteomeMaps.txt', help='Genome Proteome Maps built from blast hits.')
    parser.add_argument('--ttype', type=str, required=False, default='BRCA', help='Tumor type for acetylomics data.')
    args = parser.parse_args()

    if args.output_dir is None:
        args.output_dir = args.input.split('.')[-2].split('/')[-1]

    makedir(args.output_dir)

    if args.acetyl_flag:
        """
        Input is a processed .csv file.

        -------------------------------------------
        Ex.:
            	patient	uniprot_id	site_position	accession_number	geneSymbol	variableSites	accession_numbers	value
            0	X11BR047	P02768	28	NP_000468.1	ALB	K28k	NP_000468.1	1
            1	X11BR047	P02768	36	NP_000468.1	ALB	K36k	NP_000468.1	1
            2	X11BR047	P02768	44	NP_000468.1	ALB	K44k	NP_000468.1	1
            3	X11BR047	P02768	65	NP_000468.1	ALB	K65k	NP_000468.1	0
            4	X11BR047	P02768	75	NP_000468.1	ALB	K75k	NP_000468.1	1
            ...

        """
        input_df = pd.read_csv(args.input, sep='\t')

        # Make frequency file in CLUMPS Format
        freq_df = input_df.groupby('patient').sum().sort_values('value').rename(columns={'value':'raw'}).drop(columns='site_position')
        freq_df['log10'] = np.log10(freq_df['raw'])

        meanmf = np.mean(freq_df['log10'])
        sdmf = np.std(freq_df['log10'])
        freq_df['zlog'] = 1.0 / (np.log10(freq_df['raw'] - meanmf) / sdmf)
        freq_df['rank_score'] = freq_df.reset_index().reset_index()['index'].astype(float).values / float(freq_df.shape[0])

        # Rename for formatting
        freq_df = freq_df.reset_index().rename(columns={'patient':'SAMPLE','raw':'MUT_COUNT', 'rank_score':'TTYPE_RANK_SCORE','zlog':'ZLOG_SCORE'})
        freq_df['TTYPE'] = args.ttype
        freq_df.loc[:,['TTYPE','SAMPLE','MUT_COUNT','TTYPE_RANK_SCORE','ZLOG_SCORE']].to_csv(os.path.join(args.output_dir,'mut_freq.txt'), sep='\t',index=None)

        # Make muts file in CLUMPS Format
        muts_df = input_df[input_df['value'] == 1].loc[:,['patient','uniprot_id','site_position','geneSymbol', 'variableSites']]
        muts_df['ttype'] = args.ttype
        muts_df['fill'] = 'na'
        muts_df['mtype'] = 'A'

        muts_df = muts_df.loc[:,['ttype','patient','fill','uniprot_id','variableSites','site_position','mtype']]
        muts_df.to_csv(os.path.join(args.output_dir,'muts.gp'), sep='\t')

        makedir(os.path.join(args.output_dir,'split_proteins'))

        # Create split protein Directory
        for uniprot in tqdm(set(muts_df['uniprot_id']), desc='Creating split protein directories.'):
            if uniprot is not np.nan:
                muts_df[muts_df['uniprot_id']==uniprot].to_csv(os.path.join(args.output_dir,'split_proteins',uniprot), sep='\t',header=None, index=None)

    else:
        """
        Input is a standard .maf file. Run through steps to map from genomic
        to protein coordinates.
        """
        # Create genome-proteome mapper object
        gpm = GPmapper(hgfile=args.hgfile, spfile=args.fasta, mapfile=args.gpmaps)

        # Make genome-proteome mutations file
        # This file maps each patient, tumor-type, chormosome, position
        # to uniprot ID and amino acid location of change
        make_gp_file(gpm, args.input, output_file=os.path.join(args.output_dir,'muts.gp'))

        # Split genome-proteome mutations file
        # Splits each protein-mapped-mutation into individual files
        split_muts_file(os.path.join(args.output_dir,'muts.gp'), split_protein_dir=os.path.join(args.output_dir, 'split_proteins'))

        # Mutational Spectra
        # Compute mtuational landscape by profiling mutation types
        calc_muts_spectra(args.input, hgfile=args.hgfile, out_file=os.path.join(args.output_dir,'mut_spectra.txt'))

        # Mutational Frequencies
        # Find number of mutations, tumor-type, rank score, zlog score
        calc_muts_freq(args.input, out_file=os.path.join(args.output_dir,'mut_freq.txt'))




if __name__ == "__main__":
	main()
