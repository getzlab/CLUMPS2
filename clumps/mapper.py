import os
import sys
from Bio.Seq import Seq
import gzip
from twobitreader import TwoBitFile
import pandas as pd
from collections import defaultdict

class GPmapper(object):
    def __init__(self, hgfile='./dat/hg19.2bit',
                 spfile='./dat/UP000005640_9606.fasta.gz',
                 mapfile='./dat/genomeProteomeMaps.txt'):
        """
        GPmapper.
        """
        self.hg = TwoBitFile(hgfile)
        self.sp = {}
        self.gen2prot = defaultdict(list)
        self.prot2gen = defaultdict(list)

        # Grab Uniprot IDs and Map to Amino Acid Sequence
        with gzip.open(spfile,'r') as f:
            for line in f:
                line = line.decode('utf-8')
                if line.startswith('>'):
                    try:
                        self.sp[uniprot_id] = aa_sequence
                    except:
                        pass
                    uniprot_id = line.strip().lstrip('>').split('|')[1]
                    aa_sequence = ''
                else:
                    aa_sequence += line.strip()

        # Parse through genome protein map file
        with open(mapfile, 'r') as f:
            for line in f:
                line = line.strip('\n').split('\t')

                # line[0]: chromsome
                # line[1]: start idx
                # line[2]: end idx
                # line[3]: exons
                # line[4]: strand
                # line[5]: uniprot
                # line[6]: idk
                # line[7]: idk

                bl = set()
                if line[7]:
                    bl = set(map(int, line[7].split(',')))

                # gen2prot: genome to proteome mapping
                # chr --> [start idx, end idx, [exon ranges], strand, uniprot, [IDK], set(IDK)]
                self.gen2prot[line[0]].append([
                    int(line[1]),
                    int(line[2]),
                    [(int(i[0]),int(i[1])) for i in map(lambda x:x.split('-'), line[3].split(','))],
                    line[4],
                    line[5],
                    [((int(i[0][:-2]), int(i[0][-1])), (int(i[1][:-2]), int(i[1][-1]))) for i in map(lambda x:x.split('-'), line[6].split(','))],
                    bl
                ])

                # prot2gen: protein to genome mapping
                # Uniprot --> [Chr, index]
                self.prot2gen[line[5]].append((line[0], len(self.gen2prot[line[0]])-1))

    def _chrom_map(self, c):
        """
        Map chromsome name.
        """
        if c == '23':
            c = 'X'
        elif c == '24':
            c = 'Y'
        elif c == 'MT':
            c = 'M'
        c = 'chr'+c
        return c

    def find_exons(self, chromosome, pos):
        """
        find exons

        Find exon and gene IDs that map to a given chromosome
        and position. Returns them as a list of each gene
        ID and exon ID.
        """
        ret = []
        for gi in range(len(self.gen2prot[chromosome])):
            g = self.gen2prot[chromosome][gi]

            if not g[0] <= pos < g[1]:
                # Not contained in gene range
                continue

            if pos in g[6]:
                # Mystery set
                continue

            for ei in range(len(g[2])):
                if not g[2][ei][0] <= pos < g[2][ei][1]:
                    # Not found in this exon
                    continue
                ret.append((gi,ei))
                break
        return ret

    def map_gen2prot(self, chromosome, pos, gi, ei, newbase):
            """
            map gen2prot

            Maps a genomic change to change(s) in protein(s).
            Either position or both gi, ei should be None.

                - chromosome
                - pos: position in chromosome
                - gi: gene index
                - ei: exon index
                - newbase: new-base change


            g = self.gen2prot[chromosome][gi]
                g[0] = start_idx
                g[1] = end_idx
                g[2] = exon ranges
                g[3] = strand
                g[4] = uniprot_id
                g[5] = ?
                g[6] = ?

            """
            pos -= 1
            chromosome = self._chrom_map(chromosome)

            if gi is None:
                exons = self.find_exons(chromosome, pos)
            else:
                exons = [(gi,ei)]

            ret = []

            for gi,ei in exons:
                g = self.gen2prot[chromosome][gi]

                if g[4] not in self.sp:
                    # Depcreated UniProt ID
                    continue
                if g[3] == '+':
                    # Flip around to account for stranding
                    offset = pos - g[2][ei][0]
                else:
                    offset = g[2][ei][1] - pos - 1

                # Offset in protein coordinates
                poffset = (offset//3, offset%3)
                qblockstart = g[5][ei][0]
                pp = (
                    qblockstart[0] + poffset[0] + (qblockstart[1] + poffset[1])//3,
                    (qblockstart[1] + poffset[1])%3
                )

                # residual
                if g[3] == '+':
                    r = pp[1]      ## residual
                else:
                    r = 2-pp[1]    ## residual

                codon = []
                fromPrevExon = min(0, pos-r - g[2][ei][0])

                if fromPrevExon:
                    codon += range(g[2][ei-1][1] + fromPrevExon, g[2][ei-1][1])

                fromCurrExon = [max(pos-r, g[2][ei][0]), min(pos+(3-r), g[2][ei][1])]
                codon += range(fromCurrExon[0], fromCurrExon[1])

                fromNextExon = max(0, pos+(3-r) - g[2][ei][1])

                if fromNextExon:
                    codon += range(g[2][ei+1][0], g[2][ei+1][0]+fromNextExon)

                codonseq = ''.join([self.hg[chromosome][i:i+1] for i in codon])  ## (reference) codon sequence
                ncodonseq = list(codonseq)                                       ## new (mutant) codon sequence

                if newbase:
                    ncodonseq[r] = newbase
                ncodonseq = ''.join(ncodonseq)

                if g[3] == '+':
                    try:
                        # amino acid according to uniprot
                        raa = self.sp[g[4]][pp[0]]
                    except:
                        print("ONO", pp[0], pp[1])
                        continue  ## index error (seq has changed)
                    taa = str(Seq(codonseq).translate())  ## amino acid resulting from translation
                    naa = str(Seq(ncodonseq).translate()) ## new (mutant) amino acid
                else:
                    try:
                        raa = self.sp[g[4]][pp[0]]            ## actual amino acid according to uniprot
                    except:
                        continue  ## index error (seq has changed)
                    taa = str(Seq(codonseq).reverse_complement().translate())  ## amino acid resulting from translation
                    naa = str(Seq(ncodonseq).reverse_complement().translate()) ## new (mutant) amino acid

                if taa != raa:
                    print('WARNING: non-matching reference and translated AA', chromosome, pos, fromPrevExon, fromCurrExon, fromNextExon, len(codon), ei, len(g[2]))
                    continue

                ret.append((g[4], 1+pp[0], taa, naa))
            return ret

def make_gp_file(gpm, input_maf, output_file='muts.gp'):
    """
    Make GP File

    Create a mapping of mutations based on protein coordinates.
    Takes in an input maf with the mutations of interest and a GPmapper object.
    """
    with open(output_file, 'w') as mfile:
        with open(input_maf, 'r') as f:
            h = f.readline().strip().split('\t')

            iChr = h.index('Chromosome')
            iPos = h.index('Start_position')
            iRefAllele = h.index('Reference_Allele')
            iTumAllele2 = h.index('Tumor_Seq_Allele2')
            iPatient = h.index('Tumor_Sample_Barcode')
            iTtype = h.index('ttype')

            mfile.write('\t'.join(['patient', 'ttype', 'chr', 'pos', 'refbase', 'newbase', 'uniprot_change']) + '\n')
            COUNTER = 0

            for idx,line in enumerate(f):
                if idx > 0:
                    maf_line = line.strip('\n').split('\t')

                    if len(maf_line[iTumAllele2]) > 1 or maf_line[iTumAllele2] == '-':
                        continue
                    if len(maf_line[iRefAllele]) > 1 or maf_line[iRefAllele] == '-':
                        continue

                    # newly mutated base
                    newbase = maf_line[iTumAllele2]

                    # protien:AA-site-newAA
                    protein_aa_change = gpm.map_gen2prot(maf_line[iChr], int(maf_line[iPos]), None, None, newbase)
                    if not protein_aa_change:
                        continue

                    protein_aa_change = '; '.join(map(lambda x:'%s:%s%d%s' % (x[0], x[2], x[1], x[3]), protein_aa_change))

                    # Join all
                    nl = [maf_line[iPatient], maf_line[iTtype], maf_line[iChr], maf_line[iPos], maf_line[iRefAllele], newbase, protein_aa_change]

                    mfile.write('\t'.join(nl) + '\n')
        print(COUNTER)

def split_muts_file(input_gp, split_protein_dir='splitByProtein', mut_types=set(['M'])):
    """
    Split Muts File

    Take in mapped mutation file and splits each of these into individual protein files.
    Only writes out mutation types contained in the mut_types set.
        - N: nonsense
        - S: synonymous
        - M: missense (defualt)
    """
    pdata = defaultdict(list)

    with open(input_gp, 'r') as f:
        hdr = f.readline().strip('\n').split('\t')
        iUniprotChange = hdr.index('uniprot_change')
        iPatient = hdr.index('patient')
        iTtype = hdr.index('ttype')

        for idx,line in enumerate(f):
            if idx > 0:
                line = line.strip('\n').split('\t')

                if not line[iUniprotChange]:
                    continue

                # List all Protein Mutations
                # Determines the mutatation type
                for um in line[iUniprotChange].split('; '):
                    u,m = um.split(':')
                    if m[-1] == '*':
                        mt = 'N'
                    elif m[-1] == m[0]:
                        mt = 'S'
                    else:
                        mt = 'M'

                    # Only consider mutations in the mut_types set
                    if mt not in mut_types:
                        continue

                    aa_mutation_index = m[1:-1]
                    while not aa_mutation_index[-1].isdigit():
                        # Some mutations are annotated as SNPs although they are not
                        aa_mutation_index = aa_mutation_index[:-1]

                    # Join data
                    dline = '\t'.join([line[iTtype], line[iPatient], 'na', u, 'p.' + m, aa_mutation_index, mt])
                    pdata[u].append(dline)

        try:
            os.mkdir(split_protein_dir)
        except:
            print('dir exists.')

        for u in pdata:
            with open(os.path.join(split_protein_dir, u), 'w') as f:
                f.write('\n'.join(pdata[u]) + '\n')
                f.close()
