import random
import time
import argparse
import gzip
import os
import contextlib
import subprocess
from glob import glob
from multiprocessing import Process, Queue
import pkg_resources
import pandas as pd
import numpy as np
import math
from tqdm import tqdm
import sys

from canine import Orchestrator
from canine.utils import ArgumentHelper

from .samplers.UniformSampler import *
from .samplers.CoverageSampler import *
from .samplers.MutspecCoverageSampler import *
from .samplers.AcetylSampler import *
from .samplers.PhosphoSampler import *

from .mapping.mapper import GPmapper
from .utils import hill, parse_resmap, wap, fwap
from .utils import get_distance_matrix, transform_distance_matrix, get_pdb_muts_overlap, map_pos_with_weights, transform_distance_matrix2
from .utils import mkdir

def main():
    parser = argparse.ArgumentParser(description='Run CLUMPS.')
    parser.add_argument(
        '-d','--muts',
        required=True,
        type=str,
        help='<Required> Directory of files titled with Uniprot IDs that have mutation information'
    )
    parser.add_argument(
        '-m','--maps',
        required=True,
        type=str,
        help='<Required> File mapping uniprot ID to PDB ID with residue-level mapping information.'
    )
    parser.add_argument(
        '-f', '--mut_freq',
        required=False,
        help='Mutational frequenices of patient samples.',
        default=None
    )
    parser.add_argument(
        '--max_rand',
        type=int,
        default=10000000,
        help='Maximum number of random samples.'
    )
    parser.add_argument(
        '--mut_types',
        default=['M'],
        nargs='+',
        help='Mutation Types (N = nonsense, S = synonymous, M = mutation)',
    )
    parser.add_argument(
        '--sampler',
        default='UniformSampler',
        help='Sampler to use for Null Model.',
        choices=('UniformSampler','CoverageSampler','MutspecCoverageSampler', 'AcetylSampler', 'PhosphoSampler')
    )
    parser.add_argument(
        '--mut_spectra',
        required=False,
        help='Mutational spectra of patient samples.',
        default=None
    )
    parser.add_argument(
        '-e','--hill_exp',
        default=4,
        type=int,
        help='Hill Exponent'
    )
    parser.add_argument(
        '-p','--pancan_factor',
        default=1.0,
        type=float,
        help='Pan Cancer factor, value between 0 - 1. 1.0 if tumor type specified.'
    )
    parser.add_argument(
        '-t','--tumor_type',
        required=False,
        default='PanCan',
        type=str,
        help='Tumor type to run clumps on. PanCan indicates running clumps over all \
              tumor types found in input .maf.'
    )
    parser.add_argument(
        '-x', '--xpo',
        default=[3, 4.5, 6, 8, 10],
        type=list,
        help='Soft threshold parameter for truncated Gaussian.'
    )
    parser.add_argument(
        '--threads',
        type=int,
        default=1,
        help='Number of threads for sampling.'
    )
    parser.add_argument(
        '--n_jobs',
        type=int,
        default=1,
        help='Number of jobs to use.'
    )
    parser.add_argument(
        '--use_provided_values',
        required=False,
        default=None,
        help="Compute mutational frequencies with provided values."
    )
    parser.add_argument(
        '--coverage_track',
        required=False,
        default=None,
        help='Coverage track for null sampler.'
    )
    parser.add_argument(
        '-o', '--out_dir',
        required=False,
        default='./results',
        type=str,
        help='Output directory.'
    )
    parser.add_argument(
        '--pdb_dir',
        default=pkg_resources.resource_filename('clumps', './dat/pdbs/ftp.wwpdb.org/pub/pdb/data/structures/divided/pdb'),
        help='PDB directory.'
    )
    parser.add_argument(
        '--hgfile', type=str,
        default=os.path.join(pkg_resources.resource_filename('clumps', './dat'), 'hg19.2bit'),
        help='2bit human genome build file.'
    )
    parser.add_argument(
        '--fasta',
        type=str,
        default=os.path.join(pkg_resources.resource_filename('clumps', './dat'), 'UP000005640_9606.fasta.gz'),
        help='Protein primary sequence fasta file.'
    )
    parser.add_argument(
        '--gpmaps',
        type=str,
        default=os.path.join(pkg_resources.resource_filename('clumps', './dat'), 'genomeProteomeMaps.txt'),
        help='Genome Proteome Maps built from blast hits.'
    )
    parser.add_argument(
        '--slurm',
        action='store_true',
        help='Use SLURM backend.'
    )

    args = parser.parse_args()

    if args.tumor_type == 'PanCan':
        args.tumor_type = None

    if args.tumor_type and args.pancan_factor != 1.0:
        print('WARNING: args.pancan_factor is not 1 althought args.tumor_type is set. Correcting to args.pancan_factor=1', file = sys.stderr)
        args.pancan_factor = 1.0

    args.mut_types = set(args.mut_types)

    if args.sampler == 'MutspecCoverageSampler':
        assert args.mut_spectra is not None, "Provide mutational spectra data."
    if args.sampler == 'CoverageSampler' or args.sampler == 'MutspecCoverageSampler':
        assert args.coverage_track is not None, "Provide coverage track for null model."
    if args.sampler == 'AcetylSampler':
        assert 'A' in args.mut_types, "Specify only Acetylation events in model."
    if args.sampler == 'PhosphoSampler':
        assert 'P' in args.mut_types, "Specify only Phosphorylation events in model."

    #----------------------------------------
    # SLURM Bindings
    #----------------------------------------
    if args.slurm:
        # TODO: test this

        pipeline = {
            'name': 'wumps',
            'script': ["sudo docker run --rm -v $CANINE_ROOT:$CANINE_ROOT gcr.io/broad-cga-sanand-gtex/clumps clumps "],
            'inputs': {},
            'resources': {
                'cpus-per-task': 1,
                'mem-per-cpu': 900,
            },
            'backend': {
                'type': 'TransientGCP',
                'name': 'slorb-clumps',
                'compute_zone': 'us-east1-b',
                'controller_type': 'n1-standard-4',
                'worker_type': 'n1-highcpu-4',
                'compute_script': 'sudo docker pull gcr.io/broad-cga-sanand-gtex/clumps'

            },
            'outputs': {
                'clumps_output': args.out_dir+"/*"
            },
            'localization': {
                'transfer_bucket': 'sanand'
            }

        }

        for key,value in vars(args).items():
            if value is not None:
                if key == 'maps':
                    map_length = int(subprocess.run("zcat {} | wc -l".format(value), shell=True, executable='/bin/bash', stdout=subprocess.PIPE).stdout.decode())
                    chonk_size = math.ceil(map_length / args.n_jobs)

                    subprocess.check_call("zcat {} | split -d -l {} - huniprot_split_maps_ && gzip huniprot_split_maps_*".format(args.maps, chonk_size), shell=True, executable='/bin/bash')
                    pipeline['inputs']['maps'] = glob("huniprot_split_maps_*")

                elif key is not 'slurm':
                    pipeline['inputs'][key] = value if key != 'xpo' else str(value)
                    pipeline['script'][0] += "--{0} ${0} ".format(key)

        batch_id, jobs, outputs, sacct = Orchestrator(pipeline).run_pipeline()

        return

    #----------------------------------------
    # CLUMPS
    #----------------------------------------
    if args.sampler == 'CoverageSampler' or args.sampler == 'MutspecCoverageSampler':
        print("Building mapper...", file = sys.stderr)
        gpm = GPmapper(hgfile=args.hgfile, spfile=args.fasta, mapfile=args.gpmaps)

    # Load mutational frequencies
    if args.mut_freq is not None:
        mfreq = pd.read_csv(args.mut_freq, sep='\t').set_index(['TTYPE','SAMPLE']).loc[:,['ZLOG_SCORE']].to_dict()['ZLOG_SCORE']
    else:
        mfreq = {}

    mkdir(args.out_dir)

    #----------------------------------------
    # Run CLUMPS
    #----------------------------------------
    with contextlib.ExitStack() as stack:
        if args.sampler == 'CoverageSampler' or args.sampler == 'MutspecCoverageSampler':
            stack.enter_context(CoverageSampler.start_jvm())

        with gzip.open(args.maps, 'r') as f:
            for idx,line in tqdm(enumerate(f), desc=args.maps.split('/')[-1].split('.')[0]):
                u1,u2,pdbch,alist,resmap = line.decode('utf-8').strip('\n').split('\t', 4)

                if os.path.isfile(os.path.join(args.muts, u1)):
                    ####### Package into class #######
                    # TODO
                    #
                    pdbch = pdbch.split('-')
                    ur,pr,_ = parse_resmap(resmap)

                    if len(ur) < 5:
                        print("Bad mapping for {}.".format(ur), file = sys.stderr)
                        continue

                    # Skip structure if there are any negative UniProt -> PDB mappings
                    # (cause unknown, but likely an unusably bad structure)
                    if (pr < 0).any():
                        print(f"WARNING: skipping structure {u1} ({pdbch}) due to negative UniProt -> PDB mappings!", file = sys.stderr)
                        continue

                    # Remove non-unique UniProt -> PDB mappings (likely due to wonky homology modeling)
                    nuidx = np.flatnonzero(np.bincount(pr) > 1)
                    if len(nuidx):
                        rmidx = np.isin(pr, nuidx)
                        pr = pr[~rmidx]
                        ur = ur[~rmidx]
                        print(f"WARNING: removed {rmidx.sum()} residues with non-unique UniProt -> PDB mappings!", file = sys.stderr)

                    # Load Protein file
                    protein_muts = map_pos_with_weights(args.muts, u1, mfreq, args.tumor_type, args.mut_types, args.use_provided_values, args.mut_freq)

                    # Load PDB data
                    ## mi: index of mutated residue
                    ## mv: normalized mutation count at each residue
                    ## mt: cancer types contributing mutations
                    mi,mv,mt = get_pdb_muts_overlap(ur, protein_muts, args.hill_exp, args.use_provided_values)
                    mv = np.c_[mv]

                    # Load AA residue coordinates
                    if len(mi) > 0:
                        try:
                            D,x,pdb_resnames = get_distance_matrix(pdbch, args.pdb_dir, pdb_resids=pr)
                            #DDt = transform_distance_matrix(D, ur, args.xpo)
                            DDt2 = np.tril(transform_distance_matrix2(D, args.xpo), -1)
                        except:
                            print("Unable to load PDB...", file = sys.stderr)
                            continue

                        # print("Sampling {} | {} - {}".format(u1, pdbch, mi))

                        # Compute matrix
                        ## matrix that holds mv[i]*mv[j] values (sqrt or not)
                        #Mmv = []
                        #mvcorr = range(len(mv))

#                        for i in range(len(mi)):
#                            mrow = np.zeros(len(mi), np.float64)
#                            for j in range(len(mi)):
#                                #mrow[j] = np.sqrt(mv[i]*mv[j])  ## geometric mean; actually does not perform better in most cases
#                                if args.pancan_factor == 1.0:
#                                    mrow[j] = mv[i]*mv[j]
#                                else:
#                                    mrow[j] = (args.pancan_factor + (1.0-args.pancan_factor)*(len(mt[i] & mt[j])>0)) * mv[i]*mv[j]          ## product
#                            Mmv.append(mrow)

                        # Compute WAP score
                        #wap_obs = wap(mi, mvcorr, Mmv, DDt)
                        wap_obs = fwap(mi, mv, DDt2)

                        # Create Null Sampler
                        rnd = 0
                        P = [0]*len(args.xpo)
                        WAP_RND = [0]*len(args.xpo)
                        mireal = [i for i in mi]

                        # Sampler
                        try:
                            if args.sampler == 'UniformSampler':
                                sam = UniformSampler(ur)
                            elif args.sampler == 'CoverageSampler':
                                sam = CoverageSampler(ur, u1, args.coverage_track, gpm)
                            elif args.sampler == 'MutspecCoverageSampler':
                                sam = MutspecCoverageSampler(ur, u1, args.coverage_track, args.mut_spectra, gpm)
                                sam.calcMutSpecProbs(protein_muts)
                            elif args.sampler == 'AcetylSampler':
                                sam = AcetylSampler(pr, pdb_resnames)
                            elif args.sampler == 'PhosphoSampler':
                                sam = PhosphoSampler(pr, pdb_resnames)

                            # test sampler
                            _ = sam.sample(mireal)
                        except:
                            print("Error initializing {} for {} {} {}.".format(args.sampler, u1, u2, pdbch), file = sys.stderr)
                            continue

                        STARTTIME=time.time()

                        def rndThread(qu):
                            def booster():
                                """
                                Implements the booster algorithm by Getz et al. that saves CPU time.
                                Returns False if the randomization should be stopped.
                                """
                                ret = False
                                for i in range(len(args.xpo)):
                                    s = (rnd - p[i] + 1.0) / ((rnd + 3)*(p[i] + 1))
                                    if s >= 0.05271:  ## =0.9/(2*1.96)]**2:
                                        ret = True
                                        break
                                return ret

                            np.random.seed()
                            p = [0]*len(args.xpo)
                            wap_rnd = [0]*len(args.xpo)
                            rnd = 0

                            exitstatus=0  ## 0 means terminated OK, 1 means it had to abort due to timeout
                            while rnd < args.max_rand/args.threads and (rnd%1000 or booster()):  ## booster is applied once per 1000 randomizations
                                if not rnd%1000 and (time.time()-STARTTIME)/7200.0 > 1:
                                    exitstatus=1
                                    break

                                x = None
                                while x is None:
                                    ## some samplers will fail to yield a sample in some (small number of) of runs due to combinatorics
                                    x = sam.sample(mireal)

                                mi_perm, mut_perm_idx = x
                                #r = wap(mi, mvcorr, Mmv, DDt)
                                r = fwap(mi_perm, mv[mut_perm_idx], DDt2)

                                for rr in range(len(args.xpo)):
                                    wap_rnd[rr] += r[rr]
                                    if r[rr] >= wap_obs[rr]:
                                        p[rr] += 1
                                rnd += 1

                            qu.put((rnd,p,wap_rnd,exitstatus))

                        queue = Queue()
                        pcs = []
                        for r in range(args.threads):
                            x = Process(target=rndThread, args=(queue,))
                            x.start()
                            pcs.append(x)

                        for x in pcs:
                            x.join()

                        totalrnd = 0
                        totalexitstatus = 0
                        for r in range(args.threads):
                            rnd,p,wap_rnd,exitstatus = queue.get()
                            totalrnd += rnd
                            totalexitstatus += exitstatus
                            for i in range(len(args.xpo)):
                                P[i] += p[i]
                                WAP_RND[i] += wap_rnd[i]

                        with open(os.path.join(args.out_dir, '%s-%d_%s_%s_%s-%s_%s' % (args.maps.rsplit('.',1)[0].rsplit('_',1)[1], idx, u1, u2, pdbch[0], pdbch[1], resmap.split(':',1)[0])), 'a') as f:
                            f.write('\t'.join(['%d/%d' % (P[i], totalrnd) for i in range(len(P))]) + '\n')
                            f.write('#%d\n' % totalexitstatus)


if __name__ == "__main__":
	main()
