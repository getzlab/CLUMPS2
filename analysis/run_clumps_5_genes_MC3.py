import pandas as pd
pd.set_option('display.max_colwidth', None)
pd.set_option('display.max_rows', None)
import re
import subprocess
from clumps.wolF.tasks import clumps_workflow_localize_maf_only
import wolf
from datetime import datetime

maf = 'MC3-14k-merge_v1.5_genes.variant_classification.maf'
sampler = 'MutspecCoverageSampler'
run_name = 'run_clumps_5_genes_MC3_100000_perms'


wolf.utils.set_log_path('./%s_%s.log' % (run_name,datetime.today()))

with wolf.WolfWorkflow(
        workflow = clumps_workflow_localize_maf_only,#(maf,sampler),
        #run_name = run_name,
        conf = { "clust_frac": 0.5 }, # if you want to use more machines in the elastic cluster
        common_task_opts = { "retry" : 5 } # will retry every task up to 5 times
) as w:
    w.run(
        maf = maf,
        sampler = sampler,
        run_name = run_name,
        #genome_2bit = "gs://sa-clumps2-ref/dat/hg19.2bit",
        #fasta = "gs://sa-clumps2-ref/dat/UP000005640_9606.fasta.gz",
        #gpmaps = "gs://sa-clumps2-ref/dat/genomeProteomeMaps.txt",
        #prot2pdb_chunks = "/mnt/nfs/ro_disk/canine-c0820e9f17cde673ca995479dc35c9ae/prot2pdb_chunks/huniprot2pdb.run18_chunks/", #"gs://sa-clumps2-ref/dat/huniprot/huniprot2pdb.run18_chunks/",
        #pdb_dir = "gs://sa-clumps2-ref/dat/pdbs/ftp.wwpdb.org/pub/pdb/data/structures/divided/pdb",
        #coverage_track = "gs://sa-clumps2-ref/dat/cov/WEx_cov.fwb",
        #coverage_index = "gs://sa-clumps2-ref/dat/cov/WEx_cov.fwi",
        #cancer_genes = "gs://sa-clumps2-ref/dat/allCancerGenes.txt",
        #uniprot_map = "gs://sa-clumps2-ref/dat/huniprot/huniprot2pdb.run18.filt.txt",
        permutations = 100000,
        threads = 16,
        pancan_factor =.5,
        hillexp = 3,
        
        )
