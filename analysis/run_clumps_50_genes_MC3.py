import pandas as pd
pd.set_option('display.max_colwidth', None)
pd.set_option('display.max_rows', None)
import re
import subprocess
from clumps.wolF.tasks import clumps_workflow_localize_maf_only
import wolf
from datetime import datetime

maf = 'gs://sa-clumps2-ref/af_ajd/MC3_MAF/MC3-14k-merge_v1.50_genes.variant_classification.maf'
sampler = 'MutspecCoverageSampler'
run_name = 'run_clumps_50_genes_MC3_10000000_perms_v2_retry_recap'


wolf.utils.set_log_path('./%s_%s.log' % (run_name,datetime.today()))

with wolf.WolfWorkflow(
        workflow = clumps_workflow_localize_maf_only,#(maf,sampler),
        #run_name = run_name,
        conf = { "clust_frac": 0.99 }, # if you want to use more machines in the elastic cluster
        common_task_opts = { "retry" : 5 } # will retry every task up to 5 times
) as w:
    w.run(
        maf = maf,
        sampler = sampler,
        run_name = run_name,
        #genome_2bit = "gs://sa-clumps2-ref/dat/hg19.2bit",
        #fasta = "gs://sa-clumps2-ref/dat/UP000005640_9606.fasta.gz",
        #gpmaps = "gs://sa-clumps2-ref/dat/genomeProteomeMaps.txt",
        prot2pdb_chunks = 'gs://sa-clumps2-ref/af_ajd/uniprot50_genes/',
        #pdb_dir = "gs://sa-clumps2-ref/dat/pdbs/ftp.wwpdb.org/pub/pdb/data/structures/divided/pdb",
        #coverage_track = "gs://sa-clumps2-ref/dat/cov/WEx_cov.fwb",
        #coverage_index = "gs://sa-clumps2-ref/dat/cov/WEx_cov.fwi",
        #cancer_genes = "gs://sa-clumps2-ref/dat/allCancerGenes.txt",
        #uniprot_map = "gs://sa-clumps2-ref/dat/huniprot/huniprot2pdb.run18.filt.txt",
        permutations = 10000000,
        threads = 16,
        pancan_factor =.5,
        hillexp = 3,
        
        )
