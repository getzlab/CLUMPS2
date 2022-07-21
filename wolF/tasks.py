import wolf

CLUMPS_DOCKER_IMAGE = "gcr.io/broad-getzlab-workflows/clumps:v55"

class clumps_prep_task(wolf.Task):
    # Preparation for clumps input files:
    # computes mutational frequencies, spectra, and identifies protein structures needed
    resources = { "mem" : "8G" }
    
    # input data for the 'prep' step is the mutation annotation file (maf)
    # <Required> Input file for CLUMPS. Default expects .maf
    inputs = {
        "maf" : None,
        "genome_2bit" : None,
        "fasta" : None, 
        "gpmaps" : None
    }

    script = """
    mkdir clumps_preprocess
    clumps-prep --input ${maf} --output_dir clumps_preprocess --hgfile ${genome_2bit} --fasta ${fasta} --gpmaps ${gpmaps}
    """

    output_patterns = {
        "prep_outdir" : "clumps_preprocess"
    }
    
    docker = CLUMPS_DOCKER_IMAGE

class clumps_run_task(wolf.Task):
    # this task is the main clumps processing/algorithm
    resources = { "mem" : "8G" }
    
    # the input files for this step are the different individual prot2pdb chunks from the huniprot2pdb_chunks folder
    # provide a list of all the individual prot2pdb chunks (or the file path to each prot2pdb chunks file)
    
    # <Required> Directory of files titled with Uniprot IDs that have mutation information
    # <Required> File mapping uniprot ID to PDB ID with residue-level mapping information.
    # coverage_track is on the gs bucket
    inputs = {
        "clumps_preprocess" : None,
        "prot2pdb_chunks" : None,
        "pdb_dir" : None,
        "coverage_track" : None,
        "coverage_track_index" : None, # not actually used as an input; just needs to be localized alongside coverage_track
        "genome_2bit" : None,
        "fasta" : None,
        "gpmaps" : None,
        "sampler" : "UniformSampler",
        "max_perms" : 10000
    }

    overrides = { "prot2pdb_chunks" : "delayed" }

    ### bash script to run the `run-step` ###
    # un-tar the `clumps_preprocess.tar` file which is the output directory from the 'clumps-prep step'
    # this will create a diretory(folder) of the same name (clumps_preprocess)
    script = """    
    clumps --muts ${clumps_preprocess}/split_proteins \
        --maps ${prot2pdb_chunks} \
        --mut_freq ${clumps_preprocess}/mut_freq.txt \
        --out_dir clumps_results \
        --sampler ${sampler} \
        --coverage_track ${coverage_track} \
        --mut_spectra ${clumps_preprocess}/mut_spectra.txt \
        --pdb_dir ${pdb_dir} \
        --hgfile ${genome_2bit} --fasta ${fasta} --gpmaps ${gpmaps} \
        --max_rand ${max_perms} \
        --threads 8
    """

    resources = { "cpus-per-task" : 8 }

    output_patterns = {
        "run_outdir" : "clumps_results.tar"
    }
    
    docker = CLUMPS_DOCKER_IMAGE

class clumps_postprocess_task(wolf.Task):
    # Generates summary files from array outputs of clumps.
    # gather task which takes an array of inputs, but rather than dispatch an individual job for each input,
    # instead gather them together.
    resources = { "mem" : "8G" }
    
    # <Required> Results directory from CLUMPS
    # <Required> Split protein directory with mutation mapping information.
    inputs = {
        "clumps_preprocess" : None,
        "clumps_results" : None, #gather_param (paths to all input files being gathered)
        "cancer_genes" : None,
        "uniprot_map" : None,
        "pdb_dir" : None
    }
    
    ### bash script to run the `post-step` ###
    # untar the 'clumps_preprocess' to access the files
    # untar the 'clumps_results' to access the input files
    script = """
    #tar xf $clumps_preprocess
    #while read -r i; do
    #  tar xf $i
    #done < $clumps_results
    clumps-postprocess --input_dir ${clumps_results} \
      --proteins_dir ${clumps_preprocess}/split_proteins \
      --cancer_genes ${cancer_genes} \
      --uniprot_map ${uniprot_map} \
      --pdb_dir ${pdb_dir} \
      --output_file clumps_output.tsv
    """
    
    # Output file from CLUMPS with list of genes
    output_patterns = {
        "clumps_output" : "clumps_output.tsv"
    }
    
    # Docker Image
    docker = CLUMPS_DOCKER_IMAGE

###### workflow
def clumps_workflow(
  maf,
  genome_2bit = "gs://sa-clumps2-ref/dat/hg19.2bit",
  fasta = "gs://sa-clumps2-ref/dat/UP000005640_9606.fasta.gz",
  gpmaps = "gs://sa-clumps2-ref/dat/genomeProteomeMaps.txt",
  prot2pdb_chunks = "gs://sa-clumps2-ref/dat/huniprot/huniprot2pdb.run18_chunks/",
  pdb_dir = "gs://sa-clumps2-ref/dat/pdbs/ftp.wwpdb.org/pub/pdb/data/structures/divided/pdb",
  coverage_track = "gs://sa-clumps2-ref/dat/cov/WEx_cov.fwb",
  cancer_genes = "gs://sa-clumps2-ref/dat/allCancerGenes.txt",
  uniprot_map = "gs://sa-clumps2-ref/dat/huniprot/huniprot2pdb.run18.filt.txt"
):
    # (step-#0): localization task
    localization = wolf.LocalizeToDisk(
      files = {
        "maf" : maf,
        "genome_2bit" : genome_2bit,
        "fasta" : fasta,
        "gpmaps" : gpmaps,
        "prot2pdb_chunks" : prot2pdb_chunks,
        "pdb_dir" : pdb_dir,
        "coverage_track" : coverage_track,
        "cancer_genes" : cancer_genes,
        "uniprot_map" : uniprot_map
      }
    )

    clumps_prep = clumps_prep_task(
      inputs = {
        "maf" : localization["maf"],
        "genome_2bit" : localization["genome_2bit"],
        "fasta" : localization["fasta"],
        "gpmaps" : localization["gpmaps"]
      }
    )

#    # (step-#2): run task
#    clumps_run = clumps_run_task(
#      inputs = {
#        "clumps_preprocess" : clumps_prep["prep_outdir"],
#        "prot2pdb_chunks" : localization["prot2pdb_chunks"] chunks_list_18,
#        "pdb_dir" : localization["pdb_dir"],
#        "coverage_track" : localization["coverage_track"],
#        "genome_2bit" : localization["genome_2bit"],
#        "fasta" : localization["fasta"],
#        "gpmaps" : localization["gpmaps"]
#      }
#    )
