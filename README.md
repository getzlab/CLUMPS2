# CLUMPS

Refactored code for CLUMPS written in Python3.

This repository is based on the original code written for CLUMPS. Please cite this paper if using this for your research:

> Atanas Kamburov, Michael S. Lawrence, Paz Polak, Ignaty Leshchiner, Kasper Lage, Todd R. Golub, Eric S. Lander, and Gad Getz. "Comprehensive assessment of cancer missense mutation clustering in protein structures" PNAS October 6, 2015 112 (40) E5486-E5495; first published September 21, 2015 https://doi.org/10.1073/pnas.1516373112

---

### Installation

To install, run the following in the current directory (pypi support in progress):

```
pip install -e .
```

### Usage

CLUMPS consists of three main modules:

#### Prep
Preparation for clumps input files: computes mutational frequencies, spectra, and identifies protein structures needed.

```
clumps-prep [-h] -i INPUT [-o OUTPUT_DIR] [-a] [--hgfile HGFILE]
     [--fasta FASTA] [--gpmaps GPMAPS]

Prepare inputs for CLUMPS.

optional arguments:
  -h, --help            show this help message and exit
  -i INPUT, --input INPUT
                        <Required> Input file for CLUMPS. Default expects
                        .maf, acetylomcs allowed with -a flag.
  -o OUTPUT_DIR, --output_dir OUTPUT_DIR
                        Output directory for CLUMPS input files.
  -a, --acetylomics     Flag to indicate input is already in protein
                        coordinates.
  --hgfile HGFILE       2bit human genome build file.
  --fasta FASTA         Protein primary sequence fasta file.
  --gpmaps GPMAPS       Genome Proteome Maps built from blast hits.
```

#### Clumps
Core clumps processing and algorithm.

```
clumps [-h] -d MUTS -m MAPS [-f MUT_FREQ] [--max_rand MAX_RAND]
      [--mut_types MUT_TYPES [MUT_TYPES ...]]
      [--sampler {UniformSampler,CoverageSampler,MutspecCoverageSampler}]
      [--mut_spectra MUT_SPECTRA] [-e HILL_EXP] [-p PANCAN_FACTOR]
      [-t TUMOR_TYPE] [-x XPO] [--threads THREADS] [--n_jobs N_JOBS]
      [--use_provided_values USE_PROVIDED_VALUES]
      [--coverage_track COVERAGE_TRACK] [-o OUT_DIR]
      [--pdb_dir PDB_DIR] [--hgfile HGFILE] [--fasta FASTA]
      [--gpmaps GPMAPS] [--slurm]

Run CLUMPS.

optional arguments:
  -h, --help            show this help message and exit
  -d MUTS, --muts MUTS  <Required> Directory of files titled with Uniprot IDs
                        that have mutation information
  -m MAPS, --maps MAPS  <Required> File mapping uniprot ID to PDB ID with
                        residue-level mapping information.
  -f MUT_FREQ, --mut_freq MUT_FREQ
                        Mutational frequenices of patient samples.
  --max_rand MAX_RAND   Maximum number of random samples.
  --mut_types MUT_TYPES [MUT_TYPES ...]
                        Mutation Types (N = nonsense, S = synonymous, M =
                        mutation)
  --sampler {UniformSampler,CoverageSampler,MutspecCoverageSampler}
                        Sampler to use for Null Model.
  --mut_spectra MUT_SPECTRA
                        Mutational spectra of patient samples.
  -e HILL_EXP, --hill_exp HILL_EXP
                        Hill Exponent
  -p PANCAN_FACTOR, --pancan_factor PANCAN_FACTOR
                        Pan Cancer factor, value between 0 - 1. 1.0 if tumor
                        type specified.
  -t TUMOR_TYPE, --tumor_type TUMOR_TYPE
                        Tumor type to run clumps on. PanCan indicates running
                        clumps over all tumor types found in input .maf.
  -x XPO, --xpo XPO     Soft threshold parameter for truncated Gaussian.
  --threads THREADS     Number of threads for sampling.
  --n_jobs N_JOBS       Number of jobs to use.
  --use_provided_values USE_PROVIDED_VALUES
                        Compute mutational frequencies with provided values.
  --coverage_track COVERAGE_TRACK
                        Coverage track for null sampler.
  -o OUT_DIR, --out_dir OUT_DIR
                        Output directory.
  --pdb_dir PDB_DIR     PDB directory.
  --hgfile HGFILE       2bit human genome build file.
  --fasta FASTA         Protein primary sequence fasta file.
  --gpmaps GPMAPS       Genome Proteome Maps built from blast hits.
  --slurm               Use SLURM backend.
```

### Post-process
Generates summary files from array outputs of `clumps`.

```
clumps-postprocess [-h] -i INPUT_DIR -d PROTEINS_DIR [-p UNIPROT_MAP]
                          [-c CANCER_GENES] [--pdb_dir PDB_DIR]
                          [-o OUTPUT_FILE]

Post-process CLUMPS.

optional arguments:
  -h, --help            show this help message and exit
  -i INPUT_DIR, --input_dir INPUT_DIR
                        <Required> Results directory from CLUMPS.
  -d PROTEINS_DIR, --proteins_dir PROTEINS_DIR
                        <Required> Split protein directory with mutation
                        mapping information.
  -p UNIPROT_MAP, --uniprot_map UNIPROT_MAP
                        Aggregate file of all uniprot mappings to PDB IDs with
                        residue-level information.
  -c CANCER_GENES, --cancer_genes CANCER_GENES
                        List of cancer genes, tab-delimited.
  --pdb_dir PDB_DIR     Directory of PDB files for parsing headers.
  -o OUTPUT_FILE, --output_file OUTPUT_FILE
                        Output file from CLUMPS.

```

---

### In-progress

- [x] Add support for mutspec and coverage samplers (issue with `.jar` files)
- [x] Write acetyolomics pipeilne into processing script
- [x] Re-wite post-processing
- [ ] Finish DockerFile
- [ ] Example data/test cases
- [ ] Support for slurm for `clumps.py` via `canine`
- [ ] Update `genomeProteomeMaps.txt` files with recent PDB database; add Hg38 support
- [ ] Create `build` module to download relevant reference files
- [ ] Create `pymol` module for visualization
