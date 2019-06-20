# getzlab-CLUMPS2

Re-factored code for CLUMPS written in Python3.

Code in `./clumps/`:
  * `mapper.py`: parses through `hg19.2bit`, `genomeProteomeMaps.txt`,`UP000005640_9606.fasta.gz` to build mapping of protein residues to genomic coordinates
  * `spectra.py`: functions for computing mutational context
  * `clumps.py`: CLUMPS engine
  * `postprocess.py`: results-parser
  * `utils.py`: util functions
  * `samples/`: other samplers for different null models of clumps


TODO:

- [ ] Add support for mutspec and coverage samplers (issue with `.jar` files)
- [ ] Write acetyolomics pipeilne into processing script
- [ ] Re-wite post-processing
- [ ] Example data/test cases
- [ ] Support for slurm for `clumps.py` via `canine`
- [ ] Update `genomeProteomeMaps.txt` files with recent PDB database
- [ ] Add convenient download for reference files
- [ ] Add support for Hg38
