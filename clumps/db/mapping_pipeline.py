import os

from mapping.ref import download_mapping_ref_files
from mapping.ref import create_sifts_db
from mapping.ref import split_human_fasta
from mapping.ref import blast_uniprot_seq

from mapping.alignment import Blaster
from mapping.alignment import PDBStore
from mapping.alignment import SIFTStore

REFERENCE_DIR = "ref"
PDB_DIR = os.path.join(REFERENCE_DIR, "pdbs")
SIFTS_DIR = os.path.join(REFERENCE_DIR, "sifts")

# -------------------------------
# Download Reference Files
#
# Specs: 1-core; subject to download speed
# TODO: parallel download each file
# -------------------------------
#download_mapping_ref_files(REFERENCE_DIR)

# -------------------------------
# Split Human Fastas into Uniprots
#
# Specs: 1-core; fast
# -------------------------------
#split_human_fasta(os.path.join(REFERENCE_DIR, "UP000005640_9606.fasta.gz"))

# -------------------------------
# Create custom Seq. Data Base
#
# Specs: 1-core; very slow
# -------------------------------
# create_sifts_db(
#     os.path.join(REFERENCE_DIR, "pdb_chain_uniprot.tsv.gz"),
#     [
#         os.path.join(REFERENCE_DIR, "uniprot_sprot.fasta.gz"),
#         os.path.join(REFERENCE_DIR, "uniprot_trembl.fasta.gz")
#     ]
# )

# -------------------------------
# Blast all Sequences
#
# Specs: 32-core high cpu (64 threads) ~ 15 mins
# -------------------------------
# blast_uniprot_seq(os.path.join(REFERENCE_DIR, "uniprot.sifts.custom.db.seq"), n_threads=64)

# -------------------------------
# Load Blast Outputs & Alignments
#
# Specs: 1-core ~12 mins
# -------------------------------
# Set Up Blaster Object
blaster = Blaster(
    os.path.join(REFERENCE_DIR, "uniprot.human.sp.blasted"),
    os.path.join(REFERENCE_DIR, "pdb_chain_uniprot.tsv.gz")
)
print("Number of uniprots: {}".format(blaster.no_uniprots))
print("Number of PDBs: {}".format(len(blaster.pdbs)))

# -------------------------------
# Download PDBs
#
# Specs: 32-core high cpu (128 threads) ~ 15 mins
# -------------------------------
pdbstore = PDBStore(PDB_DIR)

# Get PDBs needed to download
pdbs_to_dl = blaster.pdbs - pdbstore.downloaded_pdbs
res = pdbstore.download_missing_pdbs(pdbs_to_dl, n_threads=128)

# One Last Attempt
res_missing = pdbstore.download_missing_pdbs(res, n_threads=128)

# Save Missing .pdbs
pd.DataFrame(list(res_missing), columns=['missing']).to_csv(os.path.join(PDB_DIR, "missing.tsv"),sep='\t')

# -------------------------------
# Download SIFTS
#
# Specs: 32-core high cpu (64 threads) ~ 15 mins
# -------------------------------
siftstore = SIFTStore(SIFTS_DIR)

# Get PDBs needed to download
pdbs_to_dl = blaster.pdbs - siftstore.downloaded_pdbs
res = siftstore.download_missing_pdbs(pdbs_to_dl, n_threads=128)

# One Last Attempt
res_missing = siftstore.download_missing_pdbs(res, n_threads=128)

# Save Missing .pdbs
pd.DataFrame(list(res_missing), columns=['missing']).to_csv(os.path.join(SIFTS_DIR, "missing.tsv"),sep='\t')

# -------------------------------
# Mapping & Cross Referencing PDBs
#
# Specs: 32-core high cpu (64 threads) ~ 15 mins
# -------------------------------
