import os
from sys import stdout
import subprocess
from tqdm import tqdm
import gzip
from typing import Union
import glob
import pandas as pd

from agutil.parallel import parallelize2

def download_mapping_ref_files(
    REFERENCE_DIR: str,
    verbose: bool = True
    ):
    """
    Download reference files
    --------------------------
    Downloads references files to provided directory.
    """
    os.makedirs(REFERENCE_DIR, exist_ok=True)

    files = (
        "ftp://ftp.ebi.ac.uk/pub/databases/msd/sifts/flatfiles/tsv/pdb_chain_uniprot.tsv.gz",
        "ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.fasta.gz",
        "ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_trembl.fasta.gz",
        "ftp://ftp.ebi.ac.uk/pub/databases/reference_proteomes/QfO/Eukaryota/UP000005640_9606.fasta.gz"
    )

    for file in files:
        if verbose:
            stdout.write("   * Downloading {}...".format(file.split('/')[-1]))
        cmd = "wget -P {} {}".format(REFERENCE_DIR, file)
        subprocess.call(cmd, executable='/bin/bash', shell=True)

def f_length(file_name: str):
    """
    Check file length
    --------------------------
    Used for iterating throughlarge files for easy tracking.
    """
    return int(subprocess.run("wc -l {}".format(file_name), shell=True, executable="/bin/bash", stdout=subprocess.PIPE).stdout.decode().strip().split(' ')[0])

def create_sifts_db(
    sifts_file: str,
    fastqs: list,
    custom_db_output_file: Union[None,str] = None,
    verbose: bool = True
    ):
    """
    Create Sifts Database.
    -------------------------
    See download_mapping_ref_files

    Args:
        * sifts_file: pdb_chain_uniprot.tsv.gz - used to select uniprot_ids
        * fastqs: uniprot fastqs to process through
        * custom_db_output_file: custom database output_file
        * verbose: verbosity

    Returns:
        None
    """
    try:
        uniprot_ids = set(pd.read_csv(sifts_file, sep='\t', skiprows=1)['SP_PRIMARY'])
    except:
        FileExistsError("Provide valid sifts-file.")

    if verbose:
        print("   * Found {} Uniprots".format(len(uniprot_ids)))

    if custom_db_output_file is None:
        custom_db_output_file = os.path.join(os.path.dirname(sifts_file), "uniprot.sifts.custom.db.seq")

    with open(custom_db_output_file,'w') as f_out:
        for in_file in fastqs:
            with gzip.open(in_file, 'r') as f:
                annot = f.readline().decode('utf-8')
                seq = list()

                _f_len = f_length(in_file)
                for line in tqdm(f, total=_f_len, desc=in_file):
                    line = line.decode('utf-8')

                    if line.startswith('>'):
                        if annot.split('|')[1] in uniprot_ids:
                            f_out.write(annot + ''.join(seq))

                        annot = line
                        seq = list()

                    else:
                        seq.append(line)

def split_human_fasta(
    up_proteome_file: str,
    outfile: Union[None,str] = None
    ):
    """
    Split Human Fasta
    ------------------------
    Download the following human fasta:
    # ! wget ftp://ftp.ebi.ac.uk/pub/databases/reference_proteomes/QfO/Eukaryota/UP000005640_9606.fasta.gz

    Splits the human proteome into a fastsa of sequences mapped with their
    uniprot IDs. Creates the an output file with all proteins <outfile> and a directory
    with individual protein sequences.

    Args:
        * up_proteome_file: path to proteome file
        * outfile: path to output_file

    Returns:
        * None
    """
    if outfile is None:
        outfile = os.path.join(os.path.dirname(up_proteome_file), 'uniprot.human.sp.fa')

    seq_output_dir = os.path.join(os.path.dirname(up_proteome_file), "uniprot.human.sp")
    os.makedirs(seq_output_dir, exist_ok=True)

    with open(outfile, 'w') as f_out:
        with gzip.open(up_proteome_file,'r') as f_in:
            annot = f_in.readline().decode('utf-8').strip()
            seq = list()

            _f_len = f_length(outfile)
            for line in tqdm(f_in, total=_f_len):
                line = line.decode('utf-8').strip()

                if line.startswith('>'):
                    if annot.startswith('>sp|'):
                        annot = annot.split('|')
                        if '-' not in annot[1]:
                            f_out.write('>%s\n' % annot[1])
                            f_out.write(''.join(seq) + '\n')

                            with open(os.path.join(seq_output_dir, '{}.seq'.format(annot[1])), 'w') as f_ind:
                                f_ind.write(''.join(seq) + '\n')
                    annot = line
                    seq = list()
                else:
                    seq.append(line)

def blast_uniprot_seq(
    seq: str,
    uniprots_dir: Union[None, str] = None,
    uniprots_blasted_dir: Union[None, str] = None,
    db_title: str = "clumps_custom_db",
    verbose: bool = True,
    n_threads=15
    ):
    """
    Blast Uniprot Sequences
    --------------------------
    First, install blast:
        sudo apt-get install ncbi-blast+

    Args:
        * seq: path to sequence database
        * uniprots_dir: directory of uniprot seq
        * uniprots_blasted_dir:
        * db_title: title of database

    Returns:
        * None
    """
    if uniprots_dir is None:
        uniprots_dir = os.path.join(os.path.dirname(seq), "uniprot.human.sp")

    if uniprots_blasted_dir is None:
        uniprots_blasted_dir = os.path.join(os.path.dirname(seq), "uniprot.human.sp.blasted")

    # Output Blasted Dir
    os.makedirs(uniprots_blasted_dir, exist_ok=True)

    # Custom Blast Database Directory
    sifts_db_dir = os.path.join(os.path.dirname(seq), "siftsdb")
    os.makedirs(sifts_db_dir, exist_ok=True)

    # Make Blast Database
    cmd = "makeblastdb -in {} -dbtype 'prot' -out {}/siftsdb -title {}".format(seq, sifts_db_dir, db_title)
    db_out = subprocess.run(cmd, executable='/bin/bash', shell=True, stdout=subprocess.PIPE).stdout.decode()

    if verbose: print(db_out)

    @parallelize2(maximum=n_threads)
    def blast(seq):
        """Blast Command."""
        uid = seq.split('/')[-1].split('.seq')[0]

        cmd = "blastp -query {} -db {} -out {} -outfmt 5 && gzip -f {} ;".format(
            seq,
            os.path.join(sifts_db_dir, "siftsdb"),
            os.path.join(uniprots_blasted_dir, "{}.blasted.seq".format(uid)),
            os.path.join(uniprots_blasted_dir, "{}.blasted.seq".format(uid))
        )
        subprocess.call(cmd, executable='/bin/bash', shell=True)

    print("   * Running blast sequences using {} threads".format(n_threads))
    seq_files = glob.glob("{}/*.seq".format(uniprots_dir))

    tmp = [blast(s) for s in seq_files]
    results = [callback() for callback in tmp]
