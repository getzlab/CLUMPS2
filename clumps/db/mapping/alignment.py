import sys
import os
import pandas as pd
from collections import defaultdict
from tqdm import tqdm
from lxml import etree
import gzip
from typing import Union

from agutil.parallel import parallelize2
import subprocess

from ...utils import gunzipper, AMINO_ACID_MAP

"""
Contains 3 classes for protein-sequence alignment.

    Blaster: blasting protein sequence data
    PDBStore: PDB-file store
    SIFTStore: SIFTS-file store
"""

class Blaster(object):
    """
    Blasted.
    ------------------
    This class is meant for organization of aligned sequence data. It requires mapping
    information from the Protein SIFTS database and results from blastp.
    """
    def __init__(self, blast_dir, sifts_file, eval_thresh=0.001, ident_thresh=0.15):
        """
        Inputs:
            * blast_dir: a directory of blasted outputs for each uniprot in ".seq.blasted.gz" format
            * sifts_file: sifts mapping of protein IDs to residue mapping
                "wget ftp://ftp.ebi.ac.uk/pub/databases/msd/sifts/flatfiles/tsv/pdb_chain_uniprot.tsv.gz"
            * eval_thresh: structure evaluation threshold
            * ident_thresh: threshold of sequence similarity

        Attributes:
            * sifts_df: pandas dataframe of all sifts inputs from "sifts_file"
            * sifts_id_maps: dictionary mapping primary uniprot ID to a pdb-chain
            * maps: mapping of primary uniprot to secondary (aligned) uniprot and PDB-Chain
            * self.human_uniprot_ids: all uniprot IDs blasted

        Properties:
            * no_uniprots: number of uniprots found
            * pdbs: number of mapped PDBs
        """
        self.blast_dir = blast_dir
        self.sifts_file = sifts_file
        self.eval_thresh = eval_thresh
        self.ident_thresh = ident_thresh

        # Uniprots
        self.human_uniprot_ids = [
            i.split('.',1)[0] for
            i in
            os.listdir(self.blast_dir)
            if i.endswith('.blasted.seq.gz')
        ]

        # Sifts
        self.sifts_id_maps = defaultdict(set)
        self.sifts_df = pd.read_csv(self.sifts_file, skiprows=1, sep='\t')

        for idx,row in tqdm(self.sifts_df.iterrows(), desc='Readings sifts', total=self.sifts_df.shape[0]):
            self.sifts_id_maps[row['SP_PRIMARY']].add((row['PDB'], row['CHAIN']))

        # Get direct & indirect uniprot maps
        self.maps = self.get_u1u2_maps()

        # Create dataframe
        _map = list()
        for u1 in self.maps:
            for u2 in self.maps[u1]:
                for pdbch in self.maps[u1][u2]:
                    _map.append((u1,u2,pdbch))

        self.maps_df = pd.DataFrame(_map).rename(columns={0:'u1',1:'u2',2:'pdbch'})

    @property
    def no_uniprots(self):
        return len(self.human_uniprot_ids)

    @property
    def pdbs(self):
        """
        Get mapped PDBs.
        ------------------
        Iterates through maps from mutated proteins to return a set of
        covered pdbs.
        """
        pdbs = set()
        for u1 in self.maps:
            for u2 in self.maps[u1]:
                for k in self.maps[u1][u2]:
                    pdbs.add(k[0])
        return pdbs

    def get_u1u2_maps(self):
        """
        Get U1-U2-PDB Mappings.
        ------------------
        Iterating through alignments can take ~10 minutes.
        TODO: optimize by multithreading.
        """
        u1u2pdbch = defaultdict(lambda: defaultdict(dict))

        # Find direct maps
        for u1 in self.human_uniprot_ids:
            if u1 in self.sifts_id_maps:
                for k in self.sifts_id_maps[u1]:
                    u1u2pdbch[u1][u1][k] = []

        # Find indirect maps
        self.alignments = defaultdict(lambda: defaultdict(list))

        for u1 in tqdm(self.human_uniprot_ids, desc='Mapping indirect maps'):
            alis = self.get_alignments(u1)

            for ali in alis:
                # Aligned Uniprot ID
                u2 = ali['Hit_def'].split('|',2)[1]
                self.alignments[u1][u2].append(ali)

                if u2 in self.sifts_id_maps:
                    for k in self.sifts_id_maps[u2]:
                        u1u2pdbch[u1][u2][k] = []

        return u1u2pdbch

    def get_alignments(self, uniprot_id):
        """
        Get alignments.
        ------------------
        Parses through protein blast output; filters above a given hsp_evalue
        threshold and a provided relative identity threshold.

        Args:
            * uniprot_id: Uniprot ID of a given protein

        Returns:
            * alignments
        """
        results = list()

        try:
            with gzip.open(os.path.join(self.blast_dir, uniprot_id + '.blasted.seq.gz'), 'r') as f:
                cur = etree.fromstring(f.read())
                qlen = 0  ## query length
                for i in cur:
                    if i.tag == 'BlastOutput_iterations':
                        cur = i
                    elif i.tag == 'BlastOutput_query-len':
                        qlen = int(i.text)

                cur = cur[-1]  ## last iteration

                for i in cur:
                    if i.tag == 'Iteration_hits':
                        cur = i
                        break

                for hit in cur:
                    hit_def = None
                    hsps = None

                    for i in hit:
                        if i.tag == 'Hit_def':
                            hit_def = i.text
                        elif i.tag == 'Hit_hsps':
                            hsps = i

                    for hsp in hsps:
                        di = {}
                        di['Hit_def'] = hit_def

                        for i in hsp:
                            di[i.tag] = i.text

                        if float(di['Hsp_evalue']) > self.eval_thresh:  ### !!! FILTER !!!
                            continue

                        rel_ident = float(di['Hsp_identity'])/float(di['Hsp_align-len'])

                        if rel_ident < self.ident_thresh:         ### !!! FILTER !!!
                            continue

                        di['relIdentity'] = rel_ident
                        results.append(di)
        except:
            pass

        return results

class PDBStore(object):
    """
    Protein Data Bank Direcotry Store
    ------------------
    This object is used to keep track of protein data bank files
    that are downloaded in a reference directory.
    """
    def __init__(self, path):
        """
        Args:
            * path: PDB directory

        """
        from prody import parsePDBStream

        self.pdb_root_dir = path
        self.pdb_dir = os.path.join(path, 'ftp.wwpdb.org/pub/pdb/data/structures/divided/pdb/')

    def __str__(self):
        return "PDBStore\n   * {} PDB files downloaded\n   * Directory Path: {}".format(len(self.downloaded_pdbs), self.pdb_root_dir)

    @property
    def downloaded_pdbs(self):
        """
        Get Downloaded Structures.
        ------------------
        Iterates through pdb structure direcotry to find what
        pdbs are already downloaded. Returns a set of downloaded PDBs.
        """
        try:
            return {
                    i[3:7]
                    for j in os.listdir(self.pdb_dir)
                    for i in os.listdir(os.path.join(self.pdb_dir, j))
                    if i.endswith('.ent.gz') and os.path.getsize(os.path.join(self.pdb_dir, j, i)) > 0
                   }
        except:
            return set()

    def missing_pdbs(self, pdbs: set):
        """
        Missing PDBs
        ------------------
        Set difference of input pdbs and all downloaded pdbs.

        Args:
            * pdbs: set of PDBs
        """
        return pdbs - self.downloaded_pdbs

    def download_missing_pdbs(self, pdbs_to_download: Union[set,list], n_threads=15):
        """
        Download Missing PDBs
        ------------------
        Args:
            * pdbs_to_download: a list or set of PDBs to download
            * n_threads: number of threads to use for downloads
        """
        if isinstance(pdbs_to_download, list):
            pdbs_to_download = set(pdbs_to_download)

        out_dir = self.pdb_root_dir

        @parallelize2(maximum=n_threads)
        def dl_pdb(pdb):
            """Download pdb."""
            try:
                cmd = "wget -qr ftp://ftp.wwpdb.org/pub/pdb/data/structures/divided/pdb/{}/pdb{}.ent.gz -P {}".format(pdb[1:3], pdb, out_dir)
                output = subprocess.check_call(cmd, executable='/bin/bash', shell=True)
                return None
            except:
                return pdb

        print("   * Downloading {} pdbs using {} threads".format(len(pdbs_to_download), n_threads))

        tmp = [dl_pdb(pdb) for pdb in pdbs_to_download]
        pdb_err = {callback() for callback in tmp}

        print("   * Downloaded {}/{} successfully.".format(len(pdbs_to_download-pdb_err) , len(pdbs_to_download)))

        return pdb_err

    def load(self, pdb, chain=None):
        """
        Load residues - amino acids.

        Returns:
            * dictionary of residue position to amino-acid
        """
        pdb_file = os.path.join(self.pdb_dir, pdb[1:3], "pdb{}.ent.gz".format(pdb))

        with gunzipper(pdb_file) as pfile:
            aa = parsePDBStream(pfile, chain=chain)

        res_map = dict(zip(aa.getResnums(),aa.getResnames()))
        return {k:v for k,v in res_map.items() if v in AMINO_ACID_MAP}

class SIFTStore(object):
    """
    Sift Store
    ------------------
    This object is used to keep track of sifts protein data bank files
    that are downloaded in a reference directory.
    """
    def __init__(self, sifts_dir):
        """
        Args:
            * sifts_dir
        """
        self.sifts_root_dir = sifts_dir
        self.sifts_dir = os.path.join(sifts_dir, 'ftp.ebi.ac.uk/pub/databases/msd/sifts/split_xml/')

    def __str__(self):
        return "SIFTStore\n   * {} PDB files downloaded\n   * Directory Path: {}".format(len(self.downloaded_pdbs), self.sifts_root_dir)

    @property
    def downloaded_pdbs(self):
        """
        Get downloaded SIFTs maps.
        ------------------
        Iterates through maps from SIFTS database to provide set of PDBs.
        """
        try:
            return {
                    i[:4]
                    for j in os.listdir(self.sifts_dir)
                    for i in os.listdir(os.path.join(self.sifts_dir, j))
                    if i.endswith('.xml.gz')
                   }
        except:
           return set()

    def download_missing_sifts(self, pdbs_to_download: Union[set,list], n_threads=15):
        """
        Download Missing PDBs
        ------------------
        Args:
            * pdbs_to_download: a list or set of PDBs to download
            * n_threads: number of threads to use for downloads
        """
        if isinstance(pdbs_to_download, list):
            pdbs_to_download = set(pdbs_to_download)

        out_dir = self.sifts_root_dir

        @parallelize2(maximum=n_threads)
        def dl_pdb(pdb):
            """Download pdb."""
            try:
                cmd = "wget -qr ftp://ftp.ebi.ac.uk/pub/databases/msd/sifts/split_xml/{}/{}.xml.gz -P {}".format(pdb[1:3], pdb, out_dir)
                output = subprocess.check_call(cmd, executable='/bin/bash', shell=True)
                return None
            except:
                return pdb

        print("   * Downloading {} pdbs using {} threads".format(len(pdbs_to_download), n_threads))

        tmp = [dl_pdb(pdb) for pdb in pdbs_to_download]
        pdb_err = {callback() for callback in tmp}

        print("   * Downloaded {}/{} successfully.".format(len(pdbs_to_download-pdb_err) , len(pdbs_to_download)))

        return pdb_err

    def get_map(self, uniprot, pdb, chain, verbose=False):
        """
        Get SIFTS Map.
        ------------------
        Parses through a SIFTS file.
        Gets the Sifts uniprot to PDB residue mapping.
        """
        try:
            with gzip.open(os.path.join(self.sifts_dir,  pdb[1:3], pdb.lower() + '.xml.gz'), 'r') as f:
                if verbose: print("Parsing {} | {}-{}".format(uniprot,pdb,chain))

                ret = {}
                tree = etree.fromstring(f.read())

                for i in range(len(tree)):
                    if tree[i].tag.split('}')[-1] != 'entity':
                        continue

                    for j in range(len(tree[i])):
                        if tree[i][j].tag.split('}')[-1] != 'segment':
                            continue

                        for k in range(len(tree[i][j])):
                            if tree[i][j][k].tag.split('}')[-1] != 'listResidue':
                                continue

                            for resid in tree[i][j][k].iterchildren():
                                pdb_id = None
                                pdb_rn = None  ## residue number
                                pdb_aa = None
                                uniprot_id = None
                                uniprot_rn = None  ## residue number
                                uniprot_aa = None

                                for ref in resid.iterchildren():
                                    if ref.tag.split('}')[-1] != "crossRefDb":
                                        continue

                                    if ref.attrib['dbSource'] == 'PDB':
                                        if chain != ref.attrib['dbChainId']:
                                            continue

                                        pdb_id = ref.attrib['dbAccessionId']
                                        pdb_chain = ref.attrib['dbChainId']
                                        pdb_aa = ref.attrib['dbResName']

                                        try:
                                            pdb_rn = int(ref.attrib['dbResNum'])
                                        except:
                                            pdb_rn = None

                                    elif ref.attrib['dbSource'] == 'UniProt':
                                        if uniprot != ref.attrib['dbAccessionId']:
                                            continue
                                        uniprot_id = ref.attrib['dbAccessionId']
                                        uniprot_aa = ref.attrib['dbResName']
                                        uniprot_rn = int(ref.attrib['dbResNum'])


                                if pdb_rn is not None and uniprot_rn is not None:
                                    if uniprot_rn in ret and ret[uniprot_rn] != pdb_rn:
                                        ## should not happen; otherwise: redundant maps
                                        print(pdbid, chain, uniprot, ret[uniprot_rn], (pdb_rn, pdb_aa))
                                    #ret[uniprot_rn] = (pdb_rn, pdb_aa)
                                    ret[uniprot_rn] = pdb_rn
            return ret
        except:
            if verbose: print("Not downloaded: {} | {}-{}".format(uniprot,pdb,chain))
            return None
