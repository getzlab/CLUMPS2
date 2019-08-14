
#EMPRINT

EMPRINT is a companion to CLUMPS which instead of looking at 3d clustering in individual protein structures looks at interactions of protein-protein boundaries (PP) or interactions between proteins and other ligands, small molecules or nucleotides (PO).


example commands to run:

For each line in emprint/interface_structure_lists/


##interfaces2.py

Script to be run for each test to be done, either looking at a single protein structure in a pdb file and its interactions with other ligands, small molecules or nucleotides:


python interfaces2.py ../set/pmbl.run2 1htv B:F8WCM5:25

Or looking at the interaction of multiple protein structures:

python interfaces2.py ../set/pmbl.run2 2hck B:P09619:803_A:P46109:124

##interfaces_postprocess.py

Script to be run after each individual test is done to aggregate results

##selectStructuresForIFanalysis.py

Parses huniprot database to map pdb structures to uniprot files.  Each line of the resulting files (current versions residing in the interface_structure_lists directory) can be used as an argument to interfaces2.py

#dat

Necessary files referenced by the scripts included.  May need to change some filepaths