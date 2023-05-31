import pandas as pd


maf = pd.read_csv('MC3-14k-merge_v1.maf',sep='\t')

maf['Variant_Classification'] = maf['Variant_Type']

#Start_Codon_SNP
maf.loc[maf['Variant_Classification']=='mis','Variant_Classification'] = 'Missense'
maf.loc[maf['Variant_Classification']== 'ncd','Variant_Classification'] = 'Intron'
maf.loc[maf['Variant_Classification']== 'non','Variant_Classification'] = 'Nonsense'
maf.loc[maf['Variant_Classification']== 'spl','Variant_Classification'] = 'Splice_site'
maf.loc[maf['Variant_Classification']== 'syn','Variant_Classification'] = 'Silent'
maf['Variant_Type'] = 'SNP'
maf['patient'] = maf['Tumor_Sample_Barcode']

maf.to_csv('MC3-14k-merge_v1.variant_classification.maf',sep='\t',index=None)

