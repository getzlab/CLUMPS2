import scipy as sp

from gzip import GzipFile
from lib import getPdbsumPPIs


## extract the native structures
fi = GzipFile('../res/huniprot2pdb.run18.filt.txt.gz')
native = set([])
while 1:
    l = fi.readline()
    if not l:
        break
    l = l.split('\t', 4)
    if l[3] == '-':
        native.add(l[2][:4])

fo = file('../dat/pdbs_forIFan_native', 'w')
fo.write('\n'.join(native) + '\n')
fo.close()

fi.close()

print 'natives done.'


## get uniprot entry_name to id
ue2id = {i[1]:i[0] for i in map(lambda x:x.strip().split(), file('../dat/uniprot_id_entryname').readlines())}

## get CPDB interactions
cpdb = set([])
fi = GzipFile('../dat/ConsensusPathDB_human_PPI.gz')
x = fi.readline()  ## header
x = fi.readline()
cnt = 0
while 1:
    l = fi.readline()
    if not l:
        break
    pp = l.strip('\n').split('\t')[2].split(',')
    ppl = len(pp)
    if ppl > 5:
        continue
    for i in xrange(ppl):
        if '.' in pp[i]:
            continue
        if pp[i] not in ue2id:
            #print pp[i], 'missing'
            continue
        a = ue2id[pp[i]]
        for j in xrange(i+1,ppl):
            if '.' in pp[j]:
                continue
            if pp[j] not in ue2id:
                #print pp[j], 'missing'
                continue
            b = ue2id[pp[j]]
            cpdb.add(tuple(sorted([a,b])))

fi.close()

## do the selection of structures
pdb2prots = {}
fi = GzipFile('../res/huniprot2pdb.run18.txt.gz')
while 1:
    l = fi.readline()
    if not l:
        break
    l = l.split('\t', 4)
    pdb,ch,res1 = l[2][:4], l[2][5], l[4].split(':',1)[0]
    iden = 100.0
    native = 1
    if l[3] != '-':
        native = 0
        for term in l[3].split():
            if term.startswith('pdb_identity:'):
                iden = float(term.rsplit(':',1)[1])
    if pdb not in pdb2prots:
        pdb2prots[pdb] = []
    pdb2prots[pdb].append((l[0],ch,iden,native,res1))

## do the selection
interlogs = {}
highiden = {}
highiden_po = set([])  ## high-identity proteins for protein-DNA/RNA/small-molecule tests
cnt = 0
for pdb in pdb2prots:
    if not cnt % 1000:
        print cnt
    cnt += 1
    ## get interface size
    pdbsum = getPdbsumPPIs(pdb)
    pp = pdb2prots[pdb]
    ppl = len(pp)
    for i in xrange(ppl):
        au,ach,aiden,anat,ares1 = pp[i]
        for j in xrange(i+1, ppl):
            bu,bch,biden,bnat,bres1 = pp[j]
            if au == bu:            ## same protein
                continue
            if ach == bch:          ## same chain mapped to two proteins
                continue
            if anat and bnat:       ## both chains are native; this case has been handled
                continue
            if min(ach,bch) not in pdbsum or max(ach,bch) not in pdbsum[min(ach,bch)]:
                continue            ## no common interface
            aubu = tuple(sorted([au,bu]))
            if aubu in cpdb:
                di = interlogs
            elif aiden >= 80 and biden >= 80:
                di = highiden
            else:
                continue
            if ach < bch:
                ## number of interface residues, weighted by corresponding pdb_identity
                ifResid = aiden*len(pdbsum[ach][bch][0]) + biden*len(pdbsum[ach][bch][1])
            else:
                ifResid = aiden*len(pdbsum[bch][ach][1]) + biden*len(pdbsum[bch][ach][0])
            if not di.has_key(aubu) or ifResid > di[aubu][3]:
                di[aubu] = (pdb, pp[i], pp[j], ifResid)
        if not anat and aiden >= 80:
            highiden_po.add((pdb,ach,au,ares1))

print 'all high-identity cases for PO analysis:', len(highiden_po)

for fn,di in [('pdbs_forIFan_interlogs', interlogs),
              ('pdbs_forIFan_highIdentity', highiden)]:
    fo = file('../dat/%s' % fn, 'w')
    for k in di:
        vv = di[k]
        l = '%s %s:%s:%s_%s:%s:%s\n' % (vv[0], vv[1][1], vv[1][0], vv[1][4], vv[2][1], vv[2][0], vv[2][4])
        fo.write(l)
        if (vv[0], vv[1][1], vv[1][0], vv[1][4]) in highiden_po:
            highiden_po.remove((vv[0], vv[1][1], vv[1][0], vv[1][4]))
        if (vv[0], vv[2][1], vv[2][0], vv[2][4]) in highiden_po:
            highiden_po.remove((vv[0], vv[2][1], vv[2][0], vv[2][4]))
    fo.close()


## add all structures that have DNA, RNA or small molecule (and are not native protein)
print 'remaining high-identity cases for PO analysis:', len(highiden_po)
fo = file('../dat/pdbs_forIFan_POhighIdentity', 'w')
for k in highiden_po:
    fo.write('%s %s:%s:%s\n' % k)

fo.close()
