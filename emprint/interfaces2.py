import os
import sys
import urllib2

from copy import deepcopy
from gzip import GzipFile
from lib import *
from samplers.MutspecCoverageSampler import *
#from samplers.UniformSampler import *


PARAM = sys.argv[1]
PDB = sys.argv[2]
if len(sys.argv) > 3:
    fix_chain2up = {i[0]:(i[1],i[2]) for i in map(lambda x:x.split(':'), sys.argv[3].split('_'))}  ## A:P12345:1_B:P67890:15  or  A:P12345:1 (in second case, only protein-DNA/RNA/small-molecule will be tested)
else:
    fix_chain2up = None

PROXIM = 5  ## max centroid proximity for interface band control

fi = file(PARAM)
for l in fi.readlines():
    l = l.strip('\n')
    exec(l)

if len(PDB) > 4:
    PDB,CH = PDB.split('-')
    CH = CH.split('_')
else:
    CH = None


if TTYPE == 'all':
    TTYPE = None

chain2up = {}
chain2mappedresid = {}  ## mapped residues in PDB coordinates
pdb2upresmaps = {}

## get all chain to uniprot maps
if not fix_chain2up:  ## all native maps are available in the filtered set and hence through the ws
    maps = [l.rsplit('\t',1)[0] for l in urllib2.urlopen(MAPURL + PDB).readlines()]  ## trim off the "filt" flag
else:  ## could be a non-native chain.  need to parse the maps file.
    maps = [l for l in GzipFile('../res/huniprot2pdb.run18.split/%s.gz' % PDB[1:3]).readlines() if l.split('\t',2)[2].startswith(PDB)]
for l in maps:
    u1,u2,pdbch,alidt,resmap = l.rstrip('\n').split('\t',5)
    pdb,ch = pdbch.split('-')
    if fix_chain2up and (ch not in fix_chain2up or (u1, resmap.split(':',1)[0]) != fix_chain2up[ch]):
        continue
    elif not fix_chain2up and alidt != '-' and u2 != '-':  ## selects only direct maps
        continue
    d = {}
    a = [] ## mapped residues in PDB coord
    for i in resmap.split():
        i = map(int,i.split(':'))
        d[i[0]] = i[1]
        a.append(i[1])
    resmap = d
    if chain2up.has_key(ch):
        print 'fusion protein. choosing the larger sequence.'
        if len(resmap) < len(chain2up[ch][1]):
            continue
    chain2up[ch] = (u1,resmap)
    chain2mappedresid[ch] = a
    pdb2upresmaps[ch] = {resmap[i]:i for i in resmap}

## get residue interaction data
dnarna = {}
for l in file('../dat/pids_with_dnarna.txt').readlines():
    l = l.strip('\n').split('\t')
    if l[0] not in dnarna:
        dnarna[l[0]] = {}
    if l[3]:
        paired = dict([tuple(i.split(':')) for i in l[3].split(',')])
    else:
        paired = {}
    dnarna[l[0]][l[1]] = (l[2], paired)

if PDB in dnarna:
    d = dnarna[PDB]
else:
    d = {}

if 1:  ## get residue interaction data from PDBsum
    ppi = getPdbsumPPIs(PDB, chain2up)
    poi = getPdbsumPOIs(PDB, chain2up, d)

ppi_direct = deepcopy(ppi)
poi_direct = deepcopy(poi)

#else:  ## get residue interaction data from Foldx. PPI-small molecule interactions are not available from Foldx.  Code works, but pdbSUM seems to be a bit better
#    ppi, poi = getFoldxPIs(PDB, chain2up, d)

##################################################
###  BAND CONTROL START
## handle the 'band' of residues in immediate neighborhood of interface residues

if 1:  ## calculate 'band' of residues near the interface residues
    ch2dist = {}
    for ch in pdb2upresmaps:
        D, pdb_resid = getResidueDistanceMatrix((PDB,ch), point='centroid')
        pdb_resid = dict(zip(pdb_resid, range(len(pdb_resid))))
        ch2dist[ch] = (D, pdb_resid)
    #DEN = 2.0 * 3**2
    ppi_band = {}
    poi_band = {}
    for pxi,pxi_band in [(ppi, ppi_band), (poi, poi_band)]:
        for c1 in pxi:
            pxi_band[c1] = {}
            for c2 in pxi[c1]:
                pxi_band[c1][c2] = [{},{}]
                for c,idx in [(c1,0), (c2,1)]:
                    if c not in ch2dist:
                        continue  ## c2 of poi is not in ch2dist
                    D, pdb_resid = ch2dist[c]
                    for i in pxi[c1][c2][idx]:
                        for j in pdb_resid:
                            if j in pxi[c1][c2][idx]:
                                continue
                            if i not in pdb_resid or j not in pdb_resid:
                                continue
                            d = D[max(pdb_resid[i],pdb_resid[j])][min(pdb_resid[i],pdb_resid[j])]
                            if d <= PROXIM:
                                #d = sp.exp(-(d**2)/DEN)  ## if we want to have a continuous factor for residues in the band
                                #if j not in pxi_band[c1][c2][idx] or pxi_band[c1][c2][idx][j] < d:
                                pxi_band[c1][c2][idx][j] = 1 # d
        if 1:   ## add band to interface
            for c1 in pxi_band:
                for c2 in pxi_band[c1]:
                    pxi[c1][c2][0].update(pxi_band[c1][c2][0])
                    pxi[c1][c2][1].update(pxi_band[c1][c2][1])    
        elif 0: ## remove band from consideration. This requires major rewriting of the script since permutations should be done interface-specifically and not globally
            pass

###  BAND CONTROL END
##################################################

## get mutations
## CAUTION: keys in md are PDB residue numbers, not uniprot residue numbers
if SAMPLEMUTFREQWEIGHT:
    mfreq = {}
    fi = file(SAMPLEMUTFREQWEIGHT)
    l = fi.readline() ## hdr
    while 1:
        l = fi.readline()
        if not l:
            break
        l = l.strip().split('\t')
        mfreq[(l[0],l[1])] = float(l[4])

chain2muts = {}
for ch in chain2up:
    md = {}
    u1,resmap = chain2up[ch]
    gm = urllib2.urlopen(MUTURL + u1).read()
    for l in map(lambda x:x.split('\t'), gm.split('\n')):
        if l == [''] or not l[5]:
            continue
        if l[6][0] not in MUTTYPES:
            continue
        ## ttype selection
        if TTYPE and not l[0].startswith(TTYPE):
            continue
        p = int(l[5])
        if p not in resmap:
            continue
        pp = resmap[p] 
        if pp not in md: ## CAUTION: keys in md are PDB residue numbers, not uniprot residue numbers
            md[pp] = [0, None, set([])]  ## the 1-st element is there just for compatibility
        l[0] = l[0].split('-')[0]
        if SAMPLEMUTFREQWEIGHT:
            md[pp][0] += mfreq[(l[0], l[1])]
        else:
            md[pp][0] += 1
        md[pp][2].add((l[1],l[0]))
    chain2muts[ch] = md

## get identities of interface residues and mutated interface residues
ifres = {}
ifresdir = {}
ifresmut = {}
for c1 in ppi:
    for c2 in ppi[c1]:
        k = (c1,c2)
        ifres[k] = [sorted([pdb2upresmaps[c1][r1] for r1 in ppi[c1][c2][0] if r1 in pdb2upresmaps[c1]]),
                    sorted([pdb2upresmaps[c2][r2] for r2 in ppi[c1][c2][1] if r2 in pdb2upresmaps[c2]])]
        ifresdir[k] = [sorted([pdb2upresmaps[c1][r1] for r1 in ppi_direct[c1][c2][0] if r1 in pdb2upresmaps[c1]]),
                       sorted([pdb2upresmaps[c2][r2] for r2 in ppi_direct[c1][c2][1] if r2 in pdb2upresmaps[c2]])]
        ifresmut[k] = [[],[]]
        for r1 in ppi[c1][c2][0]:
            if r1 in chain2muts[c1]:
                ifresmut[k][0].append((pdb2upresmaps[c1][r1],chain2muts[c1][r1][0]))
        for r2 in ppi[c1][c2][1]:
            if r2 in chain2muts[c2]:
                ifresmut[k][1].append((pdb2upresmaps[c2][r2],chain2muts[c2][r2][0]))

for c1 in poi:
    for c2 in poi[c1]:
        k = (c1,c2)
        ifresmut[k] = [[],[]]
        ifres[k] = [sorted([pdb2upresmaps[c1][r1] for r1 in poi[c1][c2][0] if r1 in pdb2upresmaps[c1]]),
                    []]
        ifresdir[k] = [sorted([pdb2upresmaps[c1][r1] for r1 in poi_direct[c1][c2][0] if r1 in pdb2upresmaps[c1]]),
                    []]
        for r1 in poi[c1][c2][0]:
            if r1 in chain2muts[c1]:
                ifresmut[k][0].append((pdb2upresmaps[c1][r1],chain2muts[c1][r1][0]))

for k in ifres:
    ifres[k][0] = ';'.join(map(str,ifres[k][0]))
    ifres[k][1] = ';'.join(map(str,ifres[k][1]))

for k in ifresdir:
    ifresdir[k][0] = ';'.join(map(str,ifresdir[k][0]))
    ifresdir[k][1] = ';'.join(map(str,ifresdir[k][1]))

for k in ifresmut:
    ifresmut[k][0] = ';'.join(['%d_%.2f' % (i[0],i[1]) for i in sorted(ifresmut[k][0])])
    ifresmut[k][1] = ';'.join(['%d_%.2f' % (i[0],i[1]) for i in sorted(ifresmut[k][1])])

## BURIAL CONTROL
##############################
"""
chain2burialgroups = {k:([],[],[]) for k in chain2up}
chain2residue2burialgroup = {}
ch2bg = {}
chain2muts_sample = {}
for ch in CH:
    chain2residue2burialgroup[ch] = {}
    chain2muts_sample[ch] = [[],[],[]]
    for l in file('../dat/sidechain_burial/%s-%s' % (PDB,ch)).readlines()[1:]:
        l = l.strip().split('\t')
        if len(l) < 4:
            continue
        i = int(l[2])
        b = float(l[3])
        if b < 0.33:
            g = 0
        elif b < 0.66:
            g = 1
        else:
            g = 2
        chain2burialgroups[ch][g].append(i)
        chain2residue2burialgroup[ch][i] = g
        if i in chain2muts[ch]:
            chain2muts_sample[ch][g].append(chain2muts[ch][i])
    ch2bg[ch] = [[chain2residue2burialgroup[ch][k] for k in chain2muts[ch]].count(0),
                 [chain2residue2burialgroup[ch][k] for k in chain2muts[ch]].count(1),
                 [chain2residue2burialgroup[ch][k] for k in chain2muts[ch]].count(2)]
    chain2muts_sample[ch] = chain2muts_sample[ch][0] + chain2muts_sample[ch][1] + chain2muts_sample[ch][2]

for ch in ppi.keys():
    if ch != CH[0]:
        del ppi[ch]

for ch in ppi[CH[0]].keys():
    if ch != CH[1]:
        del ppi[CH[0]][ch]
"""
##############################


#### DO THE TESTS
def interfaceTest():
    ret = [{},{}]  ## ppis, pois
    ## ppi
    for c1 in ppi:
        for c2 in ppi[c1]:
            nmut_c1 = 0
            for r in ppi[c1][c2][0]:
                if r in chain2muts[c1]:
                    nmut_c1 += chain2muts[c1][r][0]
                    #nmut_c1 += chain2muts[c1][r][0] * ppi[c1][c2][0][r]
            nmut_c2 = 0
            for r in ppi[c1][c2][1]:
                if r in chain2muts[c2]:
                    nmut_c2 += chain2muts[c2][r][0]
                    #nmut_c2 += chain2muts[c2][r][0] * ppi[c1][c2][1][r]
            nmut_total = nmut_c1 + nmut_c2
            ret[0][(c1,c2)] = (nmut_c1, nmut_c2, nmut_total)
    ## poi
    for c1 in poi:
        for c2 in poi[c1]:
            nmut_c1 = 0
            for r in poi[c1][c2][0]:
                if r in chain2muts[c1]:
                    nmut_c1 += chain2muts[c1][r][0]
                    #nmut_c1 += chain2muts[c1][r][0] * poi[c1][c2][0][r]
            ret[1][(c1,c2)] = (nmut_c1,)
    return ret

#raise 7

nmut_real = interfaceTest()

def bster():
    for i in [0,1]:
        for k in P[i]:
            for v in P[i][k]:
                if booster(v, rnd):
                    return True
    return False
    
P = [{k:[0]*len(nmut_real[0][k]) for k in nmut_real[0]}, {k:[0]*len(nmut_real[1][k]) for k in nmut_real[1]}]
samplers = {}
for ch in chain2up:
    u1,resmap = chain2up[ch]
    resmap = sorted(resmap.items())
    availUPresid = [i[0] for i in resmap]
    availPDBresid = [i[1] for i in resmap]
    md = {pdb2upresmaps[ch][pp]:chain2muts[ch][pp] for pp in chain2muts[ch]}
    mireal = []
    c = 0
    for i in availUPresid:
        if i in md:
            mireal.append(c)
        c += 1
    if SAMPLER == 'UniformSampler':
        sam = UniformSampler(availUPresid)
    elif SAMPLER == 'CoverageSampler':
        sam = CoverageSampler(availUPresid, u1, COVERAGETRACK)
    elif SAMPLER == 'MutspecCoverageSampler':
        sam = MutspecCoverageSampler(availUPresid, u1, COVERAGETRACK, PATMUTSPECTRA)
        sam.calcMutSpecProbs(md)
    samplers[ch] = (sam, availUPresid, availPDBresid, mireal, [md[availUPresid[i]] for i in mireal])

rnd = 0
while rnd < MAXRAND and (rnd%10000 or bster()):
    if not rnd % 10000:
        print rnd
    ok = 1
    for ch in chain2muts:
        s = samplers[ch][0].sample(samplers[ch][3])
        if s is None:
            ok = 0
            break
        mi,idx = s
        chain2muts[ch] = dict(zip([samplers[ch][2][i] for i in mi], [samplers[ch][4][i] for i in idx]))
    if not ok:
        continue
    """
#################################
    for ch in CH:
        ll = []
        for g in [0,1,2]:
            ll += random.sample(chain2burialgroups[ch][g],ch2bg[ch][g])
        chain2muts[ch] = dict(zip(ll, chain2muts_sample[ch]))
#################################
    """
    rr = interfaceTest()
    for i in [0,1]:
        for k in rr[i]:
            for j in xrange(len(nmut_real[i][k])):
                if rr[i][k][j] >= nmut_real[i][k][j]:
                    P[i][k][j] += 1
    rnd += 1

for i in [0,1]:
    for k in P[i]:
        for j in xrange(len(P[i][k])):
            P[i][k][j] = max(1, P[i][k][j])/float(rnd)

## write output
if not os.path.exists('../res/interfaces-%s.%s/' % ((TTYPE or 'PanCan'), RESFLAG)):
    os.makedirs('../res/interfaces-%s.%s/' % ((TTYPE or 'PanCan'), RESFLAG))

if fix_chain2up:
    ad = '-'.join(['_'.join((it[0], it[1][0], it[1][1])) for it in sorted(fix_chain2up.items())])
    ofn = '%s-%s' % (PDB, ad)
else:
    ofn = PDB

fo = file('../res/interfaces-%s.%s/%s' % ((TTYPE or 'PanCan'), RESFLAG, ofn), 'w')
for k in P[0]:
    c1,c2 = k
    ln = ['P', PDB, c1, c2, chain2up[c1][0], chain2up[c2][0], 
          '', '', '', '',                       ## genes and genelists
          '',                                   ## reference seq1 length
          '%d' % len(chain2mappedresid[c1]),    ## mapped residues of chain 1
          '%d' % len([i for i in ppi[c1][c2][0] if i in pdb2upresmaps[c1]]), ## mapped interface residues of chain 1
          '',                                   ## reference seq2 length
          '%d' % len(chain2mappedresid[c2]),    ## mapped residues of chain 2
          '%d' % len([i for i in ppi[c1][c2][1] if i in pdb2upresmaps[c2]]), ## mapped interface residues of chain 2
          '%.2f' % nmut_real[0][(c1,c2)][2],    ## mutations in interface
          '%g' % P[0][(c1,c2)][2],              ## P_mutations in interface
          '',                                   ## q-value(full): stays empty
          '',                                   ## q-value(constr): stays empty
          '%.2f' % nmut_real[0][(c1,c2)][0],    ## mutations in partner 1
          '%g' % P[0][(c1,c2)][0],              ## P_mutations in partner 1
          '%.2f' % nmut_real[0][(c1,c2)][1],    ## mutations in partner 2
          '%g' % P[0][(c1,c2)][1],              ## P_mutations in partner 2
          ifres[k][0],                          ## identities of residues in interface of chain 1 (in UniProt coordinates)
          ifres[k][1],                          ## identities of residues in interface of chain 2 (in UniProt coordinates)
          ifresdir[k][0],                       ## identities of residues in direct interface of chain 1 (in UniProt coordinates)
          ifresdir[k][1],                       ## identities of residues in direct interface of chain 2 (in UniProt coordinates)
          ifresmut[k][0],                       ## identities of mutations in interface of chain 1 (in UniProt coordinates)
          ifresmut[k][1]]                       ## identities of mutations in interface of chain 2 (in UniProt coordinates)
    fo.write('\t'.join(ln) + '\n')

for k in P[1]:
    c1,c2 = k
    if c2[1] == 0:  ## DNA or RNA
        if c2[2] == 'D':
            tp = 'D'
        elif c2[2] == 'R':
            tp = 'R'
        else:
            tp = 'L'
    else:
        tp = 'L'
    ln = [tp, PDB, c1, '%s_%d_%s' % c2, chain2up[c1][0], c2[2], 
          '', '', '', '',                       ## genes and genelists
          '',                                   ## reference seq1 length
          '%d' % len(chain2mappedresid[c1]),    ## mapped residues of chain 1
          '%d' % len([i for i in poi[c1][c2][0] if i in pdb2upresmaps[c1]]), ## mapped interface residues of chain 1
          '-',                                  ## reference seq2 length
          '-',                                  ## mapped residues of chain 2
          '-',                                  ## interface residues of chain 2
          '%.2f' % nmut_real[1][(c1,c2)][0],    ## mutations in interface (index is correct for poi)
          '%g' % P[1][(c1,c2)][0],              ## P_mutations in interface (index is correct for poi)
          '',                                   ## q-value(full): stays empty
          '',                                   ## q-value(constr): stays empty
          '%.2f' % nmut_real[1][(c1,c2)][0],    ## mutations in partner 1
          '%g' % P[1][(c1,c2)][0],              ## P_mutations in partner 1
          '-',                                  ## mutations in partner 2
          '-',                                  ## P_mutations in partner 2
          ifres[k][0],                          ## identities of residues in interface of chain 1 (in UniProt coordinates)
          '-',                                  ## identities of residues in interface of chain 2
          ifresdir[k][0],                       ## identities of residues in direct interface of chain 1 (in UniProt coordinates)
          '-',                                  ## identities of residues in direct interface of chain 2
          ifresmut[k][0],                       ## identities of mutations in interface of chain 1 (in UniProt coordinates)
          '-']                                  ## identities of mutations in interface of chain 2
    fo.write('\t'.join(ln) + '\n')

fo.close()
