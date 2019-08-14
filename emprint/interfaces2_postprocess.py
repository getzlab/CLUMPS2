import os
import random
import re
import rpy2.robjects as R
import urllib2

from gzip import GzipFile


## wget http://ligand-expo.rcsb.org/dictionaries/cc-counts-extra.tdd

RDIR = '../res/interfaces-PanCan.pc14mc3.run1/'
TESTED = '../dat/pdbs_forIFan_all'

MAPFILE = '../res/huniprot2pdb.run18.txt.gz'

AVOIDLIGANDS = set(['SO4', 'GOL', 'CL', 'EDO', 'CA', 'NAG', 'NA', 'PO4', 'ACT', 'PEG',
                    'K', 'NDG', 'IOD', 'MPD', 'PG4', 'BME', 'DMS', 'FUC', 'CIT', 'EPE',
                    'ACY', 'MES', 'TRS', 'MAN', 'FMT', 'PGE', 'IPA', 'NO3', 'FLC', '1PE',
                    'CO', 'BMA', 'IMD', 'DTT', 'SCN', 'GAL', 'BR', 'BGC', 'MLI', 'GLC',
                    'CAC', 'P6G', 'OLA', 'AKG', 'MRD', 'OGA'])
## comes from the following top degrees: [('SO4', 1040), ('GOL', 840), ('ZN', 676), ('CL', 667), ('EDO', 520), ('MG', 494), ('CA', 402), ('NAG', 394), ('DNA/RNA', 356), ('NA', 349), ('PO4', 345), ('ACT', 250), ('UNX', 212), ('DNA', 210), ('RNA', 157), ('PEG', 152), ('ADP', 138), ('K', 112), ('MN', 101), ('NI', 99), ('GDP', 85), ('NDG', 82), ('IOD', 81), ('ANP', 81), ('MPD', 80), ('PG4', 76), ('BME', 75), ('DMS', 74), ('FUC', 74), ('CIT', 73), ('ATP', 71), ('EPE', 70), ('ACY', 69), ('MES', 67), ('TRS', 66), ('MAN', 62), ('FMT', 56), ('CD', 55), ('PGE', 53), ('IPA', 52), ('NAD', 50), ('HEM', 49), ('NO3', 48), ('FLC', 46), ('1PE', 44), ('NAP', 43), ('CO', 43), ('BMA', 41), ('GNP', 41), ('AMP', 41), ('IMD', 40), ('SAH', 39), ('DTT', 36), ('GSH', 35), ('HG', 34), ('SCN', 34), ('GAL', 33), ('STU', 32), ('BR', 32), ('SAM', 32), ('FE', 31), ('BGC', 30), ('FAD', 29), ('MLI', 29), ('GTP', 28), ('TLA', 27), ('GLC', 26), ('CAC', 25), ('COA', 25), ('NDP', 25), ('PLP', 25), ('P6G', 24), ('FE2', 22), ('OLA', 22), ('AKG', 22), ('MRD', 21), ('CU', 21), ('OGA', 20), ('DIO', 19), ('GLY', 19)]

fl = os.listdir(RDIR)

dat = {'PPI':{},
       'PLI':{},
       'PdI':{},
       'PDI':{},
       'PRI':{}}

testedstruct = {}
for l in file(TESTED).readlines():
    l = l.strip().split()
    if len(l) > 1:
        add = map(lambda x:x.replace(':','_'), sorted(l[1].split('_')))
        n = '-'.join([l[0]] + add)
    else:
        n = l[0]
    testedstruct[n] = 0
        
print len(testedstruct), 'structures should have been tested.'

cnt = 0
for fn in fl:
    cnt += 1
    if not cnt % 100:
        print cnt
    testedstruct[fn] = 1
    fi = file(RDIR + fn)
    while 1:
        l = fi.readline()
        if not l:
            break
        l = l.strip('\n').split('\t')
        ## fix the mutation masses
        m1 = 0
        m2 = 0
        if l[28] and l[28] != '-':
            for i in l[28].split(';'):
                m1 += float(i.split('_')[1])
        if l[29] and l[29] != '-':
            for i in l[29].split(';'):
                m2 += float(i.split('_')[1])
        l[16] = '%.2f' % (m1+m2)
        l[20] = '%.2f' % (m1)
        l[22] = '%.2f' % (m2)
        if l[0] == 'P':
            if l[4] == l[5]:  ## homodimers not used for now
                continue
            ab = tuple(sorted(l[4:6]))
            for i in [11,12,14,15 ,17]:  ## covered length, interface residues and p_interface are float
                l[i] = float(l[i])
            if ab not in dat['PPI']:
                dat['PPI'][ab] = []
            dat['PPI'][ab].append(l)
        else:
            if (l[3].endswith('_0_DNA') or l[3].endswith('_0_DNA/RNA')) and ',' not in l[3]:
                typ = 'd'
            elif l[3].endswith('_0_DNA/RNA') and ',' in l[3]:
                typ = 'D'
            elif l[3].endswith('_0_RNA'):
                typ = 'R'
            else:
                typ = 'L'
                if l[5] in AVOIDLIGANDS:
                    continue
            di = dat['P%sI' % typ]
            if l[4] not in di:
                di[l[4]] = []
            di[l[4]].append(l)

print 'this number should match:', len(testedstruct)
nores = [i for i in testedstruct if not testedstruct[i]]
c = len(nores)
print 'number of structures that were not tested', c
if c < 20:
    print 'namely:', ' '.join(map(str,nores))
##### 
#random.shuffle(nores)
#for i in nores:
#    print i
    

## get gene names
up2gn = {}
up2eg = {}
up2len = {}
fi = GzipFile('../dat/uniprot.sprot.human.genename_entrez_length.gz')
x = fi.readline()
while 1:
    l = fi.readline()
    if not l:
        break
    l = l.strip('\n').split('\t')
    u = l[0]
    if l[1]:
        up2gn[u] = l[1].split(' ',1)[0].rstrip(';')
    up2len[u] = int(l[3])
    if l[2]:
        up2eg[u] = map(int,l[2].strip(';').split(';'))

up2len['P08107'] = 641

## get cancer gene list
eg2can = {int(i[0]):i[4] for i in map(lambda x:x.strip().split('\t'), file('../dat/allCancerGenes.txt').readlines()[1:])}

up2can = {}
for u in up2eg:
    for e in up2eg[u]:
        if e in eg2can and u not in up2can:
            up2can[u] = eg2can[e]

## get original upids (u2) and %identity of all structures
pdbch2u2 = {}
fi = GzipFile(MAPFILE)
while 1:
    l = fi.readline()
    if not l:
        break
    l = l.split('\t',4)
    iden = 100.0
    if l[3] != '-':
        for s in l[3].split():
            if s.startswith('pdb_identity:'):
                iden = float(s.split(':')[1])
    pdbch2u2[(l[0], l[2])] = [l[1], iden]

## do selection of instances

## PPI:  protein length and interface size, then take best p-value
out = []
for ab in dat['PPI']:
    ll = dat['PPI'][ab]
    print ll
    maxLA = float(max(map(lambda x:x[11], ll)))
    maxIA = float(max(map(lambda x:x[12], ll)))
    maxLB = float(max(map(lambda x:x[14], ll)))
    maxIB = float(max(map(lambda x:x[15], ll)))
    if min(maxIA, maxIB) == 0:  ## e.g. 1ymm, where the interacting residue (pdb:104) is before the protein start (pdb:119 -> up:1)
        continue
    for vv in ll:
        if vv[4] == ab[0]:  ## order of the proteins (could be different than in ab)
            sim1 = pdbch2u2[(vv[4], vv[1]+'-'+vv[2])][1]
            sim2 = pdbch2u2[(vv[5], vv[1]+'-'+vv[3])][1]
            score = sim1 * (float(vv[12])/maxIA + float(vv[11])/maxLA) + sim2 * (float(vv[15])/maxIB +  + float(vv[14])/maxLB)
        else:
            sim1 = pdbch2u2[(vv[5], vv[1]+'-'+vv[3])][1]
            sim2 = pdbch2u2[(vv[4], vv[1]+'-'+vv[2])][1]
            score = sim1 * (float(vv[15])/maxIA + float(vv[14])/maxLA) + sim2 * (float(vv[12])/maxIB + float(vv[11])/maxLB)
        vv.append(score)
    ll.sort(key=lambda x:x[-1], reverse=True)
    ties = [ll[0]]
    s = ll[0][-1]
    for vv in ll[1:]:
        if vv[-1] >= 0.8 * s:
            ties.append(vv)
        else:
            break
    ties.sort(key=lambda x:x[17])  ## sort by p
    sel = ties[0]  ## best p
    #if float(sel[16]) < 0.01:  ## no mutation mass at the interface
    #    continue
    out.append(sel)
    
##########################################################

if 0:
    up_ord = prot2ifclusters.keys()
    pvals = []
    for u in up_ord:
        ifc = prot2ifclusters[u]
        for i in xrange(len(ifc)):
            pvals.append(ifc[i][0][1])

    qvals = list(R.r['p.adjust'](R.FloatVector(pvals), "fdr"))

    idx = 0
    for u in up_ord:
        ifc = prot2ifclusters[u]
        for i in xrange(len(ifc)):
            ifc[i][0][2][18+ifc[i][0][3]] = '%g' % qvals[idx] 
            idx += 1

else:
    pvals = map(lambda x:float(x[17]), out)
    qvals = list(R.r['p.adjust'](R.FloatVector(pvals), "fdr"))
    for i in xrange(len(out)):
        sel = out[i]
        sel[18] = '%g' % qvals[i]
        sel[19] = '-'

for i in xrange(len(out)):
    sel = out[i]
    if sel[4] in up2gn:
        sel[6] = up2gn[sel[4]]
    if sel[4] in up2can:
        sel[8] = up2can[sel[4]]
    if sel[5] in up2gn:
        sel[7] = up2gn[sel[5]]
    if sel[5] in up2can:
        sel[9] = up2can[sel[5]]
    if sel[4] in up2len:
        sel[10] = up2len[sel[4]]
    if sel[5] in up2len:
        sel[13] = up2len[sel[5]]
    sel[11] = int(sel[11])
    sel[12] = int(sel[12])
    sel[14] = int(sel[14])
    sel[15] = int(sel[15])
    sel += pdbch2u2[(sel[4], sel[1]+'-'+sel[2])]
    sel += pdbch2u2[(sel[5], sel[1]+'-'+sel[3])]

## restricted analysis q-values
canidx = [i for i in xrange(len(out)) if out[i][8] or out[i][9]]
qvals = list(R.r['p.adjust'](R.FloatVector([out[i][17] for i in canidx]), "fdr"))
for i in xrange(len(canidx)):
    out[canidx[i]][19] = '%g' % qvals[i]

fo = file('../res/' + RDIR.strip('/').split('/')[-1] + '.PPI.csv', 'w')
for sel in out:
    fo.write('\t'.join(map(str, sel)) + '\n')

fo.close()

##  do a per-protein ranking
prot2lines = {}
for i in xrange(len(out)):
    sel = out[i]
    for j in [0,1]:
        if sel[4+j] not in prot2lines:
            prot2lines[sel[4+j]] = []
        prot2lines[sel[4+j]].append((sel[21+2*j], j, i, sel[8+j]))

x = map(lambda x:x.sort(key=lambda x:float(x[0])), prot2lines.values())
prot2lines = prot2lines.items()
prot2lines.sort(key=lambda x:float(x[1][0][0]))
qq = list(R.r['p.adjust'](R.FloatVector([x[1][0][0] for x in prot2lines]), "fdr"))
rqq = list(R.r['p.adjust'](R.FloatVector([x[1][0][0] for x in prot2lines if x[1][0][3]]), "fdr"))
p2restrq = dict(zip([x[0] for x in prot2lines if x[1][0][3]], rqq))

fo = file('../res/' + RDIR.strip('/').split('/')[-1] + '.PPI_perProtein.csv', 'w')
for i in xrange(len(prot2lines)):
    q = '%g' % qq[i]
    u,aa = prot2lines[i]
    if u in p2restrq:
        rq = '%g' % p2restrq[u]
    else:
        rq = '-'
    ll = out[aa[0][2]]
    nl = [u, ll[6+aa[0][1]], ll[20+2*aa[0][1]], aa[0][0], q, rq]
    fo.write('\t'.join(nl) + '\n')

fo.close()

#########################################################################

ligNames = {l[0]:l[2] for l in map(lambda x:x.strip('\n').split('\t'), file('../dat/cc-counts-extra.tdd').readlines())}

## POI:  combine nested interfaces; p-value of biggest interface
pdbsumPO = {}
for itp in ['L','d','D','R']:
    di = dat['P%sI' % itp]
    upid2ifaces = {}
    for upid in di:
        if not upid2ifaces.has_key(upid):
            upid2ifaces[upid] = []
        vv = di[upid]
        for i in vv:
            if not i[24]:
                continue
            i[24] = set(map(int,i[24].split(';')))
        vv.sort(key=lambda x:len(x[24]), reverse=True)
        for v in vv:
            #if float(v[16]) < 0.01:  ## no mutation mass at the interface
            #    continue
            newif = True
            for o in upid2ifaces[upid]:
                match = 0
                for i in v[24]:
                    if i in o[0][24]:
                        match += 1
                if match >= 0.8*len(v[24]):
                    newif = False
                    o[1].add(v[1])
                    o[2].add(v[5])
                    break
            if newif:
                upid2ifaces[upid].append([v,set([v[1]]), set([v[5]])])
    out = []
    for upid in upid2ifaces:
        for iface in upid2ifaces[upid]:
            ln = iface[0]
            ln[24] = ';'.join(map(str, sorted(list(ln[24]))))
            if len(iface[1]) > 1:
                ln[1] = ln[1] + '; ' + ','.join(sorted([i for i in iface[1] if i != ln[1]]))
            if len(iface[2]) > 1:
                ln[5] = ln[5] + '; ' + ','.join(sorted([i for i in iface[2] if i != ln[5]]))
            out.append(ln)
    for l in out:
        if l[4] in up2can:
            l[8] = up2can[l[4]]
        l[9] = '-'
        if l[4] not in up2len:
            print '!!! l[4] not in up2len', l[4]
            l[10] = 'NA'
        else:
            l[10] = up2len[l[4]]
        allmolids = []
        ll = l[5].split('; ')
        allmolids.append(ll[0])
        if len(ll) > 1:
            allmolids += ll[1].split(',')
        l += pdbch2u2[(l[4], l[1].split(';')[0]+'-'+l[2])]
        if itp == 'L':
            lnm = []
            for i in allmolids:
                if i not in ligNames:
                    print i, 'missing in ligNames, retrieving from PDB.'
                    html = urllib2.urlopen('http://www.rcsb.org/pdb/files/ligand/%s.cif' % i).read()
                    name = i
                    for x in html.split('\n'):
                        if x.startswith('_chem_comp.name '):
                            name = x.split(' ',1)[1].strip()
                            break
                    ligNames[i] = name
                lnm.append(ligNames[i])
            l[7] = lnm[0]
            if len(lnm) > 1:
                l[7] += '; ' + ','.join(lnm[1:])
        else:
            l[7] = l[3].rsplit('_',1)[-1]
    ## q-values 
    out.sort(key=lambda x:float(x[17]))
    qvals = list(R.r['p.adjust'](R.FloatVector([float(x[17]) for x in out]), "fdr"))
    for i in xrange(len(out)):
        sel = out[i]
        sel[18] = '%g' % qvals[i]
        sel[19] = '-'
    canidx = [i for i in xrange(len(out)) if out[i][8]]
    qvals = list(R.r['p.adjust'](R.FloatVector([float(out[i][17]) for i in canidx]), "fdr"))
    for i in xrange(len(canidx)):
        out[canidx[i]][19] = '%g' % qvals[i]
    fo = file('../res/' + RDIR.strip('/').split('/')[-1] + '.P%sI.csv' % itp, 'w')
    for sel in out:
        if sel[4] in up2gn:
            sel[6] = up2gn[sel[4]]
        fo.write('\t'.join(map(str, sel)) + '\n')
    fo.close()
    #raise Exception('Done.')
