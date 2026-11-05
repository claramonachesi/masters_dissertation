#!/usr/bin/env python3
# This script is addapted from https://github.com/monteirotorres/docking-scripts/blob/main/python/Consensus.py

import sys
import multiprocessing as mp

#Path dos ligantes (rescores)
sdffile = sys.argv[1] #sdf file
outfile = sys.argv[2] #sdf file

with open(outfile, 'w'):
    pass

rfscores = []
vinascores = []
plecscores = []
with open(sdffile, 'r') as f:
    for line in f:
        if '>  <rfscore_v3>' in line:
            rfscores.append(float(next(f).strip()))
        if '>  <vina_affinity>' in line:
            vinascores.append(float(next(f).strip()))
        if '>  <PLEClinear_p5_l1_s65536>' in line:
            plecscores.append(float(next(f).strip()))

# faixa dos scores
maxrf = max(rfscores)
minrf = min(rfscores)
rfrange = maxrf - minrf

maxvina = max(vinascores)
minvina = min(vinascores)
vinarange = maxvina - minvina

maxplec = max(plecscores)
minplec = min(plecscores)
plecrange = maxplec - minplec

def getnext(file):
    return next(file).strip()

normalizados_rf = []
normalizados_plec = []
normalizados_vina = []
molecules = set() #OBS1
consensus_scores = {}
with open(sdffile, 'r') as f:
    outstring = []
    getmolname = True
    for line in f:
        if getmolname is True:
            getmolname = False
            molname = line.strip()
        if '>  <rfscore_v3>' in line:
            rfscore = float(getnext(f))
            relative_rfscore = (rfscore - minrf) / rfrange
            normalizados_rf.append(relative_rfscore)
            outstring.append(line)
            outstring.append(str(rfscore)+'\n')
            outstring.append('\n')
            outstring.append('>  <relative_rfscore_v3>\n')
            outstring.append(str(relative_rfscore)+'\n')
            continue

        if '>  <PLEClinear_p5_l1_s65536>' in line:
            plecscore = float(getnext(f))
            relative_plecscore = (plecscore - minplec) / plecrange
            normalizados_plec.append(relative_plecscore)
            outstring.append(line)
            outstring.append(str(plecscore)+'\n')
            outstring.append('\n')
            outstring.append('>  <relative_plecscore>\n')
            outstring.append(str(relative_plecscore)+'\n')
            continue

        if '>  <vina_affinity>' in line:
            vinascore = float(getnext(f))
            relative_vinascore = 1 - ((vinascore - minvina) / vinarange)
            normalizados_vina.append(relative_vinascore)
            outstring.append(line)
            outstring.append(str(vinascore)+'\n')
            outstring.append('\n')
            outstring.append('>  <relative_vinascore>\n')
            outstring.append(str(relative_vinascore)+'\n')
            continue

        if '$$$$' not in line:
            outstring.append(line)
            continue
        if '$$$$' in line:
            consensus_vina_rf = (relative_rfscore + relative_vinascore) /2
            consensus_vina_rf_plec = (relative_rfscore + relative_vinascore + relative_plecscore) /3
            outstring.append('>  <relative_vina_affinity>\n')
            outstring.append(str(relative_vinascore)+'\n')
            outstring.append('\n')
            outstring.append('>  <consensus_vina_rf>\n')
            outstring.append(str(consensus_vina_rf)+'\n')
            outstring.append('\n')
            outstring.append('>  <consensus_vina_rf_plec>\n')
            outstring.append(str(consensus_vina_rf_plec)+'\n')
            outstring.append('\n')
            outstring.append('>  <Enamine_Catalog_ID>\n')
            outstring.append(molname.split('/')[0]+'\n')
            outstring.append('\n')
            if molname.startswith('decoy'):
                outstring.append('>  <type>\n')
                outstring.append('decoy\n')
                outstring.append('\n')
            elif molname.startswith('active'):
                outstring.append('>  <type>\n')
                outstring.append('active\n')
                outstring.append('\n')
            outstring.append(line)
            if vinascore < 0:
                with open(outfile, 'a') as outf:
                    outf.write(''.join(outstring))
                    consensus_scores[outstring[0]] = consensus_vina_rf_plec
                    molecules.add(molname.split('_')[0])
            outstring = []
            getmolname = True
            continue
