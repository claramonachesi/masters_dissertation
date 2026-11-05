#!/usr/bin/env python
# Counts how many of each interaction type from PLIP (Protein-Ligand Interaction Profiler)

import sys
import re

reportfile = sys.argv[1]
outfile_txt = sys.argv[2]
outfile_tsv = sys.argv[3]

molecule_name = re.split('/', reportfile)[0]

with open(outfile_txt, 'a'):
    pass
with open(outfile_tsv, 'a'):
    pass

## counting the interactions
with open(reportfile, 'r') as file:
    lines = [l.strip() for l in file.readlines()]
    hydrophobic_interactions = False
    hydrogen_bonds = False
    pi_stacking = False
    pi_cation = False
    salt_bridges = False
    water_bridges = False
    halogen_bonds = False
    metal_complexes = False

    count_hydrophobic = -1
    count_hydrogen = -1 #-1 because of the header of the table
    count_pi_stacking = -1
    count_pi_cation = -1
    count_salt_bridges = -1
    count_water_bridges = -1
    count_halogen = -1
    count_metal = -1

    for line in lines:
        if "**Hydrophobic Interactions**" in line:
            hydrophobic_interactions = True
        elif "**Hydrogen Bonds**" in line:
            hydrophobic_interactions = False
            hydrogen_bonds = True
        elif "**pi-Stacking**" in line:
            hydrophobic_interactions = False
            hydrogen_bonds = False
            pi_stacking = True
        elif "**pi-Cation Interactions**" in line:
            hydrophobic_interactions = False
            hydrogen_bonds = False
            pi_stacking = False
            pi_cation = True
        elif "**Salt Bridges**" in line:
            hydrophobic_interactions = False
            hydrogen_bonds = False
            pi_stacking = False
            pi_cation = False
            salt_bridges = True
        elif "**Water Bridges**" in line:
            hydrophobic_interactions = False
            hydrogen_bonds = False
            pi_stacking = False
            pi_cation = False
            salt_bridges = False
            water_bridges = True
        elif "**Halogen Bonds**" in line:
            hydrophobic_interactions = False
            hydrogen_bonds = False
            pi_stacking = False
            pi_cation = False
            salt_bridges = False
            water_bridges = False
            halogen_bonds = True
        elif "**Metal Complexes**" in line:
            hydrophobic_interactions = False
            hydrogen_bonds = False
            pi_stacking = False
            pi_cation = False
            salt_bridges = False
            water_bridges = False
            halogen_bonds = False
            metal_complexes = True

        elif hydrophobic_interactions:
            if "|" in line:
                count_hydrophobic += 1
            elif "+" in line:
                next_item = iter(line)
                next(next_item)
        elif hydrogen_bonds:
            if "|" in line:
                count_hydrogen += 1
            elif "+" in line:
                next_item = iter(line)
                next(next_item)
        elif pi_stacking:
            if "|" in line:
                count_pi_stacking += 1
            elif "+" in line:
                next_item = iter(line)
                next(next_item)
        elif pi_cation:
            if "|" in line:
                count_pi_cation += 1
            elif "+" in line:
                next_item = iter(line)
                next(next_item)
        elif salt_bridges:
            if "|" in line:
                count_salt_bridges += 1
            elif "+" in line:
                next_item = iter(line)
                next(next_item)
        elif water_bridges:
            if "|" in line:
                count_water_bridges += 1
            elif "+" in line:
                next_item = iter(line)
                next(next_item)
        elif halogen_bonds:
            if "|" in line:
                count_halogen += 1
            elif "+" in line:
                next_item = iter(line)
                next(next_item)
        elif metal_complexes:
            if "|" in line:
                count_metal += 1
            elif "+" in line:
                next_item = iter(line)
                next(next_item)

## writing the interactions in a tsv
with open(outfile_tsv, 'w') as f:
    f.write("Molecule Name\tHydrophobic Interactions\tHydrogen Bonds\tPi-Stacking\tPi-Cation Interactions\tSalt Bridges\tWater Bridges\tHalogen Bonds\tMetal Complexes\n")
    f.write(f"{molecule_name}\t")
    if count_hydrophobic != -1:
        f.write(f"{count_hydrophobic}\t")
    if count_hydrophobic == -1:
        f.write("\t")
    if count_hydrogen != -1:
        f.write(f"{count_hydrogen}\t")
    if count_hydrogen == -1:
        f.write("\t")
    if count_pi_stacking != -1:
        f.write(f"{count_pi_stacking}\t")
    if count_pi_stacking == -1:
        f.write("\t")
    if count_pi_cation != -1:
        f.write(f"{count_pi_cation}\t")
    if count_pi_cation == -1:
        f.write("\t")
    if count_salt_bridges != -1:
        f.write(f"{count_salt_bridges}\t")
    if count_salt_bridges == -1:
        f.write("\t")
    if count_water_bridges != -1:
        f.write(f"{count_water_bridges}\t")
    if count_water_bridges == -1:
        f.write("\t")
    if count_halogen != -1:
        f.write(f"{count_halogen}\t")
    if count_halogen == -1:
        f.write("\t")
    if count_metal != -1:
        f.write(f"{count_metal}")
    if count_metal == -1:
        f.write("\t")

