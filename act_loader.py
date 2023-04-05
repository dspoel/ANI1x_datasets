#!/usr/bin/env python3

import argparse, os, sys
# ANI-1x code
import dataloader as dl

# ACT stuff
from elements import *
from molprops import *
from get_mol_dict import *
from gaff_to_alexandria import *
from mol_csv_api import *

# Datasets from:
# 1. Original ANI-1x data (https://doi.org/10.1063/1.5023802)
# 2. CHNO portion of the data set used in AIM-Net (https://doi.org/10.1126/sciadv.aav6490)
# 3. The coupled cluster ANI-1ccx data set (https://doi.org/10.1038/s41467-019-10827-4)
datasets = { 1: "wb97x_dz", 2: "wb97x_tz", 3: "ccsd(t)_cbs" }

# TODO: implement atomization energies for DFT.
atomization = { 1: None,
                2: None, 
                3: { "H": -0.5991501324919538, "C": -38.03750806057356, "N": -54.67448347695333, "O": -75.16043537275567 }
                }

def parseArguments():
    parser = argparse.ArgumentParser(description=
    """
    Convert data from the ANI-1x library to molprop compatible with the Alexandria
    Chemistry Toolkit. https://github.com/dspoel/ACT
    """)
    # Path to the ANI-1x data set
    h5file  = "ani1x.h5"
    parser.add_argument("-f",     "--inputfile",  help="Input file for reading, default "+h5file,   type=str,    default=h5file)
    molprop = "ani1x.xml"
    parser.add_argument("-o",     "--outputfile", help="Output file for writing, default "+molprop,  type=str,    default=molprop)
    parser.add_argument("-inchi", "--inchi",      help="Output file containing InChi identifers for all compounds", type=str, default=None)
    dataset = 3
    myhelp  = ""
    for d in datasets:
        myhelp += ( " %d: %s" % ( d, datasets[d] ))
    myhelp += (" default: %d" % dataset)
    parser.add_argument("-set", "--set", help=myhelp, type=int, default=dataset)
    args = parser.parse_args()
    return args

if __name__ == "__main__":
    args = parseArguments()

    # The coupled cluster ANI-1ccx data set (https://doi.org/10.1038/s41467-019-10827-4)
    data_keys = [ datasets[args.set] + ".energy" ]
    have_forces = args.set < 3
    if have_forces:
        data_keys.append(datasets[args.set] + ".forces")

    mp = Molprops()
    mp.open(args.outputfile)
    pretty_print  = True
    g2a           = GaffToAlexandria()
    method, basis = datasets[args.set].split("_")
    
    # https://static-content.springer.com/esm/art%3A10.1038%2Fs41467-019-10827-4/MediaObjects/41467_2019_10827_MOESM1_ESM.pdf
    atom_energy = atomization[args.set]

    # Get the alexandria.csv file
    M = Molecules()
    M.read_default()

    if args.inchi:
        inchi = open(args.inchi, "w")
    # Example for extracting DFT/DZ energies and forces
    counter = 0
    numener = 0
    for data in dl.iter_data_buckets(args.inputfile,keys=data_keys):
        X = data['coordinates']
        Z = data['atomic_numbers']
        E = data[data_keys[0]]
        if len(E) == 0:
            continue
        if have_forces:
            F = data[data_keys[1]]
        elements = []
        for zzz in Z:
            elements.append(AtomNumberToAtomName(zzz))
        MD = MoleculeDict()
        if MD.from_coords_elements(elements, X[0]):
            mymol = M.find_inchi(MD.inchi)
            if mymol and mymol.iupac:
                mp1 = Molprop(mymol.iupac)
            else:
                mp1 = Molprop(MD.inchi)
            if args.inchi:
                inchi.write("%s\n" % MD.inchi)
            for b in MD.bonds.keys():
                mp1.add_bond(b[0], b[1], MD.bonds[b])
            fatoms = []
            for j in range(len(elements)):
                fatoms.append(j+1)
            frag = Fragment(MD.inchi, 0, 1, 1, fatoms, MD.mol_weight, MD.formula)
            mp1.add_fragment(frag)
            jobtype = "Opt"
        
            for i in range(len(E)):
                myexp  = Experiment("Theory", "Smith2019a", "program", method, basis, "sp",
                                    jobtype, args.inputfile, have_forces)
                jobtype = "SP"
                myener = E[i]
                for atom in range(len(MD.atoms)):
                    obtype = g2a.rename(MD.atoms[atom+1]["obtype"])
                    fi = [ 0, 0, 0 ]
                    if have_forces:
                        fi = F[i][atom]
                    myexp.add_atom(elements[atom], obtype, atom+1, "Angstrom",
                                   X[i][atom][0], X[i][atom][1], X[i][atom][2],
                                   "Hartree/Bohr", fi[0], fi[1], fi[2])
                    if elements[atom] in atom_energy:
                        myener -= atom_energy[elements[atom]]
                myexp.add_energy("DeltaE0", "Hartree", 0, "Gas", myener)
                mp1.add_experiment(myexp)
            mp.add_molecule(mp1, pretty_print)
            counter += 1
            numener += len(E)
        else:
            print("Cannot interpret compound %d" % counter)
            print(Z)
    
mp.close()
if args.inchi:
    inchi.close()

print("There are %d compounds and %d energies in the data set %s" % (counter, numener, datakeys[0]))

