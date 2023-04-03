#!/usr/bin/env python3

import os, sys        
import dataloader as dl

# ACT stuff
from elements import *
from molprops import *
from get_mol_dict import *
from gaff_to_alexandria import *

# Path to the ANI-1x data set
path_to_h5file = 'ani1x.h5'

# List of keys to point to requested data
#data_keys = ['wb97x_dz.energy','wb97x_dz.forces'] # Original ANI-1x data (https://doi.org/10.1063/1.5023802)
#data_keys = ['wb97x_tz.energy','wb97x_tz.forces'] # CHNO portion of the data set used in AIM-Net (https://doi.org/10.1126/sciadv.aav6490)

data_keys = ['ccsd(t)_cbs.energy'] # The coupled cluster ANI-1ccx data set (https://doi.org/10.1038/s41467-019-10827-4)

#data_keys = ['wb97x_dz.dipoles'] # A subset of this data was used for training the ACA charge model (https://doi.org/10.1021/acs.jpclett.8b01939)

mp = Molprops()
mp.open("ani1x.xml")
pretty_print = True
g2a  = GaffToAlexandria()

# https://static-content.springer.com/esm/art%3A10.1038%2Fs41467-019-10827-4/MediaObjects/41467_2019_10827_MOESM1_ESM.pdf
atom_energy = { "H": -0.5991501324919538, "C": -38.03750806057356, "N": -54.67448347695333, "O": -75.16043537275567 }

# Example for extracting DFT/DZ energies and forces
counter = 0
numener = 0
for data in dl.iter_data_buckets(path_to_h5file,keys=data_keys):
    X = data['coordinates']
    Z = data['atomic_numbers']
    E = data['ccsd(t)_cbs.energy']
    elements = []
    for zzz in Z:
        elements.append(AtomNumberToAtomName(zzz))
    mp1      = None
    firstmol = True
    for i in range(len(E)):
        MD = MoleculeDict()
        if MD.from_coords_elements(elements, X[i]):
            jobtype = "SP"
            if firstmol:
                mp1 = Molprop(MD.inchi)
                for b in MD.bonds.keys():
                    mp1.add_bond(b[0], b[1], MD.bonds[b])
                fatoms = []
                for j in range(len(elements)):
                    fatoms.append(j+1)
                frag = Fragment(MD.inchi, 0, 1, 1, fatoms, MD.mol_weight, MD.formula)
                mp1.add_fragment(frag)
                jobtype = "Opt"
                firstmol = False
            myexp  = Experiment("Theory", "smith2019a", "program", "ccsd(t)", "cbs", "sp", jobtype, path_to_h5file, False)
            myener = E[i]
            for atom in range(len(MD.atoms)):
                obtype = g2a.rename(MD.atoms[atom+1]["obtype"])
                myexp.add_atom(elements[atom], obtype, atom+1, "Angstrom",
                               X[i][atom][0], X[i][atom][1], X[i][atom][2],
                               "Hartree/Bohr", 0, 0, 0)
                if elements[atom] in atom_energy:
                    myener -= atom_energy[elements[atom]]
            myexp.add_energy("DeltaE0", "Hartree", 0, "Gas", myener)
            mp1.add_experiment(myexp)
    if not firstmol:
        mp.add_molecule(mp1, pretty_print)
    counter += 1
    numener += len(E)
    if counter == 10:
        break
    
print("There are %d compounds and %d energies in the data set %s" % (counter, numener, data_keys[0]))

mp.close()
