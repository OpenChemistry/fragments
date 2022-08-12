#!/usr/bin/env python

import sys
import os

from openbabel import openbabel
from openbabel import pybel


# repeat through all the files on the command-line
# we can change this to use the glob module as well
#  e.g., find all the files in a set of folders
for argument in sys.argv[1:]:
    filename, extension = os.path.splitext(argument)

    # read all the molecules from the supplied file
    for mol in pybel.readfile(extension[1:], argument):
        print("Running", mol.title)

        # change the dummy atom to a metal atom
        for atom in mol.atoms:
            if atom.OBAtom.GetAtomicNum() == 0:
                atom.OBAtom.SetAtomicNum(46) # Pd
        
        # generate coordinates
        pybel._builder.Build(mol.OBMol)
        mol.addh()
        
        ff = pybel._forcefields["uff"]
        success = ff.Setup(mol.OBMol)
        if not success:
            sys.exit("Cannot set up forcefield")

        ff.ConjugateGradients(500, 1.0e-4)

        ff.GetCoordinates(mol.OBMol)

        # change the Pd atom back to a dummy atom
        for atom in mol.atoms:
            if atom.OBAtom.GetAtomicNum() == 46:
                atom.OBAtom.SetAtomicNum(0)

        output = mol.title + ".sdf"

        mol.write("sdf", output, overwrite=True)

