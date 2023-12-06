#!/usr/bin/env python

import sys

from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem.Draw import rdMolDraw2D
from rdkit.Chem import rdDepictor
rdDepictor.SetPreferCoordGen(True)

import cairosvg

def svgDepict(mol):
    # disable FreeType rendering to enable bold font-weight
    # this is the version for ligands (e.g., show a red metal circle)
    d2d = rdMolDraw2D.MolDraw2DSVG(256,256,-1,-1,noFreetype=True)
    opts = d2d.drawOptions()
    opts.explicitMethyl = True
    opts.addStereoAnnotation = True
    opts.atomHighlightsAreCircles = True

    mc = Chem.Mol(mol.ToBinary())
    Chem.Kekulize(mc)
    rdDepictor.Compute2DCoords(mc)

    # position of the "*" atom
    posStar = 0
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 0:
            posStar = atom.GetIdx()

    d2d.DrawMolecule(mc, highlightAtoms=[posStar])
    d2d.FinishDrawing()
    svg = d2d.GetDrawingText()
    return svg

# repeat through all the files on the command-line
# we can change this to use the glob module as well
#  e.g., find all the files in a set of folders

for argument in sys.argv[1:]:
    # each of these files should have a bunch of SMILES
    with open(argument) as smiles_file:
        for line in smiles_file:

            # offer comment lines
            if str(line).startswith("#"):
                continue

            smiles = line.split()[0]
            name = "_".join(line.split()[1:])
            print("Running", name)

            if "*" not in smiles:
                smiles = "*" + smiles

            mol = Chem.MolFromSmiles(smiles, sanitize=False)
            svg = svgDepict(mol).replace("*", "")

            # save the SVG
            with open(name+'.svg', 'w') as svg_file:
                svg_file.write(svg)

            # save a PNG
            cairosvg.svg2png(bytestring=svg, write_to=name+".png")
