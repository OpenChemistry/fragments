#!/usr/bin/env python

import sys
import glob

from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem.Draw import rdMolDraw2D
from rdkit.Chem import rdDepictor
rdDepictor.SetPreferCoordGen(True)

import cairosvg

def svgDepict(mol):
    # disable FreeType rendering to enable bold font-weight
    d2d = rdMolDraw2D.MolDraw2DSVG(256,256,-1,-1,noFreetype=True)
    opts = d2d.drawOptions()
    opts.dummiesAreAttachments = True
    opts.explicitMethyl = True
    opts.addStereoAnnotation = True

    mc = Chem.Mol(mol.ToBinary())
    Chem.Kekulize(mc)
    rdDepictor.Compute2DCoords(mc)

    d2d.DrawMolecule(mc)
    d2d.FinishDrawing()
    svg = d2d.GetDrawingText()
    return svg

# repeat through all the files on the command-line
# we can change this to use the glob module as well
#  e.g., find all the files in a set of folders
for filename in glob.iglob("*.smi"):

    smiles = ''
    with open(filename) as f:
        smiles = f.read().strip()
        smiles = '*' + smiles + '*'

    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        print("Error reading SMILES: %s" % smiles, filename)
        continue
    svg = svgDepict(mol)
    svg = svg.replace("font-weight:normal", "font-weight:bold").replace("font-size:16px", "font-size:20px")

    # save the SVG
    name = filename[:-4]
    with open(name+'.svg', 'w') as svg_file:
        svg_file.write(svg)

    # save a PNG
    cairosvg.svg2png(bytestring=svg, scale=2.0, write_to=name+".png")
