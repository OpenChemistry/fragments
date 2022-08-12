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
    d2d = rdMolDraw2D.MolDraw2DSVG(256,256,-1,-1,True)
    opts = d2d.drawOptions()
    opts.dummiesAreAttachments = True
    opts.explicitMethyl = True
    opts.addStereoAnnotation = True
    d2d.DrawMolecule(mol)
    d2d.FinishDrawing()
    svg = d2d.GetDrawingText()
    svg.replace("font-weight:normal", "font-weight:bold")
    svg.replace("font-size:16px", "font-size:20px")
    return svg

# repeat through all the files on the command-line
# we can change this to use the glob module as well
#  e.g., find all the files in a set of folders
embed = AllChem.ETKDGv3()                                                                                                        
embed.useRandomCoords = True
embed.useSmallRingTorsions = True

for argument in sys.argv[1:]:
    # each of these files should have a bunch of SMILES
    with open(argument) as smiles_file:
        for line in smiles_file:
            smiles = line.split()[0]
            name = "_".join(line.split()[1:])
            print("Running", name)
            # add an attachment center if needed
            if not smiles.startswith('*'):
                smiles = '*' + smiles
            mol = Chem.MolFromSmiles(smiles)
            svg = svgDepict(mol)

            # save the SVG
            with open(name+'.svg', 'w') as svg_file:
                svg_file.write(svg)

            # save a PNG
            cairosvg.svg2png( bytestring=svg, write_to=name+".png" )

            # save an SDF with coordinates
            m2 = Chem.AddHs(mol)
            AllChem.EmbedMolecule(m2, params=embed)
            AllChem.MMFFOptimizeMoleculeConfs(m2, maxIters=2000)
            with open(name + ".mol", 'w') as mol_file:
                mol_file.write(Chem.MolToMolBlock(m2))