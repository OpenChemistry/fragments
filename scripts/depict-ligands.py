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


def is_transition_metal(at):
    """define transition metals by their atomic number

    For the purpose of a motif in the template library of ligands, the
    dummy atom `*` equally should be processed as if it were a transition
    metal.  By convention, its atomic number is 0."""
    n = at.GetAtomicNum()
    return (n>=22 and n<=29) or (n>=40 and n<=47) or (n>=72 and n<=79) or (n==0)


def set_dative_bonds(mol, fromAtoms=(7,8)):
    """ convert some bonds to dative

    Replaces some single bonds between metals and atoms with atomic numbers in fomAtoms
    with dative bonds. The replacement is only done if the atom has "too many" bonds.

    Returns the modified molecule.

    """
    pt = Chem.GetPeriodicTable()
    rwmol = Chem.RWMol(mol)
    rwmol.UpdatePropertyCache(strict=False)
    metals = [at for at in rwmol.GetAtoms() if is_transition_metal(at)]
    for metal in metals:
        for nbr in metal.GetNeighbors():
            if nbr.GetAtomicNum() in fromAtoms and \
               nbr.GetExplicitValence()>pt.GetDefaultValence(nbr.GetAtomicNum()) and \
               rwmol.GetBondBetweenAtoms(nbr.GetIdx(),metal.GetIdx()).GetBondType() == Chem.BondType.SINGLE:
                rwmol.RemoveBond(nbr.GetIdx(),metal.GetIdx())
                rwmol.AddBond(nbr.GetIdx(),metal.GetIdx(),Chem.BondType.DATIVE)
    return rwmol


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
            mol = set_dative_bonds(mol)
            svg = svgDepict(mol).replace("*", "")

            # save the SVG
            with open(name+'.svg', 'w') as svg_file:
                svg_file.write(svg)

            # save a PNG
            cairosvg.svg2png(bytestring=svg, write_to=name+".png")
