#!/usr/bin/env python3
"""This script provides .svg and .png ligand template previews of Avogadro."""

import argparse

import cairosvg

from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem.Draw import rdMolDraw2D
from rdkit.Chem import rdDepictor
from rdkit.Chem import rdCIPLabeler

rdDepictor.SetPreferCoordGen(True)


def get_args():
    """Get all command-line arguments"""

    parser = argparse.ArgumentParser(
        description="write .svg and .png ligand template previews for Avogadro",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )

    parser.add_argument(
        "file",
        help="one or multiple input file(s) to process",
        metavar="FILE",
        nargs="+",
        type=argparse.FileType("rt"),
    )

    return parser.parse_args()


def record_chopper(files_read):
    """read input file(s) to return a list of records

    A valid report is a line to contain a SMILES string, a label, and
    the (chemical) name (later used to build the file name of the .svg
    and .png preview); each separated from each other by white space."""
    temporary_content = []
    for file_read in files_read:
        temporary_content += file_read.readlines()

    record_list = [i.strip() for i in temporary_content if len(i.split()) >= 3]

    return record_list


def svgDepict(mol):
    """define the more general script parameters

    - disable FreeType rendering to enable bold font-weight
    - this is the version for ligands (e.g., show a red metal circle)
    """
    d2d = rdMolDraw2D.MolDraw2DSVG(256, 256, -1, -1, noFreetype=True)
    opts = d2d.drawOptions()
    opts.explicitMethyl = True
    opts.addStereoAnnotation = True
    opts.atomHighlightsAreCircles = True

    metal_complex = Chem.Mol(mol.ToBinary())
    Chem.Kekulize(metal_complex)
    rdDepictor.Compute2DCoords(metal_complex)
    rdCIPLabeler.AssignCIPLabels(metal_complex)

    # position of the "*" atom
    posStar = 0
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 0:
            posStar = atom.GetIdx()

    d2d.DrawMolecule(metal_complex, highlightAtoms=[posStar])
    d2d.FinishDrawing()
    svg = d2d.GetDrawingText()
    return svg


def is_transition_metal(atom):
    """define transition metals by their atomic number

    For the purpose of a motif in the template library of ligands, the
    dummy atom indicated by `*` equally should be processed as if it
    were a transition metal.  By convention, its atomic number is 0.
    The pattern initially seen in RDKit's cookbook was edited."""
    n = atom.GetAtomicNum()
    return (22 <= n <= 29) or (40 <= n <= 47) or (72 <= n <= 79) or (n == 0)


def reset_dative_bonds(mol):
    """edit some "dative bonds"

    Bonds between atoms of transition metals and typical Lewis base /
    donor atoms will be marked as single bonds.  Initially the function
    by the RDKit Cookbook[1] depicting an example with pointy arrows, a
    subsequent discussion in RDKit's user forum[2] convinced nbehrnd to
    drop this approach in favor of plain bonds.

    [1] http://rdkit.org/docs/Cookbook.html#organometallics-with-dative-bonds
    [2] https://github.com/rdkit/rdkit/discussions/6995

    Returns the modified molecule.
    """
    from_atoms = [6, 7, 8, 15, 16]  # i.e., C, N, O, P, S
    pt = Chem.GetPeriodicTable()
    rwmol = Chem.RWMol(mol)
    rwmol.UpdatePropertyCache(strict=False)
    metals = [at for at in rwmol.GetAtoms() if is_transition_metal(at)]
    for metal in metals:
        for nbr in metal.GetNeighbors():
            if (
                nbr.GetAtomicNum() in from_atoms
                and rwmol.GetBondBetweenAtoms(
                    nbr.GetIdx(), metal.GetIdx()
                ).GetBondType()
                == Chem.BondType.SINGLE
            ):
                rwmol.RemoveBond(nbr.GetIdx(), metal.GetIdx())
                rwmol.AddBond(nbr.GetIdx(), metal.GetIdx(), Chem.BondType.SINGLE)
    return rwmol


def generate_previews(line):
    """provide per record a .svg and .png preview"""
    smiles = line.split()[0]
    name = "_".join(line.split()[2:])
    print("Running", name)

    mol = Chem.MolFromSmiles(smiles, sanitize=False)
    mol = reset_dative_bonds(mol)
    svg = svgDepict(mol).replace("*", "")

    with open(name + ".svg", "w", encoding="utf-8") as svg_file:
        svg_file.write(svg)

    cairosvg.svg2png(bytestring=svg, write_to=name + ".png")


def main():
    """join the functionalities"""
    args = get_args()
    record_list = record_chopper(args.file)
    process_manually = []
    process_skipped = []

    for record in record_list:
        try:
            split_record = record.split()
            if split_record[1] == "a":
                generate_previews(record)

            elif split_record[1] == "m":
                process_manually.append(record)

            elif split_record[1] == "#":
                process_skipped.append(record)
            else:
                process_skipped.append(record)
        except OSError:
            print(f"error to process {record}")
        except Exception as e:
            print(f"non-anticipated error: {e}")

    if process_manually:
        print("\nentries to be processed manually (label `m`):")
        print(*(entry for entry in process_manually), sep="\n")

    if process_skipped:
        print("\nentries commented out (`#`), or with an unknown label:")
        print(*(entry for entry in process_skipped), sep="\n")


if __name__ == "__main__":
    main()
