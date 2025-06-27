#!/usr/bin/env python3
# SPDX-License-Identifier: BSD-3-Clause
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
    """Collect command-line arguments"""

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

    parser.add_argument(
        "-a",
        "--atom_palette",
        help="select the atom palette used by RDKit",
        type=str,
        default="rdkit",
        choices=["avalon", "bw", "cdk", "rdkit"],
    )

    return parser.parse_args()


def record_chopper(files_read):
    """Read input file(s) to return a list of records

    A valid record is a line to contain a SMILES string, a label, and
    a name, each separated by one space.  By

    - label`a`, RDKit attempts to provide a preview
    - label `m`, a preview is to be prepared manually (e.g., ChemDraw)
    - label `#` comments out to skip this entry

    The name of a structure is appended by `.svg` and `.png` for each of the
    preview files.  If present, spaces in the name (e.g., `chiral ligand`)
    are converted into an underscore (`chiral_ligand.png`)."""
    temporary_content = []
    for file_read in files_read:
        temporary_content += file_read.readlines()

    record_list = [i.strip() for i in temporary_content if len(i.split()) >= 3]

    return record_list


def select_atom_palette(atom_palette):
    """Adjust the atom palette to be used by RDKit

    This allows to select one of the four atom palettes presented in
    Greg Landrum's blog post by May 26, 2023, or to implicitly use
    RDKit's default.

    https://greglandrum.github.io/rdkit-blog/posts/2023-05-26-drawing-options-explained.html
    """
    choices = {
        "avalon": "opts.useAvalonAtomPalette()",
        "bw": "opts.useBWAtomPalette()",
        "cdk": "opts.useCDKAtomPalette()",
        "rdkit": "",  # i.e. the implicit default
    }
    selection = ""
    selection = choices.get(atom_palette, None)

    return selection if selection is not None else ""


def svgDepict(mol, colors):
    """Define the more general script parameters

    - disable FreeType rendering to enable bold font-weight
    - this is the version for ligands (e.g., show a red metal circle)
    """
    d2d = rdMolDraw2D.MolDraw2DSVG(256, 256, -1, -1, noFreetype=True)
    opts = d2d.drawOptions()
    opts.explicitMethyl = True
    opts.addStereoAnnotation = True
    opts.atomHighlightsAreCircles = True

    # a safeguard if using RDKit's default atom palette
    if colors:
        eval(colors)

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
    """Define transition metals by their atomic number

    The example of the RDKit Cookbook[1] is extended to assign a dummy
    atom indicated by `*` the atomic number of 0.  To prepare the previews
    about the ligands, this equally is processed as a transition metal."""
    n = atom.GetAtomicNum()
    return (22 <= n <= 29) or (40 <= n <= 47) or (72 <= n <= 79) or (n == 0)


def reset_dative_bonds(mol):
    """Edit some "dative bonds"

    Contrasting to the example in the RDKit Cookbook,[1] and based on a
    discussion in the RDKit user forum[2] about this topic, the previews of
    Avogadro2 use plain single bonds instead of 'dative bonds' (depicted by
    pointy arrows) between ligand(s) and transition metal/dummy atom.

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


def generate_previews(line, colors):
    """Provide each record a .svg and .png preview"""
    smiles = line.split()[0]
    name = "_".join(line.split()[2:])
    print("Running", name)

    mol = Chem.MolFromSmiles(smiles, sanitize=False)
    mol = reset_dative_bonds(mol)
    svg = svgDepict(mol, colors).replace("*", "")

    with open(name + ".svg", "w", encoding="utf-8") as svg_file:
        svg_file.write(svg)

    cairosvg.svg2png(bytestring=svg, write_to=name + ".png")


def main():
    """Join the functionalities"""
    args = get_args()
    record_list = record_chopper(args.file)
    colors = select_atom_palette(args.atom_palette)

    process_manually = []
    process_skipped = []

    for record in record_list:
        try:
            split_record = record.split()
            if split_record[1] == "a":
                generate_previews(record, colors)

            elif split_record[1] == "m":
                process_manually.append(record)

            elif split_record[1] == "#":
                process_skipped.append(record)
            else:
                process_skipped.append(record)
        except OSError as e:
            print(f"error to process {record}, {e}")
        except Exception as e:
            print(f"unexpected error for {record}: {e}")

    if process_manually:
        print("\nentries to be processed manually (label `m`):")
        print(*(entry for entry in process_manually), sep="\n")

    if process_skipped:
        print("\nentries commented out (`#`), or with an unknown label:")
        print(*(entry for entry in process_skipped), sep="\n")


if __name__ == "__main__":
    main()
