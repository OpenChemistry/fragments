#!/usr/bin/env python3
# SPDX-License-Identifier: BSD-3-Clause
"""Generate the amino acid z-matrices used by Avogadro's peptide builder.

Atom names, elements, bonds, formal charges and ideal coordinates all come
from the PDB Chemical Component Dictionary (CCD), so the names match the
residue bond dictionary in avogadrolibs (avogadro/core/residuedata.h).

Each row of a .zmat file is

    Element Name Length Angle Dihedral RefA RefB RefC

with 1-based references (0 = none), angles in degrees, and the dihedral as
the IUPAC signed torsion Atom-RefA-RefB-RefC.  There are no bonds or charges
in the file: the builder takes bonds from the residue dictionary and makes
any nitrogen with four bonds (bond orders summed) +1.

The peptide builder depends on the row layout, which this script enforces:

- rows 1-4 are N, CA, C, O; O is placed from (C, CA, N) because the
  builder resets its dihedral from psi
- the backbone H is placed from (N, CA, C) because the builder resets its
  angle and dihedral from phi
- every atom bonded to CA is placed from (CA, N, C), which preserves the
  chirality when the builder rebuilds the backbone from phi/psi/omega
- side-chain atoms never reference O, H, OXT or HXT, which the builder
  moves or drops
- H2 is left out (the builder adds the N-terminal hydrogens); the other
  leaving atoms (OXT, HXT, and PRO's H) come last, because the builder
  skips them for interior residues

Usage:

    python scripts/getamino.py              # regenerate amino/*.zmat
    python scripts/getamino.py TRP ARG      # just these residues
    python scripts/getamino.py --check      # validate the existing files
"""

import argparse
import math
import os
import shlex
import sys
import urllib.request

CCD_URL = "https://files.rcsb.org/ligands/download/{}.cif"

RESIDUES = [
    # the 20 common amino acids
    "ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS", "ILE",
    "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL",
    # selenocysteine and pyrrolysine (one-letter codes U and O)
    "SEC", "PYL",
    # alpha-aminoisobutyric acid
    "AIB",
]  # fmt: skip

BACKBONE = ("N", "CA", "C", "O")
# atoms the builder moves or drops, so nothing else may be placed from them
NOT_A_REFERENCE = {"O", "H", "H2", "OXT", "HXT"}
# the builder adds the N-terminal hydrogens itself
OMIT = {"H2"}

# --- small vector helpers (avoids a numpy dependency) ---


def sub(p, q):
    return [p[0] - q[0], p[1] - q[1], p[2] - q[2]]


def add(p, q):
    return [p[0] + q[0], p[1] + q[1], p[2] + q[2]]


def scale(p, s):
    return [p[0] * s, p[1] * s, p[2] * s]


def dot(p, q):
    return p[0] * q[0] + p[1] * q[1] + p[2] * q[2]


def cross(p, q):
    return [
        p[1] * q[2] - p[2] * q[1],
        p[2] * q[0] - p[0] * q[2],
        p[0] * q[1] - p[1] * q[0],
    ]


def norm(p):
    return math.sqrt(dot(p, p))


def unit(p):
    return scale(p, 1.0 / norm(p))


def distance(p, q):
    return norm(sub(p, q))


def angle(p, q, r):
    """Angle p-q-r in degrees."""
    u, v = unit(sub(p, q)), unit(sub(r, q))
    return math.degrees(math.acos(max(-1.0, min(1.0, dot(u, v)))))


def dihedral(p0, p1, p2, p3):
    """IUPAC signed torsion p0-p1-p2-p3 in degrees."""
    b0 = sub(p0, p1)
    b1 = unit(sub(p2, p1))
    b2 = sub(p3, p2)
    v = sub(b0, scale(b1, dot(b0, b1)))
    w = sub(b2, scale(b1, dot(b2, b1)))
    return math.degrees(math.atan2(dot(cross(b1, v), w), dot(v, w)))


def signed_volume(center, p, q, r):
    return dot(sub(p, center), cross(sub(q, center), sub(r, center)))


# --- CCD mmCIF ---


def read_loops(text):
    """Return {category: [row dict]} for every loop_ in an mmCIF file."""
    lines = text.splitlines()
    loops = {}
    i = 0
    while i < len(lines):
        if lines[i].strip() != "loop_":
            i += 1
            continue
        i += 1
        keys = []
        while i < len(lines) and lines[i].startswith("_"):
            keys.append(lines[i].strip())
            i += 1
        category = keys[0].split(".")[0]
        fields = [k.split(".")[1] for k in keys]
        rows, tokens = [], []
        while i < len(lines) and not lines[i].startswith(("_", "loop_", "#")):
            tokens += shlex.split(lines[i])
            i += 1
            while len(tokens) >= len(fields):
                rows.append(dict(zip(fields, tokens)))
                tokens = tokens[len(fields) :]
        loops[category] = rows
    return loops


def element_symbol(symbol):
    return symbol[0].upper() + symbol[1:].lower()


def read_component(code, ccd_dir):
    """Read (or download and cache) one CCD component."""
    path = os.path.join(ccd_dir, code + ".cif") if ccd_dir else None
    if path and os.path.exists(path):
        with open(path) as handle:
            text = handle.read()
    else:
        with urllib.request.urlopen(CCD_URL.format(code)) as response:
            text = response.read().decode()
        if path:
            os.makedirs(ccd_dir, exist_ok=True)
            with open(path, "w") as handle:
                handle.write(text)

    loops = read_loops(text)
    atoms = {}
    order = []
    for row in loops["_chem_comp_atom"]:
        xyz = [row["pdbx_model_Cartn_%s_ideal" % k] for k in "xyz"]
        if "?" in xyz:  # no ideal coordinates, so use the model ones
            xyz = [row["model_Cartn_%s" % k] for k in "xyz"]
        name = row["atom_id"]
        atoms[name] = {
            "element": element_symbol(row["type_symbol"]),
            "charge": int(row["charge"]) if row["charge"] not in "?." else 0,
            "leaving": row["pdbx_leaving_atom_flag"] == "Y",
            "xyz": [float(v) for v in xyz],
        }
        order.append(name)

    neighbors = {name: [] for name in order}
    for row in loops["_chem_comp_bond"]:
        a, b = row["atom_id_1"], row["atom_id_2"]
        neighbors[a].append(b)
        neighbors[b].append(a)
    orders = {}
    for row in loops["_chem_comp_bond"]:
        bond = {"SING": 1, "DOUB": 2, "TRIP": 3}[row["value_order"]]
        orders[(row["atom_id_1"], row["atom_id_2"])] = bond
        orders[(row["atom_id_2"], row["atom_id_1"])] = bond
    for name in order:
        atoms[name]["valence"] = sum(orders[(name, n)] for n in neighbors[name])
    return atoms, order, neighbors


# --- z-matrix construction ---


def layout(atoms, order, neighbors):
    """Choose the row order and the parent (RefA) of every atom."""
    kept = [n for n in order if n not in OMIT]
    for name in BACKBONE:
        if name not in kept:
            raise ValueError("no %s: not an amino acid backbone" % name)

    tail = [n for n in ("OXT", "HXT") if n in kept]
    tail += [n for n in kept if atoms[n]["leaving"] and n not in tail]

    rows = list(BACKBONE)
    parent = {"CA": "N", "C": "CA", "O": "C", "OXT": "C"}
    if "H" in kept and "H" not in tail:
        rows.append("H")
        parent["H"] = "N"

    # heavy side-chain atoms, breadth first out from CA
    queue = ["CA"]
    while queue:
        current = queue.pop(0)
        for other in neighbors[current]:
            if (
                other in rows
                or other in tail
                or other in OMIT
                or atoms[other]["element"] == "H"
            ):
                continue
            rows.append(other)
            parent[other] = current
            queue.append(other)

    # then the hydrogens, in the order of the atoms they sit on
    for heavy in list(rows):
        for other in neighbors[heavy]:
            if other in rows or other in tail or other in OMIT:
                continue
            rows.append(other)
            parent[other] = heavy

    for name in tail:
        if name not in parent:
            parent[name] = neighbors[name][0]
        rows.append(name)

    missing = set(kept) - set(rows)
    if missing:
        raise ValueError("unplaced atoms %s" % sorted(missing))
    return rows, parent


def references(name, rows, parent, xyz):
    """Pick (RefA, RefB, RefC) for one atom, following the layout rules."""
    if name in ("O", "OXT"):
        return "C", "CA", "N"
    if name == "H":
        return "N", "CA", "C"
    if name == "HXT":
        return "OXT", "C", "CA"
    a = parent[name]
    if a == "CA":
        return "CA", "N", "C"

    placed = rows[: rows.index(name)]
    usable = [n for n in placed if n not in NOT_A_REFERENCE and n != a]

    # walk back up the tree first (a standard chain z-matrix), then fall
    # back to any usable atom, skipping near-linear choices
    b_choices = [parent.get(a)] + usable
    for b in b_choices:
        if b is None or b not in usable:
            continue
        if not 5.0 < angle(xyz[name], xyz[a], xyz[b]) < 175.0:
            continue
        c_choices = [parent.get(b)] + usable
        for c in c_choices:
            if c is None or c not in usable or c == b:
                continue
            if 5.0 < angle(xyz[a], xyz[b], xyz[c]) < 175.0:
                return a, b, c
    raise ValueError("no usable references for %s" % name)


def build_zmatrix(code, ccd_dir):
    atoms, order, neighbors = read_component(code, ccd_dir)
    rows, parent = layout(atoms, order, neighbors)
    xyz = {n: atoms[n]["xyz"] for n in rows}
    index = {n: i + 1 for i, n in enumerate(rows)}

    zmat = []
    for i, name in enumerate(rows):
        a = b = c = None
        length = bend = torsion = 0.0
        if i == 1:
            a = "N"
        elif i == 2:
            a, b = "CA", "N"
        elif i > 2:
            a, b, c = references(name, rows, parent, xyz)
        if a:
            length = distance(xyz[name], xyz[a])
        if b:
            bend = angle(xyz[name], xyz[a], xyz[b])
        if c:
            torsion = dihedral(xyz[name], xyz[a], xyz[b], xyz[c])
        zmat.append(
            {
                "element": atoms[name]["element"],
                "name": name,
                "length": length,
                "angle": bend,
                "dihedral": torsion,
                "refs": [index[r] if r else 0 for r in (a, b, c)],
            }
        )
    return zmat


def format_zmatrix(zmat):
    lines = []
    for row in zmat:
        lines.append(
            "%-2s %-4s %7.4f %9.4f %9.4f %3d %3d %3d"
            % (
                row["element"],
                row["name"],
                row["length"],
                row["angle"],
                row["dihedral"],
                *row["refs"],
            )
        )
    return "\n".join(lines) + "\n"


# --- reading and checking ---


def read_zmatrix(path):
    """Parse a .zmat file the same way the peptide builder does."""
    zmat = []
    with open(path) as handle:
        for line in handle:
            parts = line.split()
            if len(parts) < 8:
                continue
            zmat.append(
                {
                    "element": parts[0],
                    "name": parts[1],
                    "length": float(parts[2]),
                    "angle": float(parts[3]),
                    "dihedral": float(parts[4]),
                    "refs": [int(v) for v in parts[5:8]],
                }
            )
    return zmat


def to_cartesian(zmat):
    """Build coordinates, matching Avogadro's internalToCartesian."""
    xyz = []
    for i, row in enumerate(zmat):
        a, b, c = (r - 1 for r in row["refs"])
        length = row["length"]
        theta = math.radians(row["angle"])
        if i == 0:
            xyz.append([0.0, 0.0, 0.0])
        elif i == 1:
            xyz.append(add(xyz[a], [length, 0.0, 0.0]))
        elif i == 2:
            ab = unit(sub(xyz[b], xyz[a]))
            perpendicular = [-ab[1], ab[0], 0.0]
            xyz.append(
                add(
                    xyz[a],
                    add(
                        scale(ab, length * math.cos(theta)),
                        scale(perpendicular, length * math.sin(theta)),
                    ),
                )
            )
        else:
            phi = math.radians(row["dihedral"])
            ab = unit(sub(xyz[a], xyz[b]))
            n = unit(cross(sub(xyz[b], xyz[c]), ab))
            m = cross(n, ab)
            offset = add(
                scale(ab, -length * math.cos(theta)),
                add(
                    scale(m, length * math.sin(theta) * math.cos(phi)),
                    scale(n, length * math.sin(theta) * math.sin(phi)),
                ),
            )
            xyz.append(add(xyz[a], offset))
    return xyz


def check_zmatrix(code, zmat, ccd_dir, exact):
    """Return a list of problems with one z-matrix.

    With exact=True (freshly generated files) every interatomic distance
    must reproduce the CCD ideal geometry.  Otherwise (hand-made or older
    files, which use other conformations) bonded lengths, stereochemistry
    and cis/trans arrangement are compared instead.
    """
    problems = []
    atoms, order, neighbors = read_component(code, ccd_dir)
    names = [row["name"] for row in zmat]
    row_of = {n: i for i, n in enumerate(names)}

    if len(set(names)) != len(names):
        dupes = sorted({n for n in names if names.count(n) > 1})
        problems.append("duplicate atom names %s" % dupes)
    expected = [n for n in order if n not in OMIT]
    for name in expected:
        if name not in row_of:
            problems.append("missing atom %s" % name)
    for row in zmat:
        name = row["name"]
        if name not in atoms:
            problems.append("%s is not a CCD atom of %s" % (name, code))
            continue
        if element_symbol(row["element"]) != atoms[name]["element"]:
            problems.append("%s should be %s" % (name, atoms[name]["element"]))
    if problems:
        return problems

    # the builder has no charges to read: it makes any nitrogen with four
    # bonds +1, so that rule has to reproduce the CCD's charges
    for name in names:
        inferred = 0
        if atoms[name]["element"] == "N" and atoms[name]["valence"] == 4:
            inferred = 1
        if inferred != atoms[name]["charge"]:
            problems.append(
                "%s: the builder would assign charge %d, CCD has %d"
                % (name, inferred, atoms[name]["charge"])
            )

    # layout rules the peptide builder depends on
    if names[:4] != list(BACKBONE):
        problems.append("rows 1-4 must be N CA C O, not %s" % names[:4])
    fixed = {"O": ("C", "CA", "N"), "OXT": ("C", "CA", "N")}
    if "H" in row_of:
        fixed["H"] = ("N", "CA", "C")
    for name, refs in fixed.items():
        if name in row_of:
            got = tuple(names[r - 1] for r in zmat[row_of[name]]["refs"])
            if got != refs:
                problems.append("%s must be placed from %s" % (name, refs))
    leaving = {n for n in names if atoms[n]["leaving"]}
    first_leaving = min((row_of[n] for n in leaving), default=len(names))
    for i, row in enumerate(zmat):
        name = row["name"]
        refs = [names[r - 1] for r in row["refs"] if r]
        if any(row_of[r] >= i for r in refs):
            problems.append("%s references a later row" % name)
        if name in BACKBONE or name in fixed:
            continue
        if refs and refs[0] == "CA" and refs != ["CA", "N", "C"]:
            problems.append("%s (on CA) must be placed from CA N C" % name)
        if name == "HXT":
            continue
        bad = [r for r in refs if r in NOT_A_REFERENCE]
        if bad:
            problems.append("%s references %s" % (name, bad))
        if i > first_leaving and any(row_of[r] >= first_leaving for r in refs):
            problems.append("%s references a row the builder may skip" % name)

    built = dict(zip(names, to_cartesian(zmat)))
    ideal = {n: atoms[n]["xyz"] for n in names}
    bonded = {(a, b) for a in names for b in neighbors[a] if b in built}

    if exact:
        worst = max(
            abs(distance(built[p], built[q]) - distance(ideal[p], ideal[q]))
            for p in names
            for q in names
        )
        if worst > 0.002:
            problems.append("geometry deviates from CCD by %.4f A" % worst)
    else:
        for a, b in bonded:
            dz, dc = distance(built[a], built[b]), distance(ideal[a], ideal[b])
            if a < b and abs(dz - dc) > 0.08:
                problems.append("bond %s-%s %.3f, CCD %.3f" % (a, b, dz, dc))
        for i, p in enumerate(names):
            for q in names[:i]:
                d = distance(built[p], built[q])
                if d < 1.6 and (p, q) not in bonded:
                    problems.append("%s..%s only %.2f A apart" % (p, q, d))
        # cis/trans about bonds between atoms with at most three neighbors
        for b, c in bonded:
            if b > c or len(neighbors[b]) > 3 or len(neighbors[c]) > 3:
                continue
            if atoms[b]["element"] == "H" or atoms[c]["element"] == "H":
                continue
            for a in neighbors[b]:
                for d in neighbors[c]:
                    if a == c or d == b or a not in built or d not in built:
                        continue
                    tc = dihedral(ideal[a], ideal[b], ideal[c], ideal[d])
                    tz = dihedral(built[a], built[b], built[c], built[d])
                    planar = min(abs(tc), 180.0 - abs(tc)) < 25.0
                    # rotamers about an O-H bond are conformation, not naming
                    if "HXT" in (a, d) or not planar:
                        continue
                    if (abs(tz) < 90.0) != (abs(tc) < 90.0):
                        problems.append(
                            "%s-%s-%s-%s is %.0f, CCD %.0f" % (a, b, c, d, tz, tc)
                        )

    # stereochemistry (and prochiral naming) at four-coordinate atoms; this
    # also catches a mirror image, which the distance check cannot
    for center in names:
        around = sorted(n for n in neighbors[center] if n in built)
        if len(around) < 4:
            continue
        p, q, r = around[:3]
        vz = signed_volume(built[center], built[p], built[q], built[r])
        vc = signed_volume(ideal[center], ideal[p], ideal[q], ideal[r])
        if (vz > 0) != (vc > 0):
            problems.append("wrong handedness or naming at %s" % center)
    return problems


def main():
    root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    parser = argparse.ArgumentParser(
        description=__doc__.split("\n")[0],
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument(
        "residues",
        nargs="*",
        default=RESIDUES,
        help="three-letter CCD codes (default: all peptide builder residues)",
    )
    parser.add_argument(
        "--output-dir",
        default=os.path.join(root, "amino"),
        help="where the .zmat files live",
    )
    parser.add_argument(
        "--ccd-dir",
        help="cache directory for CCD .cif files (default: download each run)",
    )
    parser.add_argument(
        "--check",
        action="store_true",
        help="validate the existing files instead of writing new ones",
    )
    args = parser.parse_args()

    failed = False
    for code in (r.upper() for r in args.residues):
        path = os.path.join(args.output_dir, code + ".zmat")
        if args.check:
            zmat = read_zmatrix(path)
        else:
            zmat = build_zmatrix(code, args.ccd_dir)
            # check what will actually be read back, including rounding
            with open(path, "w") as handle:
                handle.write(format_zmatrix(zmat))
            zmat = read_zmatrix(path)
        problems = check_zmatrix(code, zmat, args.ccd_dir, not args.check)
        print("%s: %s" % (code, "ok" if not problems else "FAILED"))
        for problem in problems:
            print("    " + problem)
        failed = failed or bool(problems)
    return 1 if failed else 0


if __name__ == "__main__":
    sys.exit(main())
