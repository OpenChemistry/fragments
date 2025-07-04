#!/usr/bin/env python3
# SPDX-License-Identifier: BSD-3-Clause
"""Script to flatten CJSON for structure templates of Avogadro2."""

import argparse
import json
from pathlib import Path

# Functions flatten_dumps and minimal contain the pairwise comments
# of `# fmt: off` and `# fmt: on`.  This locally disables an automated
# reformat by Black or Ruff which otherwise would overwrite the functions'
# intentional patterns of indentation, and use of single/double quotes.


def recursive_search(path: Path):
    file_list = []
    file_list.extend([x for x in path.iterdir() if x.is_file()])
    for d in [x for x in path.iterdir() if x.is_dir()]:
        file_list.extend(recursive_search(d))
    return file_list


def flatten_arrays(data: dict) -> dict:
    """Turn any lists of simple items (not dicts or lists) into strings."""
    if isinstance(data, list):
        # Turn simple lists into flat strings
        if all(not isinstance(i, (dict, list)) for i in data):
            return json.dumps(data)
        # Recursively flatten any nested lists
        else:
            items = [flatten_arrays(i) for i in data]
            return items
    elif isinstance(data, dict):
        # Recursively flatten all entries
        new = {k: flatten_arrays(v) for k, v in data.items()}
        return new
    else:
        return data


def flatten_dumps(data: dict) -> str:
    """Do the same as json.dumps() but write simple lists on a single line."""
    # fmt: off
    flattened = flatten_arrays(data)
    # Lists are now strings, remove quotes to turn them back into lists
    output = json.dumps(flattened, indent=2).replace('"[', '[').replace(']"', ']')
    # Any strings within lists will have had their quotes escaped, so get rid of escapes
    output = output.replace(r'\"', '"')
    return output
    # fmt: on


def minimal(cjson: dict) -> dict:
    """Reduce a CJSON to core geometry data.

    This retains the atoms with their coordinates, the bonds, and
    the overall charge."""
    # fmt: off
    minimal_cjson = {
        "chemicalJson": cjson.get("chemicalJson", 1),
        "atoms": {
            "coords": {
                "3d": cjson["atoms"]["coords"]["3d"]
            },
            "elements": {
                "number": cjson["atoms"]["elements"]["number"]
            }
        },
        "bonds": {
            "connections": {
                "index": cjson["bonds"]["connections"]["index"]
            },
            "order": cjson["bonds"]["order"]
        }
    }
    # Formal charges are useful but may or may not be there
    if "formalCharges" in cjson["atoms"]:
        minimal_cjson["atoms"]["formalCharges"] = cjson["atoms"]["formalCharges"]
    # Keep total charge if present
    if "properties" in cjson:
        minimal_cjson["properties"] = {}
        for prop in ["totalCharge"]:
            if prop in cjson["properties"]:
                minimal_cjson["properties"][prop] = cjson["properties"][prop]

    return minimal_cjson
    # fmt: on


def round_coords(cjson: dict, places: int) -> dict:
    """Round off the atomic coordinates in a CJSON to a specified number of decimal places."""
    coords = cjson["atoms"]["coords"]["3d"]
    rounded = [round(c, places) for c in coords]
    cjson["atoms"]["coords"]["3d"] = rounded
    return cjson


def flatten_all(
    cjson_list: list[Path],
    minimize: bool,
    round_coords_places: int | None = None,
    validate: bool = False,
):
    """Flatten a list of CJSON files according to the parameters set."""
    # checks = {}

    # Read then write each cjson
    for file in cjson_list:
        with open(file, "r", encoding="utf-8") as f:
            cjson = json.load(f)
        # print(cjson)
        if minimize:
            cjson = minimal(cjson)
        # print(cjson)
        if round_coords_places:
            cjson = round_coords(cjson, round_coords_places)
        # print(cjson)
        flattened = flatten_dumps(cjson)
        # print(cjson)
        # print(flattened)
        with open(file, "w", encoding="utf-8") as f:
            f.write(flattened)


#        if validate:
#            # Test we get the same object back as we originally read
#            check = (cjson == json.loads(flattened))
#            checks[file] = check
#            if check is False:
#                print(f"{file} was not validated")

#    if validate:
#        print(checks)


def get_args():
    """Collect command-line arguments"""
    parser = argparse.ArgumentParser(
        description="""
Script to flatten CJSON files of structure templates for Avogadro2""",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument(
        "directory",
        help="Directory of files to process",
        type=Path,
    )

    parser.add_argument(
        "-m",
        "--minimize",
        help="Reduce CJSON to retain only atoms, bonds, and charges",
        action="store_true",
    )

    parser.add_argument(
        "-r",
        "--round_coords",
        metavar="",
        help="""
Define explicitly the number of decimals of atomic coordinates to round
to.  If you choose 0 (zero), the coordinates are not rounded.""",
        type=int,
        default=4,
    )

    return parser.parse_args()


if __name__ == "__main__":
    args = get_args()
    # print(args)

    # Get all CJSON files in dir
    file_list = recursive_search(args.directory)
    cjson_list = [f for f in file_list if f.suffix == ".cjson"]

    flatten_all(cjson_list, args.minimize, args.round_coords)
