import argparse
import json
from pathlib import Path


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
    flattened = flatten_arrays(data)
    # Lists are now strings, remove quotes to turn them back into lists
    output = json.dumps(flattened, indent=2).replace('"[', '[').replace(']"', ']')
    # Any strings within lists will have had their quotes escaped, so get rid of escapes
    output = output.replace(r'\"', '"')
    return output


def minimal(cjson: dict) -> dict:
    """Reduce a CJSON to core geometry data."""
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


def round_coords(cjson: dict, places: int) -> dict:
    """Round off the atomic coordinates in a CJSON to the specified number of decimal
    places."""
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
    checks = {}

    # Read then write each cjson
    for file in cjson_list:
        with open(file) as f:
            cjson = json.load(f)
        print(cjson)
        if minimize:
            cjson = minimal(cjson)
        print(cjson)
        if round_coords_places:
            cjson = round_coords(cjson, round_coords_places)
        print(cjson)
        flattened = flatten_dumps(cjson)
        print(cjson)
        #print(flattened)
        with open(file, "w") as f:
            f.write(flattened)
        if validate:
            # Test we get the same object back as we originally read
            check = (cjson == json.loads(flattened))
            checks[file] = check
            if check is False:
                print(f"{file} was not validated")
    
    if validate:
        print(checks)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("directory", type=Path)
    parser.add_argument("-m", "--minimize", action="store_true")
    parser.add_argument("-r", "--round_coords", nargs="?", type=int, const=5, default=None)
    args = parser.parse_args()
    print(args)

    # Get all CJSON files in dir
    file_list = recursive_search(args.directory)
    cjson_list = [f for f in file_list if f.suffix == ".cjson"]

    flatten_all(cjson_list, args.minimize, args.round_coords)
