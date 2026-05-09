#!/usr/bin/env python3
"""Quote MAT Newick labels that contain Newick control characters.

Some public datasets use sample names with semicolons. Those names are valid
sample identifiers, but unquoted semicolons make the embedded MAT Newick invalid
for parsers such as treeswift/taxoniumtools.
"""

import argparse
import gzip
import shutil
from pathlib import Path

from taxoniumtools import parsimony_pb2


def quote_label(label: str) -> str:
    return "'" + label.replace("'", "''") + "'"


def quote_semicolon_labels(newick: str) -> str:
    """Quote node labels after '(' or ',' when the label contains ';'."""
    output = []
    i = 0
    length = len(newick)

    while i < length:
        char = newick[i]
        output.append(char)
        i += 1

        if char not in "(,":
            continue

        if i >= length or newick[i] in "(),:;":
            continue

        if newick[i] == "'":
            continue

        start = i
        while i < length and newick[i] not in ":,)":
            i += 1

        label = newick[start:i]
        output.append(quote_label(label) if ";" in label else label)

    return "".join(output)


def open_maybe_gzip(path: Path, mode: str):
    if path.suffix == ".gz":
        return gzip.open(path, mode)
    return open(path, mode)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()

    with open_maybe_gzip(args.input, "rb") as input_handle:
        tree_data = parsimony_pb2.data()
        tree_data.ParseFromString(input_handle.read())

    sanitized_newick = quote_semicolon_labels(tree_data.newick)
    if sanitized_newick == tree_data.newick:
        shutil.copyfile(args.input, args.output)
        return

    tree_data.newick = sanitized_newick
    with open_maybe_gzip(args.output, "wb") as output_handle:
        output_handle.write(tree_data.SerializeToString())


if __name__ == "__main__":
    main()
