#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
import uproot


def get_tree(root_file):
    candidate_names = ["events", "Events", "tree", "Tree", "t", "T"]

    for name in candidate_names:
        if name in root_file:
            return root_file[name]

    for key in root_file.keys():
        obj = root_file[key]
        if hasattr(obj, "arrays"):
            return obj

    raise RuntimeError("ROOT dosyasında okunabilir TTree bulunamadı.")


def main():
    parser = argparse.ArgumentParser(
        description="Check available branches in FCC-ee event-display ROOT file."
    )

    parser.add_argument(
        "-i",
        "--input",
        default="data/event_display_alphaS_0118.root",
        help="Input ROOT file",
    )

    args = parser.parse_args()

    with uproot.open(args.input) as f:
        tree = get_tree(f)
        branches = list(tree.keys())

    print("----------------------------------------")
    print("All branches")
    print("----------------------------------------")
    for b in branches:
        print(b)

    print("----------------------------------------")
    print("Possible jet-ID branches")
    print("----------------------------------------")

    jet_like = [
        b for b in branches
        if (
            "jet" in b.lower()
            or "durham" in b.lower()
            or "antikt" in b.lower()
            or "anti" in b.lower()
            or "cambridge" in b.lower()
            or b.lower().startswith("ca_")
        )
    ]

    if len(jet_like) == 0:
        print("No jet-like branch found.")
    else:
        for b in jet_like:
            print(b)

    print("----------------------------------------")


if __name__ == "__main__":
    main()