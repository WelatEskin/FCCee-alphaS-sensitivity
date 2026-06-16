#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
import numpy as np
import uproot


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-i",
        "--input",
        default="data/event_display_alphaS_0118_multialg.root",
        help="Input ROOT file",
    )
    parser.add_argument(
        "--event-id",
        type=int,
        default=10570,
        help="Event ID to inspect",
    )

    args = parser.parse_args()

    branches = [
        "event_id",
        "pt",
        "E",
        "durham_jet_id",
        "antikt_jet_id",
        "kt_jet_id",
        "cambridge_jet_id",
    ]

    with uproot.open(args.input) as f:
        tree = f["events"]
        arrays = tree.arrays(branches, library="np")

    mask = arrays["event_id"] == args.event_id

    if not np.any(mask):
        raise RuntimeError(f"Event not found: {args.event_id}")

    print("----------------------------------------")
    print(f"Event ID: {args.event_id}")
    print(f"Particles in event: {np.sum(mask)}")
    print("----------------------------------------")

    for b in [
        "durham_jet_id",
        "antikt_jet_id",
        "kt_jet_id",
        "cambridge_jet_id",
    ]:
        vals = arrays[b][mask]
        unique, counts = np.unique(vals, return_counts=True)

        print(b)
        print("  unique values:", unique)
        print("  counts       :", counts)
        print("  n(-1)        :", np.sum(vals < 0))
        print("  n(>=0)       :", np.sum(vals >= 0))
        print("  N_jets       :", len(np.unique(vals[vals >= 0])))
        print("----------------------------------------")


if __name__ == "__main__":
    main()