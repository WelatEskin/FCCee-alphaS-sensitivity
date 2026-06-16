#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import argparse
import numpy as np
import matplotlib.pyplot as plt
import uproot


# ============================================================
# Physics / label settings
# ============================================================

SQRT_S_GEV = 91.2
ALPHA_S = 0.118


# ============================================================
# Plot style settings
# ============================================================

FIGSIZE = (13.0, 9.0)
FIG_DPI = 400

Y_LIM = (-2.1, 2.4)
PHI_LIM = (-3.1, 3.1)
PT_LIM = (0.0, 10.0)

BAR_DY = 0.085
BAR_DPHI = 0.085
BAR_ALPHA = 0.98
BAR_LINEWIDTH = 0.30

MAX_PARTICLES_SHOWN = 180
MIN_PT_DISPLAY = 0.05

TITLE_SIZE = 25
LABEL_SIZE = 18
TICK_SIZE = 13
INFOBOX_SIZE = 13.5


# ============================================================
# Algorithms to plot
# ============================================================

ALGORITHMS = [
    {
        "key": "durham",
        "branch": "durham_jet_id",
        "title": r"Durham $e^+e^-$ $k_T$",
        "output": "event_display_durham_prism",
    },
    {
        "key": "antikt",
        "branch": "antikt_jet_id",
        "title": r"anti-$k_T$-like",
        "output": "event_display_antikt_prism",
    },
    {
        "key": "kt",
        "branch": "kt_jet_id",
        "title": r"$k_T$-like",
        "output": "event_display_kt_prism",
    },
    {
        "key": "cambridge",
        "branch": "cambridge_jet_id",
        "title": r"Cambridge/Aachen-like",
        "output": "event_display_cambridge_prism",
    },
]


# ============================================================
# Utilities
# ============================================================

def normalize_phi(phi):
    return (phi + np.pi) % (2.0 * np.pi) - np.pi


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


def palette():
    return [
        "#e41a1c",  # red
        "#0047ff",  # blue
        "#24b33b",  # green
        "#f0a202",  # orange
        "#00cfd1",  # cyan
        "#d800d8",  # magenta
        "#8c6d00",  # olive / brown
        "#9e9e9e",  # grey
        "#ffd400",  # yellow
        "#7b2cbf",  # purple
        "#a52a2a",  # brown
        "#f4b6c2",  # pink
        "#00a087",  # teal
        "#ff7f00",  # bright orange
        "#4d4d4d",  # dark grey
        "#80b1d3",  # soft blue
    ]


def colors_from_labels(labels):
    cols = palette()
    output = []

    for lab in labels:
        lab = int(lab)

        if lab < 0:
            output.append("#cfcfcf")
        else:
            output.append(cols[lab % len(cols)])

    return output


def count_jets(labels):
    labels = np.asarray(labels, dtype=int)
    labels = labels[labels >= 0]

    if len(labels) == 0:
        return 0

    return len(np.unique(labels))


def leading_jet_energy(E, labels):
    E = np.asarray(E, dtype=float)
    labels = np.asarray(labels, dtype=int)

    valid = labels >= 0

    if not np.any(valid):
        return 0.0

    max_energy = 0.0

    for lab in np.unique(labels[valid]):
        energy_sum = np.sum(E[labels == lab])

        if energy_sum > max_energy:
            max_energy = energy_sum

    return float(max_energy)


# ============================================================
# 3D style
# ============================================================

def apply_3d_style(ax):
    ax.set_facecolor("white")

    # Side panes
    ax.xaxis.pane.set_facecolor((0.98, 0.98, 0.98, 1.0))
    ax.yaxis.pane.set_facecolor((0.98, 0.98, 0.98, 1.0))

    # Floor / base plane
    ax.zaxis.pane.set_facecolor((0.84, 0.84, 0.84, 1.0))

    ax.xaxis.pane.set_edgecolor((0.65, 0.65, 0.65, 1.0))
    ax.yaxis.pane.set_edgecolor((0.65, 0.65, 0.65, 1.0))
    ax.zaxis.pane.set_edgecolor((0.55, 0.55, 0.55, 1.0))

    ax.grid(True)

    for axis in [ax.xaxis, ax.yaxis, ax.zaxis]:
        axis._axinfo["grid"]["color"] = (0.68, 0.68, 0.68, 1.0)
        axis._axinfo["grid"]["linewidth"] = 0.85
        axis._axinfo["grid"]["linestyle"] = "-"
        axis._axinfo["axisline"]["linewidth"] = 1.15
        axis._axinfo["axisline"]["color"] = (0.25, 0.25, 0.25, 1.0)

    ax.tick_params(axis="both", which="major", labelsize=TICK_SIZE, width=1.2)
    ax.tick_params(axis="z", which="major", labelsize=TICK_SIZE, width=1.2)


# ============================================================
# ROOT reading
# ============================================================

def read_flat_particle_tree(input_file):
    with uproot.open(input_file) as f:
        tree = get_tree(f)
        branches = list(tree.keys())

        print("----------------------------------------")
        print("Available branches:")
        for b in branches:
            print(" ", b)
        print("----------------------------------------")

        required = [
            "event_id",
            "pt",
            "phi",
            "E",
        ]

        for alg in ALGORITHMS:
            required.append(alg["branch"])

        for r in required:
            if r not in branches:
                raise KeyError(f"Gerekli branch eksik: {r}")

        if "rapidity" in branches:
            y_branch = "rapidity"
        elif "eta" in branches:
            y_branch = "eta"
        else:
            raise KeyError("Ne 'rapidity' ne de 'eta' branch'i bulundu.")

        read_branches = [
            "event_id",
            "pt",
            "phi",
            "E",
            y_branch,
        ] + [alg["branch"] for alg in ALGORITHMS]

        arrays = tree.arrays(read_branches, library="np")

    data = {
        "event_id": np.asarray(arrays["event_id"], dtype=int),
        "pt": np.asarray(arrays["pt"], dtype=float),
        "phi": normalize_phi(np.asarray(arrays["phi"], dtype=float)),
        "E": np.asarray(arrays["E"], dtype=float),
        "y": np.asarray(arrays[y_branch], dtype=float),
    }

    for alg in ALGORITHMS:
        data[alg["key"]] = np.asarray(arrays[alg["branch"]], dtype=int)

    return data


def build_event(data, event_id=None):
    all_ids = np.unique(data["event_id"])

    if event_id is None:
        best_id = None
        best_score = -1.0e99

        for eid in all_ids:
            mask = data["event_id"] == eid

            pt = data["pt"][mask]

            n_particles = np.sum(pt > MIN_PT_DISPLAY)
            lead_pt = np.max(pt) if len(pt) > 0 else 0.0

            n_durham = count_jets(data["durham"][mask])
            n_antikt = count_jets(data["antikt"][mask])
            n_kt = count_jets(data["kt"][mask])
            n_cambridge = count_jets(data["cambridge"][mask])

            # Prefer a visually rich but not pathological event.
            score = (
                1.00 * min(n_particles, 160)
                + 7.50 * min(n_durham, 8)
                + 4.00 * min(n_antikt, 18)
                + 4.00 * min(n_kt, 18)
                + 4.00 * min(n_cambridge, 18)
                + 3.00 * lead_pt
            )

            if 4 <= n_durham <= 8:
                score += 80.0

            if n_particles < 50:
                score -= 1000.0

            if score > best_score:
                best_score = score
                best_id = eid

        event_id = int(best_id)

    if event_id not in all_ids:
        raise ValueError(f"Event bulunamadı: {event_id}")

    mask = data["event_id"] == event_id

    ev = {
        "event_id": int(event_id),
        "pt": data["pt"][mask],
        "phi": data["phi"][mask],
        "E": data["E"][mask],
        "y": data["y"][mask],
    }

    for alg in ALGORITHMS:
        ev[alg["key"]] = data[alg["key"]][mask]

    finite = (
        np.isfinite(ev["pt"])
        & np.isfinite(ev["phi"])
        & np.isfinite(ev["y"])
        & np.isfinite(ev["E"])
        & (ev["pt"] > MIN_PT_DISPLAY)
    )

    base_keys = ["pt", "phi", "E", "y"]
    alg_keys = [alg["key"] for alg in ALGORITHMS]

    for key in base_keys + alg_keys:
        ev[key] = ev[key][finite]

    return ev


def downselect(ev, labels):
    pt = ev["pt"]
    phi = ev["phi"]
    y = ev["y"]
    E = ev["E"]
    labels = np.asarray(labels, dtype=int)

    if len(pt) <= MAX_PARTICLES_SHOWN:
        return y, phi, pt, E, labels

    order = np.argsort(pt)[::-1]
    keep = np.sort(order[:MAX_PARTICLES_SHOWN])

    return y[keep], phi[keep], pt[keep], E[keep], labels[keep]


# ============================================================
# Plotting
# ============================================================

def add_info_box(ax, ev, njet, elead):
    n_visible = len(ev["pt"])

    text = (
        r"$\bf{FCC{-}ee\ Simulation}$" + "\n"
        + rf"$\sqrt{{s}} = {SQRT_S_GEV:.1f}\ \mathrm{{GeV}},$ "
        + rf"$\bf{{Event\ {ev['event_id']}}},$ "
        + rf"$\alpha_s = {ALPHA_S:.3f}$" + "\n"
        + rf"$N_{{jet}} = {njet},\ E_{{lead}} = {elead:.1f}\ \mathrm{{GeV}}$" + "\n"
        + rf"$\bf{{Visible\ particles\ shown:\ {n_visible}}}$"
    )

    ax.text2D(
        0.060,
        0.840,
        text,
        transform=ax.transAxes,
        fontsize=INFOBOX_SIZE,
        ha="left",
        va="top",
        linespacing=1.03,
        bbox=dict(
            boxstyle="round,pad=0.35",
            facecolor="white",
            edgecolor="black",
            linewidth=1.0,
            alpha=0.90,
        ),
    )


def add_manual_zlabel(ax):
    ax.set_zlabel("")

    ax.text(
        Y_LIM[1] + 0.18,
        PHI_LIM[1] + 0.05,
        PT_LIM[1] + 0.25,
        r"$p_T\ [\mathrm{GeV}]$",
        fontsize=LABEL_SIZE,
        fontweight="bold",
        ha="center",
        va="bottom",
        clip_on=False,
    )


def plot_one_event(ev, algorithm_key, title, output_base):
    labels_full = ev[algorithm_key]
    y, phi, pt, E, labels = downselect(ev, labels_full)

    njet = count_jets(labels)
    elead = leading_jet_energy(E, labels)
    colors = colors_from_labels(labels)

    x_pos = y - BAR_DY / 2.0
    y_pos = phi - BAR_DPHI / 2.0
    z_pos = np.zeros_like(pt)

    dx = np.full_like(pt, BAR_DY)
    dy = np.full_like(pt, BAR_DPHI)
    dz = np.clip(pt, PT_LIM[0], PT_LIM[1])

    fig = plt.figure(figsize=FIGSIZE)
    ax = fig.add_subplot(111, projection="3d")

    apply_3d_style(ax)

    ax.bar3d(
        x_pos,
        y_pos,
        z_pos,
        dx,
        dy,
        dz,
        color=colors,
        edgecolor="black",
        linewidth=BAR_LINEWIDTH,
        alpha=BAR_ALPHA,
        shade=True,
        zsort="average",
    )

    ax.set_xlim(*Y_LIM)
    ax.set_ylim(*PHI_LIM)
    ax.set_zlim(*PT_LIM)

    ax.set_xlabel(
        r"$y$",
        fontsize=LABEL_SIZE,
        fontweight="bold",
        labelpad=10,
    )

    ax.set_ylabel(
        r"$\phi$",
        fontsize=LABEL_SIZE,
        fontweight="bold",
        labelpad=12,
    )

    add_manual_zlabel(ax)

    ax.set_title(
        title,
        fontsize=TITLE_SIZE,
        fontweight="bold",
        pad=16,
    )

    ax.set_xticks([-2, -1, 0, 1, 2])
    ax.set_yticks([-3, -2, -1, 0, 1, 2, 3])
    ax.set_zticks([0, 2, 4, 6, 8, 10])

    ax.view_init(elev=24, azim=-58)

    try:
        ax.set_box_aspect((1.45, 1.25, 0.80))
    except Exception:
        pass

    add_info_box(ax, ev, njet, elead)

    fig.subplots_adjust(
        left=0.02,
        right=0.98,
        bottom=0.02,
        top=0.92,
    )

    png_path = output_base + ".png"
    pdf_path = output_base + ".pdf"

    fig.savefig(
        png_path,
        dpi=FIG_DPI,
        bbox_inches="tight",
        pad_inches=0.35,
        facecolor="white",
    )

    fig.savefig(
        pdf_path,
        bbox_inches="tight",
        pad_inches=0.35,
        facecolor="white",
    )

    plt.close(fig)

    return png_path, pdf_path, njet, elead


# ============================================================
# Main
# ============================================================

def main():
    parser = argparse.ArgumentParser(
        description="Multi-algorithm FCC-ee prism-style event display."
    )

    parser.add_argument(
        "-i",
        "--input",
        default="data/event_display_alphaS_0118_multialg.root",
        help="Input ROOT file",
    )

    parser.add_argument(
        "-o",
        "--output-dir",
        default="figures_multialg",
        help="Output directory",
    )

    parser.add_argument(
        "--event-id",
        type=int,
        default=None,
        help="Specific event_id to plot",
    )

    args = parser.parse_args()

    input_file = args.input
    output_dir = args.output_dir

    print("----------------------------------------")
    print("FCC-ee multi-algorithm prism event display")
    print("----------------------------------------")
    print(f"Input file : {input_file}")
    print(f"Output dir : {output_dir}")
    print("----------------------------------------")

    if not os.path.exists(input_file):
        raise FileNotFoundError(
            f"Input ROOT dosyası bulunamadı:\n{input_file}\n\n"
            "Komutu proje ana klasöründen çalıştır veya -i ile tam yolu ver."
        )

    os.makedirs(output_dir, exist_ok=True)

    data = read_flat_particle_tree(input_file)
    ev = build_event(data, event_id=args.event_id)

    print("Selected event summary")
    print("----------------------------------------")
    print(f"event_id               = {ev['event_id']}")
    print(f"nParticles             = {len(ev['pt'])}")

    for alg in ALGORITHMS:
        njet = count_jets(ev[alg["key"]])
        elead = leading_jet_energy(ev["E"], ev[alg["key"]])
        print(f"nJets_{alg['key']:<12} = {njet}")
        print(f"Elead_{alg['key']:<12} [GeV] = {elead:.4f}")

    print("----------------------------------------")

    written = []

    for alg in ALGORITHMS:
        output_base = os.path.join(output_dir, alg["output"])

        png, pdf, njet, elead = plot_one_event(
            ev=ev,
            algorithm_key=alg["key"],
            title=alg["title"],
            output_base=output_base,
        )

        written.append((png, pdf))

    print("Output written:")
    for png, pdf in written:
        print(f"  {png}")
        print(f"  {pdf}")

    print("----------------------------------------")


if __name__ == "__main__":
    main()