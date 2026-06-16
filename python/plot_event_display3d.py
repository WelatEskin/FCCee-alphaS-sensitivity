#!/usr/bin/env python3
# plot_event_display_3d.py
#
# 3D event-display style plot for FCC-ee alphaS study.
#
# Input:
#   data/event_display_alphaS_0118.root
#
# Output:
#   figures/event_display_3d_durham_vs_antikt_event17.png
#   figures/event_display_3d_durham_vs_antikt_event17.pdf
#
# Usage:
#   python3 scripts/plot_event_display_3d.py

import os
import numpy as np
import matplotlib.pyplot as plt

from mpl_toolkits.mplot3d import Axes3D  # noqa: F401

try:
    import uproot
except ImportError:
    raise ImportError(
        "uproot is required. Install it with:\n"
        "  pip install uproot awkward matplotlib numpy\n"
    )


INPUT_FILE = "data/event_display_alphaS_0118.root"
OUTPUT_PNG = "figures/event_display_3d_durham_vs_antikt_event17.png"
OUTPUT_PDF = "figures/event_display_3d_durham_vs_antikt_event17.pdf"

TARGET_EVENT = 17

# Use rapidity for the horizontal coordinate.
# For an e+e- event display this is only an illustrative projection.
X_VARIABLE = "rapidity"   # alternatives: "eta"
Z_VARIABLE = "pt"         # alternatives: "E"

COLORS = [
    "tab:red",
    "tab:blue",
    "tab:green",
    "tab:purple",
    "tab:orange",
    "tab:cyan",
    "tab:pink",
    "tab:brown",
    "tab:olive",
    "tab:gray",
]


def jet_color(jet_id: int) -> str:
    if jet_id < 0:
        return "lightgray"
    return COLORS[jet_id % len(COLORS)]


def marker_size(pt):
    size = 12.0 + 4.0 * np.sqrt(max(pt, 0.0))
    return min(size, 80.0)


def load_event():
    with uproot.open(INPUT_FILE) as f:
        event_info = f["EventInfo"].arrays(library="np")
        particles = f["Particles"].arrays(library="np")

    event_mask = event_info["event_id"] == TARGET_EVENT

    if not np.any(event_mask):
        raise RuntimeError(f"Event {TARGET_EVENT} was not found in EventInfo.")

    info = {key: val[event_mask][0] for key, val in event_info.items()}

    particle_mask = particles["event_id"] == TARGET_EVENT
    event_particles = {key: val[particle_mask] for key, val in particles.items()}

    return info, event_particles


def setup_3d_axes(ax, title):
    ax.set_title(title, fontsize=17, pad=12, fontweight="bold")

    ax.set_xlabel(r"$y$", labelpad=9, fontsize=13)
    ax.set_ylabel(r"$\phi$", labelpad=9, fontsize=13)
    ax.set_zlabel(r"$p_T$ [GeV]", labelpad=9, fontsize=13)

    ax.set_xlim(-5.5, 5.5)
    ax.set_ylim(-3.2, 3.2)
    ax.set_zlim(0.0, 30.0)

    ax.view_init(elev=24, azim=-58)

    # Light grid, white background
    ax.xaxis.pane.set_facecolor((1.0, 1.0, 1.0, 1.0))
    ax.yaxis.pane.set_facecolor((1.0, 1.0, 1.0, 1.0))
    ax.zaxis.pane.set_facecolor((1.0, 1.0, 1.0, 1.0))

    ax.grid(True, linewidth=0.4, alpha=0.35)

    ax.tick_params(axis="both", which="major", labelsize=10)


def draw_event(ax, particles, algorithm="durham"):
    if algorithm == "durham":
        jet_branch = "durham_jet_id"
    elif algorithm == "antikt":
        jet_branch = "antikt_jet_id"
    else:
        raise ValueError("algorithm must be 'durham' or 'antikt'")

    x = particles[X_VARIABLE]
    phi = particles["phi"]
    z = particles[Z_VARIABLE]
    pt = particles["pt"]
    jet_ids = particles[jet_branch]

    # Avoid pathological rapidity/eta values
    finite_mask = np.isfinite(x) & np.isfinite(phi) & np.isfinite(z)
    range_mask = (x > -5.5) & (x < 5.5) & (phi > -3.2) & (phi < 3.2)
    mask = finite_mask & range_mask

    x = x[mask]
    phi = phi[mask]
    z = z[mask]
    pt = pt[mask]
    jet_ids = jet_ids[mask]

    # Draw vertical bars and top markers
    for xi, phii, zi, pti, jid in zip(x, phi, z, pt, jet_ids):
        c = jet_color(int(jid))

        ax.plot(
            [xi, xi],
            [phii, phii],
            [0.0, zi],
            color=c,
            linewidth=1.2,
            alpha=0.90,
        )

        ax.scatter(
            xi,
            phii,
            zi,
            color=c,
            s=marker_size(pti),
            edgecolor="black",
            linewidth=0.25,
            alpha=0.95,
            depthshade=True,
        )

    return len(x)


def add_text(ax, lines, x2d=0.04, y2d=0.92):
    text = "\n".join(lines)
    ax.text2D(
        x2d,
        y2d,
        text,
        transform=ax.transAxes,
        fontsize=11,
        va="top",
        ha="left",
        bbox=dict(
            boxstyle="round,pad=0.35",
            facecolor="white",
            edgecolor="none",
            alpha=0.80,
        ),
    )


def main():
    os.makedirs("figures", exist_ok=True)

    info, particles = load_event()

    n_durham = int(info["nJets_durham"])
    n_antikt = int(info["nJets_antikt"])
    e_durham = float(info["leadingJetEnergy_durham"])
    e_antikt = float(info["leadingJetEnergy_antikt"])
    y23 = float(info["y23"])

    fig = plt.figure(figsize=(15.5, 7.2))

    ax1 = fig.add_subplot(1, 2, 1, projection="3d")
    ax2 = fig.add_subplot(1, 2, 2, projection="3d")

    setup_3d_axes(ax1, r"Durham $e^+e^-\, k_T$")
    setup_3d_axes(ax2, r"anti-$k_T$-like")

    n_drawn_durham = draw_event(ax1, particles, algorithm="durham")
    n_drawn_antikt = draw_event(ax2, particles, algorithm="antikt")

    add_text(
        ax1,
        [
            r"FCC-ee Simulation",
            rf"Event {TARGET_EVENT}, $\alpha_s=0.118$",
            rf"$N_{{jet}}={n_durham}$, $E_{{lead}}={e_durham:.1f}$ GeV",
            rf"$y_{{23}}={y23:.4f}$, $y_{{cut}}=0.002$",
            rf"Visible particles: {n_drawn_durham}",
        ],
        x2d=0.05,
        y2d=0.93,
    )

    add_text(
        ax2,
        [
            r"FCC-ee Simulation",
            rf"Event {TARGET_EVENT}, $\alpha_s=0.118$",
            rf"$N_{{jet}}={n_antikt}$, $E_{{lead}}={e_antikt:.1f}$ GeV",
            r"$e^+e^-$ anti-$k_T$-like, $E_{jet}>1$ GeV",
            rf"Visible particles: {n_drawn_antikt}",
        ],
        x2d=0.05,
        y2d=0.93,
    )

    fig.suptitle(
        r"Same-event jet clustering comparison at generator level",
        fontsize=18,
        fontweight="bold",
        y=0.98,
    )

    plt.tight_layout(rect=[0.0, 0.0, 1.0, 0.95])

    fig.savefig(OUTPUT_PNG, dpi=300)
    fig.savefig(OUTPUT_PDF)

    print("Selected event:", TARGET_EVENT)
    print("Durham Njet       =", n_durham)
    print("anti-kT-like Njet =", n_antikt)
    print("Durham Elead      =", f"{e_durham:.4f} GeV")
    print("anti-kT Elead     =", f"{e_antikt:.4f} GeV")
    print("y23               =", f"{y23:.6f}")
    print("Output written:")
    print(" ", OUTPUT_PNG)
    print(" ", OUTPUT_PDF)


if __name__ == "__main__":
    main()