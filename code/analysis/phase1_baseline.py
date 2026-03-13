#!/usr/bin/env python3
"""
ZMET-NID Phase 1 baseline analysis
====================================

Produces the first analytical artifacts from the selected Z→μμ sample:
  1. Inclusive dimuon mass + MET (diagnostic + baseline)
  2. MET stratified by pT(Z)           — 4 bins from protocol §4
  3. MET stratified by jet multiplicity — 0 / 1 / ≥2 jets (pT > 30 GeV)

Outputs
-------
  figures/inclusive_fit/met_inclusive.png
  figures/stratified_tests/met_by_zpt.png
  figures/stratified_tests/met_by_njets.png
  data/phase0_outputs/baseline_counts.json

Usage
-----
  python code/analysis/phase1_baseline.py --file data/raw/test.root
  python code/analysis/phase1_baseline.py --file data/raw/test.root --max-events 200000
"""

import sys
import json
import argparse
from pathlib import Path

import uproot
import awkward as ak
import numpy as np
import vector
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

vector.register_awkward()

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from utils.event_selection import select_good_muons, PT_LEAD_MIN, ZMASS_LOW, ZMASS_HIGH

FIGURES_INCL   = Path("figures/inclusive_fit")
FIGURES_STRAT  = Path("figures/stratified_tests")
OUTPUTS_DIR    = Path("data/phase1_outputs")

# Protocol §4 pT(Z) bin edges [GeV]
ZPT_EDGES = [0, 10, 30, 60, float("inf")]
ZPT_LABELS = ["pT(Z) < 10", "10 ≤ pT(Z) < 30", "30 ≤ pT(Z) < 60", "pT(Z) ≥ 60"]

# MET histogram range and bins
MET_RANGE = (0, 200)
MET_BINS  = 40


# ---------------------------------------------------------------------------
# Event selection returning the full Z→μμ event record
# ---------------------------------------------------------------------------

def build_zmumu_sample(events,
                       id_branch="Muon_tightId",
                       iso_branch="Muon_pfRelIso04_all"):
    """
    Apply Z→μμ selection (protocol §2) and return a flat record with one row
    per selected event containing:
      dimuon_mass  (GeV)
      dimuon_pt    (GeV)
      met_pt       (GeV)
      jet_mult     (int, jets with pT > 30 GeV)

    Returns
    -------
    sample : dict of numpy arrays, same length
    n_before : total events before selection
    """
    n_before = len(events["nMuon"])

    # -- per-muon quality cuts
    good_muons = select_good_muons(events, id_branch, iso_branch)

    # -- build all opposite-sign pairs per event
    mu1, mu2 = ak.unzip(ak.combinations(good_muons, 2))
    os = mu1.charge != mu2.charge
    mu1, mu2 = mu1[os], mu2[os]

    pair = mu1 + mu2          # Momentum4D four-vector addition
    masses  = pair.mass
    pts_dz  = pair.pt
    lead_pts = ak.where(mu1.pt >= mu2.pt, mu1.pt, mu2.pt)

    # pick pair closest to Z pole per event
    best_idx  = ak.argmin(np.abs(masses - 91.1876), axis=1, keepdims=True)
    best_mass = ak.firsts(masses[best_idx])
    best_pt   = ak.firsts(pts_dz[best_idx])
    best_lead = ak.firsts(lead_pts[best_idx])

    # retain selected muon eta/phi for jet cleaning (fill None → -999 so dR≫0.4)
    mu1_eta = ak.fill_none(ak.firsts(mu1.eta[best_idx]), -999.0)
    mu1_phi = ak.fill_none(ak.firsts(mu1.phi[best_idx]), -999.0)
    mu2_eta = ak.fill_none(ak.firsts(mu2.eta[best_idx]), -999.0)
    mu2_phi = ak.fill_none(ak.firsts(mu2.phi[best_idx]), -999.0)

    # event-level mask: has a valid OS pair with lead pT ≥ 25 GeV
    has_pair = (~ak.is_none(best_mass)) & (ak.fill_none(best_lead, 0.0) >= PT_LEAD_MIN)

    # Z mass window
    m_filled = ak.fill_none(best_mass, -1.0)
    in_window = (m_filled >= ZMASS_LOW) & (m_filled <= ZMASS_HIGH)

    sel = has_pair & in_window

    # -- jet multiplicity: pT > 30 GeV, cleaned against Z muons (ΔR > 0.4)
    # Jets overlapping the Z muons are NOT counted — they are the muons.
    j_eta = events["Jet_eta"]
    j_phi = events["Jet_phi"]
    j_pt  = events["Jet_pt"]

    def _dphi(a, b):
        d = a - b
        return ak.where(d > np.pi, d - 2*np.pi, ak.where(d < -np.pi, d + 2*np.pi, d))

    # broadcast per-event muon scalars to match ragged jet axis
    j_eta_bc, mu1_eta_bc = ak.broadcast_arrays(j_eta, mu1_eta)
    j_phi_bc, mu1_phi_bc = ak.broadcast_arrays(j_phi, mu1_phi)
    dR1 = np.sqrt((j_eta_bc - mu1_eta_bc)**2 + _dphi(j_phi_bc, mu1_phi_bc)**2)

    j_eta_bc, mu2_eta_bc = ak.broadcast_arrays(j_eta, mu2_eta)
    j_phi_bc, mu2_phi_bc = ak.broadcast_arrays(j_phi, mu2_phi)
    dR2 = np.sqrt((j_eta_bc - mu2_eta_bc)**2 + _dphi(j_phi_bc, mu2_phi_bc)**2)

    clean_jet = (dR1 > 0.4) & (dR2 > 0.4)
    jet_count = ak.sum((j_pt > 30.0) & clean_jet, axis=1)

    # -- assemble flat arrays for selected events
    dimuon_mass = ak.to_numpy(ak.fill_none(best_mass[sel], 0.0))
    dimuon_pt   = ak.to_numpy(ak.fill_none(best_pt[sel],   0.0))
    met_pt      = ak.to_numpy(events["MET_pt"][sel])
    jet_mult    = ak.to_numpy(jet_count[sel])

    sample = {
        "dimuon_mass": dimuon_mass,
        "dimuon_pt":   dimuon_pt,
        "met_pt":      met_pt,
        "jet_mult":    jet_mult,
    }

    # Optional MET vector components for Family C (included if branches present)
    for branch, key in [
        ("MET_phi",                  "met_phi"),
        ("MET_MetUnclustEnUpDeltaX", "met_delta_x"),
        ("MET_MetUnclustEnUpDeltaY", "met_delta_y"),
    ]:
        if branch in ak.fields(events):
            sample[key] = ak.to_numpy(events[branch][sel])

    return sample, n_before


# ---------------------------------------------------------------------------
# Z peak estimate (histogram mode in 85–97 GeV)
# ---------------------------------------------------------------------------

def estimate_z_peak_mode(masses, window=(85.0, 97.0), n_bins=24):
    """
    Return the bin centre with the maximum count in the given mass window.
    This is a robust, assumption-free peak locator for diagnostic purposes.
    """
    m = masses[(masses >= window[0]) & (masses <= window[1])]
    if len(m) < 10:
        return None
    counts, edges = np.histogram(m, bins=n_bins, range=window)
    peak_idx  = np.argmax(counts)
    peak_mode = 0.5 * (edges[peak_idx] + edges[peak_idx + 1])
    return float(peak_mode)


# ---------------------------------------------------------------------------
# Figures
# ---------------------------------------------------------------------------

def plot_inclusive(sample, out_path):
    """Two-panel: dimuon mass (with mode peak annotation) + inclusive MET."""
    fig, axes = plt.subplots(1, 2, figsize=(12, 5))

    # -- panel 1: dimuon mass
    ax = axes[0]
    masses = sample["dimuon_mass"]
    ax.hist(masses, bins=40, range=(75, 105), histtype="stepfilled",
            color="steelblue", alpha=0.75, label=f"N={len(masses):,}")
    mode = estimate_z_peak_mode(masses)
    if mode is not None:
        ax.axvline(mode, color="red", ls="--", lw=1.5,
                   label=f"Mode = {mode:.1f} GeV")
    ax.axvline(91.1876, color="gray", ls=":", lw=1.0, alpha=0.7,
               label="PDG Z (91.2 GeV)")
    ax.set_xlabel("Dimuon invariant mass (GeV)", fontsize=11)
    ax.set_ylabel("Events / bin", fontsize=11)
    ax.set_title("Z→μμ mass (selected events)", fontsize=12)
    ax.legend(fontsize=9)

    # -- panel 2: inclusive MET
    ax = axes[1]
    met = sample["met_pt"]
    ax.hist(met, bins=MET_BINS, range=MET_RANGE, histtype="stepfilled",
            color="darkorange", alpha=0.75, label=f"N={len(met):,}")
    ax.set_xlabel("MET$_\\mathrm{pt}$ (GeV)", fontsize=11)
    ax.set_ylabel("Events / bin", fontsize=11)
    ax.set_title("Inclusive MET in Z→μμ", fontsize=12)
    ax.legend(fontsize=9)

    fig.suptitle("Phase 1 baseline — inclusive Z→μμ sample", fontsize=13, y=1.01)
    fig.tight_layout()
    fig.savefig(out_path, dpi=130, bbox_inches="tight")
    plt.close(fig)
    print(f"Wrote: {out_path}")


def plot_met_by_zpt(sample, out_path):
    """MET distributions in 4 pT(Z) bins (protocol §4 binning)."""
    zpt = sample["dimuon_pt"]
    met = sample["met_pt"]

    n_bins  = len(ZPT_LABELS)
    fig, axes = plt.subplots(1, n_bins, figsize=(4 * n_bins, 4.5), sharey=True)

    colors = ["#2166ac", "#4dac26", "#d6604d", "#762a83"]

    for i, (lo, hi, label) in enumerate(
        zip(ZPT_EDGES[:-1], ZPT_EDGES[1:], ZPT_LABELS)
    ):
        mask = (zpt >= lo) & (zpt < hi)
        n    = int(mask.sum())
        ax   = axes[i]
        if n > 0:
            ax.hist(met[mask], bins=MET_BINS, range=MET_RANGE,
                    histtype="stepfilled", color=colors[i], alpha=0.75,
                    label=f"N={n:,}")
            med = float(np.median(met[mask]))
            ax.axvline(med, color="k", ls="--", lw=1.2, label=f"Median={med:.0f}")
        ax.set_title(f"{label} GeV", fontsize=10)
        ax.set_xlabel("MET$_\\mathrm{pt}$ (GeV)", fontsize=9)
        if i == 0:
            ax.set_ylabel("Events / bin", fontsize=9)
        ax.legend(fontsize=8)

    fig.suptitle("MET stratified by pT(Z)", fontsize=13)
    fig.tight_layout()
    fig.savefig(out_path, dpi=130, bbox_inches="tight")
    plt.close(fig)
    print(f"Wrote: {out_path}")


def plot_met_by_njets(sample, out_path):
    """MET distributions by jet multiplicity: 0-jet / 1-jet / ≥2-jet."""
    nj  = sample["jet_mult"]
    met = sample["met_pt"]

    labels  = ["0 jets (pT>30 GeV)", "1 jet", "≥2 jets"]
    masks   = [nj == 0, nj == 1, nj >= 2]
    colors  = ["#1b7837", "#762a83", "#b35806"]

    fig, axes = plt.subplots(1, 3, figsize=(13, 4.5), sharey=True)

    for i, (label, mask, color) in enumerate(zip(labels, masks, colors)):
        n  = int(mask.sum())
        ax = axes[i]
        if n > 0:
            ax.hist(met[mask], bins=MET_BINS, range=MET_RANGE,
                    histtype="stepfilled", color=color, alpha=0.75,
                    label=f"N={n:,}")
            med = float(np.median(met[mask]))
            ax.axvline(med, color="k", ls="--", lw=1.2, label=f"Median={med:.0f}")
        ax.set_title(label, fontsize=10)
        ax.set_xlabel("MET$_\\mathrm{pt}$ (GeV)", fontsize=9)
        if i == 0:
            ax.set_ylabel("Events / bin", fontsize=9)
        ax.legend(fontsize=8)

    fig.suptitle("MET stratified by jet multiplicity", fontsize=13)
    fig.tight_layout()
    fig.savefig(out_path, dpi=130, bbox_inches="tight")
    plt.close(fig)
    print(f"Wrote: {out_path}")


# ---------------------------------------------------------------------------
# Event count table
# ---------------------------------------------------------------------------

def print_count_table(sample, n_before):
    nsel = len(sample["met_pt"])
    zpt  = sample["dimuon_pt"]
    nj   = sample["jet_mult"]

    rows = [
        ("Inclusive Z→μμ selected",    nsel,               nsel),
        ("  pT(Z) < 10 GeV",           (zpt < 10).sum(),   nsel),
        ("  10 ≤ pT(Z) < 30 GeV",      ((zpt >= 10) & (zpt < 30)).sum(), nsel),
        ("  30 ≤ pT(Z) < 60 GeV",      ((zpt >= 30) & (zpt < 60)).sum(), nsel),
        ("  pT(Z) ≥ 60 GeV",           (zpt >= 60).sum(),  nsel),
        ("  0 jets (pT>30 GeV)",        (nj == 0).sum(),    nsel),
        ("  1 jet",                     (nj == 1).sum(),    nsel),
        ("  ≥2 jets",                   (nj >= 2).sum(),    nsel),
    ]

    print(f"\n{'='*60}")
    print(f"  Events processed:  {n_before:>10,}")
    print(f"  Events selected:   {nsel:>10,}  ({100*nsel/max(n_before,1):.1f}%)")
    print(f"{'='*60}")
    print(f"  {'Stratum':<35}  {'N':>8}  {'frac':>6}")
    print(f"  {'-'*35}  {'-'*8}  {'-'*6}")
    for label, n, base in rows:
        frac = 100.0 * n / max(base, 1)
        print(f"  {label:<35}  {n:>8,}  {frac:>5.1f}%")
    print(f"{'='*60}\n")

    return {row[0].strip(): int(row[1]) for row in rows}


# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description="Phase 1 baseline: inclusive + stratified MET in Z→μμ"
    )
    parser.add_argument("--file",       required=True,
                        help="Path to NanoAOD ROOT file")
    parser.add_argument("--max-events", type=int, default=None,
                        help="Cap on events read (default: all)")
    args = parser.parse_args()

    # ------------------------------------------------------------------
    # Open and check branches
    # ------------------------------------------------------------------
    print(f"Opening: {args.file}")
    try:
        f = uproot.open(args.file)
    except Exception as exc:
        print(f"ERROR: cannot open file — {exc}")
        sys.exit(1)

    if "Events" not in f:
        print("ERROR: 'Events' tree not found.")
        sys.exit(1)

    avail = set(f["Events"].keys())

    # resolve ID / iso branches
    id_branch  = next(
        (b for b in ["Muon_tightId", "Muon_mediumId", "Muon_looseId"] if b in avail),
        None,
    )
    iso_branch = next(
        (b for b in ["Muon_pfRelIso04_all", "Muon_pfRelIso03_all",
                     "Muon_miniPFRelIso_all"] if b in avail),
        None,
    )
    if id_branch is None or iso_branch is None:
        print("ERROR: no usable muon ID or isolation branch found.")
        sys.exit(1)
    if id_branch != "Muon_tightId":
        print(f"  [warn] Using {id_branch} as Muon_tightId substitute.")
    if iso_branch != "Muon_pfRelIso04_all":
        print(f"  [warn] Using {iso_branch} as Muon_pfRelIso04_all substitute.")

    branches_needed = [
        "Muon_pt", "Muon_eta", "Muon_phi", "Muon_charge", "Muon_mass",
        "nMuon", id_branch, iso_branch,
        "MET_pt", "Jet_pt", "Jet_eta", "Jet_phi",
    ]
    missing = [b for b in branches_needed if b not in avail]
    if missing:
        print(f"ERROR: missing required branches: {missing}")
        sys.exit(1)

    # ------------------------------------------------------------------
    # Read arrays
    # ------------------------------------------------------------------
    n_entries = f["Events"].num_entries
    cap = args.max_events if args.max_events else n_entries
    print(f"Reading {min(cap, n_entries):,} / {n_entries:,} events…")

    events = f["Events"].arrays(
        branches_needed,
        entry_stop=args.max_events,
        library="ak",
    )
    f.close()

    # normalise substituted branch names for event_selection utilities
    if id_branch != "Muon_tightId":
        events["Muon_tightId"] = events[id_branch]
    if iso_branch != "Muon_pfRelIso04_all":
        events["Muon_pfRelIso04_all"] = events[iso_branch]

    # ------------------------------------------------------------------
    # Build Z→μμ sample
    # ------------------------------------------------------------------
    print("Applying Z→μμ selection…")
    sample, n_before = build_zmumu_sample(events)
    counts = print_count_table(sample, n_before)

    # Z peak mode estimate (diagnostic)
    mode = estimate_z_peak_mode(sample["dimuon_mass"])
    if mode is not None:
        print(f"  Z peak mode (histogram, 85–97 GeV):  {mode:.2f} GeV")
    else:
        print("  Z peak mode: insufficient statistics in 85–97 GeV window")

    if len(sample["met_pt"]) == 0:
        print("ERROR: no events survived selection — check input file.")
        sys.exit(1)

    # ------------------------------------------------------------------
    # Produce figures
    # ------------------------------------------------------------------
    FIGURES_INCL.mkdir(parents=True, exist_ok=True)
    FIGURES_STRAT.mkdir(parents=True, exist_ok=True)
    OUTPUTS_DIR.mkdir(parents=True, exist_ok=True)

    plot_inclusive(sample,     FIGURES_INCL  / "met_inclusive.png")
    plot_met_by_zpt(sample,    FIGURES_STRAT / "met_by_zpt.png")
    plot_met_by_njets(sample,  FIGURES_STRAT / "met_by_njets.png")

    # ------------------------------------------------------------------
    # Save counts
    # ------------------------------------------------------------------
    counts_out = {
        "input_file":      args.file,
        "events_processed": n_before,
        "events_selected":  int(len(sample["met_pt"])),
        "z_peak_mode_GeV": mode,
        "strata":          counts,
    }
    counts_path = OUTPUTS_DIR / "baseline_counts.json"
    with open(counts_path, "w") as fh:
        json.dump(counts_out, fh, indent=2)
    print(f"Wrote: {counts_path}")

    print("\nDone.")


if __name__ == "__main__":
    main()
