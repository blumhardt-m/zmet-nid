#!/usr/bin/env python3
"""
fit_minimal_nuisance_families.py — minimal two-family identifiability test
===========================================================================

Demonstrates the core identifiability question:

  Given a target inclusive MET distribution, can two structurally distinct
  nuisance families both explain it equally well in the inclusive view
  while making different predictions in stratified views?

Setup
-----
A target distribution (h_target) is constructed by injecting Family A at a
fixed alpha (default 0.15).  Both families are then fit to h_target using only
the INCLUSIVE histogram.  Stratified predictions from the best-fit parameters
are then compared — without refitting.

  If chi²_A ≈ chi²_B on inclusive MET  →  inclusive is non-identifiable.
  If the families diverge in pT(Z) strata  →  stratified tests can separate them.

Family A — recoil-like scale distortion (event-level MET warp)
  MET' = MET * (1 + alpha * clip(zpt, 0, 100) / 100)
  Effect is zero at low pT(Z), maximal at pT(Z) = 100 GeV.

Family B — topology-mixture reweighting (changes sample composition)
  Event weights: 1 / 1+beta / 1+2*beta for 0-jet / 1-jet / >=2-jet.
  MET values are unchanged; the relative jet-stratum contribution shifts.

Outputs
-------
  figures/inclusive_fit/met_familyA_vs_familyB.png
  figures/stratified_tests/met_familyA_vs_familyB_by_zpt.png
  data/processed/phase1/minimal_fit_summary.json

Usage
-----
  python code/analysis/fit_minimal_nuisance_families.py --file data/raw/test.root
  python code/analysis/fit_minimal_nuisance_families.py --file data/raw/test.root \\
        --inject-alpha 0.20 --max-events 200000
"""

import sys
import json
import argparse
from pathlib import Path

import numpy as np
from scipy.stats import chi2 as chi2_dist
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))   # code/
sys.path.insert(0, str(Path(__file__).resolve().parent))       # code/analysis/

from phase1_baseline import build_zmumu_sample

FIGURES_INCL  = Path("figures/inclusive_fit")
FIGURES_STRAT = Path("figures/stratified_tests")
OUT_DIR       = Path("data/processed/phase1")

# MET bins (variable-width; finer at low MET)
MET_BINS = np.array([0, 5, 10, 15, 20, 30, 40, 50, 60, 80, 100, 150, 200],
                    dtype=float)

# pT(Z) strata — protocol §4
ZPT_EDGES  = [0, 10, 30, 60, float("inf")]
ZPT_LABELS = ["pT(Z)<10", "10≤pT(Z)<30", "30≤pT(Z)<60", "pT(Z)≥60"]

# Parameter grids
ALPHAS  = np.linspace(-0.3, 0.3, 61)
BETAS   = np.linspace(-0.4, 0.4, 81)
GAMMAS  = np.linspace(-3.0, 3.0, 121)   # Family C: unclustered-energy scale factor


# ---------------------------------------------------------------------------
# Family transforms
# ---------------------------------------------------------------------------

def warp_family_a(met, zpt, alpha):
    """Scale MET by a pT(Z)-dependent factor.  Returns (met', unit weights)."""
    scale = 1.0 + alpha * np.clip(zpt, 0.0, 100.0) / 100.0
    return met * scale, np.ones(len(met))


def warp_family_b(met, njets, beta):
    """
    Reweight events by jet multiplicity; MET values unchanged.
    Returns (met, weights) or (None, None) if weights go non-positive.
    """
    w = np.ones(len(met), dtype=float)
    w[njets == 1] *= (1.0 + beta)
    w[njets >= 2] *= (1.0 + 2.0 * beta)
    if w.min() <= 0.0:
        return None, None
    w *= len(w) / w.sum()
    return met, w


def warp_family_c(met_pt, met_phi, delta_x, delta_y, gamma):
    """
    Shift the MET vector by gamma * CMS unclustered-energy delta vector.
      MET_x' = MET_pt * cos(phi) + gamma * delta_x
      MET_y' = MET_pt * sin(phi) + gamma * delta_y
    gamma = 1 corresponds to the full CMS MetUnclustEnUp shift.
    Returns (met_pt_new, unit weights).
    """
    met_x = met_pt * np.cos(met_phi) + gamma * delta_x
    met_y = met_pt * np.sin(met_phi) + gamma * delta_y
    return np.sqrt(met_x ** 2 + met_y ** 2), np.ones(len(met_pt))


# ---------------------------------------------------------------------------
# Histogram and fit helpers
# ---------------------------------------------------------------------------

def histo(x, weights=None, bins=None):
    if bins is None:
        bins = MET_BINS
    h, _ = np.histogram(x, bins=bins, weights=weights)
    return h.astype(float)


def chisq(obs, model):
    """χ² with (obs+1) denominator — handles empty bins without blowup."""
    return float(np.sum((obs - model) ** 2 / (obs + 1.0)))


def fit_family_a(met, zpt, h_target):
    best_alpha, best_chi2 = 0.0, np.inf
    for alpha in ALPHAS:
        met_w, ww = warp_family_a(met, zpt, alpha)
        c = chisq(h_target, histo(met_w, ww))
        if c < best_chi2:
            best_chi2, best_alpha = c, alpha
    return best_alpha, best_chi2


def fit_family_b(met, njets, h_target):
    best_beta, best_chi2 = 0.0, np.inf
    for beta in BETAS:
        met_w, ww = warp_family_b(met, njets, beta)
        if met_w is None:
            continue
        c = chisq(h_target, histo(met_w, ww))
        if c < best_chi2:
            best_chi2, best_beta = c, beta
    return best_beta, best_chi2


def fit_family_c(met_pt, met_phi, delta_x, delta_y, h_target):
    best_gamma, best_chi2 = 0.0, np.inf
    for gamma in GAMMAS:
        met_w, ww = warp_family_c(met_pt, met_phi, delta_x, delta_y, gamma)
        c = chisq(h_target, histo(met_w, ww))
        if c < best_chi2:
            best_chi2, best_gamma = c, gamma
    return best_gamma, best_chi2


def compute_stratum_chi2(met_target, met_a, wa, met_b, wb, mask):
    """Return (chi2_A, chi2_B) evaluated within the events selected by mask."""
    h_t = histo(met_target[mask])
    h_a = histo(met_a[mask], wa[mask])
    h_b = histo(met_b[mask], wb[mask])
    return chisq(h_t, h_a), chisq(h_t, h_b)


def compute_binning_robustness(target_met, met, zpt, njets, alpha_best, beta_best):
    """
    Re-evaluate the inclusive fit quality under three nearby MET binnings.
    Returns a dict keyed by binning name with chi2_A, chi2_B, p_value_B per binning.
    dof convention: (number of histogram bins) - 1 fitted parameter.
    """
    binning_sets = {
        "baseline": np.array([0, 5, 10, 15, 20, 30, 40, 50, 60, 80, 100, 150, 200],
                             dtype=float),
        "coarse":   np.array([0, 10, 20, 30, 40, 60, 80, 100, 150, 200],
                             dtype=float),
        "tail":     np.array([0, 5, 10, 15, 20, 30, 40, 50, 70, 90, 120, 160, 200],
                             dtype=float),
    }
    met_a, wa = warp_family_a(met, zpt, alpha_best)
    met_b, wb = warp_family_b(met, njets, beta_best)

    results = {}
    for name, bins in binning_sets.items():
        n_bins = len(bins) - 1
        dof    = n_bins - 1          # one fitted parameter (beta)
        h_t    = histo(target_met, bins=bins)
        h_a    = histo(met_a, wa,  bins=bins)
        h_b    = histo(met_b, wb,  bins=bins)
        c_a    = chisq(h_t, h_a)
        c_b    = chisq(h_t, h_b)
        p_b    = float(chi2_dist.sf(c_b, dof))
        results[name] = {
            "n_bins":    n_bins,
            "dof":       dof,
            "chi2_A":    round(c_a, 4),
            "chi2_B":    round(c_b, 4),
            "p_value_B": round(p_b, 6),
        }
    return results


# ---------------------------------------------------------------------------
# Plotting
# ---------------------------------------------------------------------------

def plot_inclusive(h_baseline, h_target, h_a, h_b,
                   alpha_inj, alpha_fit, beta_fit,
                   chi2_a, chi2_b, out_path,
                   h_c=None, gamma_fit=None, chi2_c=None):
    """
    Two-panel inclusive figure:
      Top    : target (injected) vs best-fit A vs best-fit B [vs best-fit C]
      Bottom : ratio model / target
    """
    bw = MET_BINS[1:] - MET_BINS[:-1]

    def dens(h): return h / bw

    fig, (ax1, ax2) = plt.subplots(
        2, 1, figsize=(8, 7),
        gridspec_kw={"height_ratios": [3, 1]},
        sharex=True,
    )

    ax1.step(MET_BINS[:-1], dens(h_baseline), where="post",
             color="gray",    lw=1.2, ls="-",  alpha=0.5,
             label="Baseline (no injection)")
    ax1.step(MET_BINS[:-1], dens(h_target),   where="post",
             color="black",   lw=2.0, ls="-",
             label=f"Target  (inject α={alpha_inj:+.2f})")
    ax1.step(MET_BINS[:-1], dens(h_a),        where="post",
             color="#d62728", lw=1.5, ls="--",
             label=f"Family A  α={alpha_fit:+.3f}  χ²={chi2_a:.1f}")
    ax1.step(MET_BINS[:-1], dens(h_b),        where="post",
             color="#1f77b4", lw=1.5, ls=":",
             label=f"Family B  β={beta_fit:+.3f}  χ²={chi2_b:.1f}")
    if h_c is not None:
        ax1.step(MET_BINS[:-1], dens(h_c), where="post",
                 color="#2ca02c", lw=1.5, ls="-.",
                 label=f"Family C  γ={gamma_fit:+.3f}  χ²={chi2_c:.1f}")

    ax1.set_ylabel("Events / GeV", fontsize=11)
    ax1.set_yscale("log")
    ax1.legend(fontsize=9)
    ax1.set_title(
        f"Inclusive MET — all families fitted to α={alpha_inj:+.2f} injection",
        fontsize=12,
    )

    safe = h_target + 1.0
    ax2.step(MET_BINS[:-1], h_a / safe, where="post",
             color="#d62728", lw=1.5, ls="--", label="A / target")
    ax2.step(MET_BINS[:-1], h_b / safe, where="post",
             color="#1f77b4", lw=1.5, ls=":",  label="B / target")
    if h_c is not None:
        ax2.step(MET_BINS[:-1], h_c / safe, where="post",
                 color="#2ca02c", lw=1.5, ls="-.", label="C / target")
    ax2.axhline(1.0, color="black", lw=0.8)
    ax2.set_ylim(0.7, 1.3)
    ax2.set_xlabel("MET$_\\mathrm{pt}$ (GeV)", fontsize=11)
    ax2.set_ylabel("Model / target", fontsize=10)
    ax2.legend(fontsize=8)

    fig.tight_layout()
    fig.savefig(out_path, dpi=130, bbox_inches="tight")
    plt.close(fig)
    print(f"Wrote: {out_path}")


def plot_stratified_zpt(met, zpt, njets,
                        alpha_inj, alpha_fit, beta_fit, out_path,
                        met_phi=None, delta_x=None, delta_y=None, gamma_fit=None):
    """
    4-panel pT(Z) strata: target vs best-fit A vs best-fit B [vs best-fit C].
    Families were fitted to inclusive only — divergence here is the signal.
    """
    met_target, _ = warp_family_a(met, zpt, alpha_inj)
    met_a,      wa = warp_family_a(met, zpt, alpha_fit)
    met_b,      wb = warp_family_b(met, njets, beta_fit)
    has_c = (met_phi is not None) and (gamma_fit is not None)
    if has_c:
        met_c, wc = warp_family_c(met, met_phi, delta_x, delta_y, gamma_fit)
    bw = MET_BINS[1:] - MET_BINS[:-1]

    fig, axes = plt.subplots(1, 4, figsize=(16, 4.8), sharey=False)

    for i, (lo, hi, label) in enumerate(
        zip(ZPT_EDGES[:-1], ZPT_EDGES[1:], ZPT_LABELS)
    ):
        mask = (zpt >= lo) & (zpt < hi)
        ax   = axes[i]
        n    = int(mask.sum())

        if n < 5:
            ax.set_title(label)
            continue

        h_t = histo(met_target[mask])
        h_a = histo(met_a[mask], wa[mask])
        h_b = histo(met_b[mask], wb[mask] if wb is not None else None)

        ax.step(MET_BINS[:-1], h_t / bw, where="post",
                color="black",   lw=2.0, label=f"Target N={n:,}")
        ax.step(MET_BINS[:-1], h_a / bw, where="post",
                color="#d62728", lw=1.5, ls="--",
                label=f"A α={alpha_fit:+.2f}")
        ax.step(MET_BINS[:-1], h_b / bw, where="post",
                color="#1f77b4", lw=1.5, ls=":",
                label=f"B β={beta_fit:+.2f}")
        if has_c:
            h_c = histo(met_c[mask], wc[mask])
            ax.step(MET_BINS[:-1], h_c / bw, where="post",
                    color="#2ca02c", lw=1.5, ls="-.",
                    label=f"C γ={gamma_fit:+.2f}")

        ax.set_yscale("log")
        ax.set_title(label + " GeV", fontsize=10)
        ax.set_xlabel("MET$_\\mathrm{pt}$ (GeV)", fontsize=9)
        if i == 0:
            ax.set_ylabel("Events / GeV", fontsize=9)
        ax.legend(fontsize=7)

    fig.suptitle(
        f"Stratified MET by pT(Z)  —  families fitted to inclusive only  "
        f"(inject α={alpha_inj:+.2f})",
        fontsize=12,
    )
    fig.tight_layout()
    fig.savefig(out_path, dpi=130, bbox_inches="tight")
    plt.close(fig)
    print(f"Wrote: {out_path}")


def plot_stratified_njets(met, zpt, njets,
                          alpha_inj, alpha_fit, beta_fit, out_path,
                          met_phi=None, delta_x=None, delta_y=None, gamma_fit=None):
    """
    3-panel jet-multiplicity strata: target vs best-fit A vs best-fit B [vs best-fit C].
    """
    met_target, _  = warp_family_a(met, zpt, alpha_inj)
    met_a,      wa = warp_family_a(met, zpt, alpha_fit)
    met_b,      wb = warp_family_b(met, njets, beta_fit)
    has_c = (met_phi is not None) and (gamma_fit is not None)
    if has_c:
        met_c, wc = warp_family_c(met, met_phi, delta_x, delta_y, gamma_fit)
    bw = MET_BINS[1:] - MET_BINS[:-1]

    nj_labels = ["0 jets (pT>30 GeV)", "1 jet", "≥2 jets"]
    nj_masks  = [njets == 0, njets == 1, njets >= 2]
    colors_b  = ["#1b7837", "#762a83", "#b35806"]

    fig, axes = plt.subplots(1, 3, figsize=(13, 4.8), sharey=False)

    for i, (label, mask, col) in enumerate(zip(nj_labels, nj_masks, colors_b)):
        ax = axes[i]
        n  = int(mask.sum())
        if n < 5:
            ax.set_title(label)
            continue

        h_t = histo(met_target[mask])
        h_a = histo(met_a[mask], wa[mask])
        h_b = histo(met_b[mask], wb[mask])

        ax.step(MET_BINS[:-1], h_t / bw, where="post",
                color="black",   lw=2.0, label=f"Target N={n:,}")
        ax.step(MET_BINS[:-1], h_a / bw, where="post",
                color="#d62728", lw=1.5, ls="--",
                label=f"A α={alpha_fit:+.2f}")
        ax.step(MET_BINS[:-1], h_b / bw, where="post",
                color=col,       lw=1.5, ls=":",
                label=f"B β={beta_fit:+.2f}")
        if has_c:
            h_c = histo(met_c[mask], wc[mask])
            ax.step(MET_BINS[:-1], h_c / bw, where="post",
                    color="#2ca02c", lw=1.5, ls="-.",
                    label=f"C γ={gamma_fit:+.2f}")

        ax.set_yscale("log")
        ax.set_title(label, fontsize=10)
        ax.set_xlabel("MET$_\\mathrm{pt}$ (GeV)", fontsize=9)
        if i == 0:
            ax.set_ylabel("Events / GeV", fontsize=9)
        ax.legend(fontsize=7)

    fig.suptitle(
        f"Stratified MET by jet multiplicity  —  families fitted to inclusive only  "
        f"(inject α={alpha_inj:+.2f})",
        fontsize=12,
    )
    fig.tight_layout()
    fig.savefig(out_path, dpi=130, bbox_inches="tight")
    plt.close(fig)
    print(f"Wrote: {out_path}")


# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description="Minimal two-family identifiability fit"
    )
    parser.add_argument("--file",         required=True,
                        help="NanoAOD ROOT file")
    parser.add_argument("--max-events",   type=int, default=None)
    parser.add_argument("--inject-alpha", type=float, default=0.15,
                        help="Family A alpha used to construct the target "
                             "distribution (default: 0.15)")
    parser.add_argument("--report-binning-robustness", action="store_true",
                        help="Repeat inclusive fit under coarse and tail-emphasised "
                             "MET binnings to verify the degeneracy is not a "
                             "binning artefact.")
    args = parser.parse_args()

    # ------------------------------------------------------------------
    # Load events
    # ------------------------------------------------------------------
    import uproot
    import awkward as ak
    import vector
    vector.register_awkward()

    print(f"Opening: {args.file}")
    f     = uproot.open(args.file)
    avail = set(f["Events"].keys())

    id_b  = next(b for b in ["Muon_tightId",  "Muon_mediumId",  "Muon_looseId"]
                 if b in avail)
    iso_b = next(b for b in ["Muon_pfRelIso04_all", "Muon_pfRelIso03_all",
                              "Muon_miniPFRelIso_all"] if b in avail)

    branches = [
        "Muon_pt", "Muon_eta", "Muon_phi", "Muon_charge", "Muon_mass",
        "nMuon", id_b, iso_b,
        "MET_pt", "Jet_pt", "Jet_eta", "Jet_phi",
    ]
    c_opt = ["MET_phi", "MET_MetUnclustEnUpDeltaX", "MET_MetUnclustEnUpDeltaY"]
    branches += [b for b in c_opt if b in avail]
    events = f["Events"].arrays(branches, entry_stop=args.max_events, library="ak")
    f.close()

    if id_b  != "Muon_tightId":        events["Muon_tightId"]        = events[id_b]
    if iso_b != "Muon_pfRelIso04_all": events["Muon_pfRelIso04_all"] = events[iso_b]

    print("Applying Z→μμ selection…")
    sample, n_before = build_zmumu_sample(events)
    met   = sample["met_pt"].astype(float)
    zpt   = sample["dimuon_pt"].astype(float)
    njets = sample["jet_mult"].astype(int)
    # Family C inputs (present only if unclustered-energy branches exist)
    has_family_c = all(k in sample for k in ("met_phi", "met_delta_x", "met_delta_y"))
    if has_family_c:
        met_phi = sample["met_phi"].astype(float)
        delta_x = sample["met_delta_x"].astype(float)
        delta_y = sample["met_delta_y"].astype(float)
        print("  Family C branches found: MET_phi, MET_MetUnclustEnUpDeltaX/Y")
    else:
        met_phi = delta_x = delta_y = None
        print("  Family C branches not found — skipping Family C.")
    print(f"  Selected {len(met):,} / {n_before:,} events.")

    # ------------------------------------------------------------------
    # Build baseline and target distributions
    # ------------------------------------------------------------------
    h_baseline = histo(met)                                 # raw observed MET

    alpha_inj   = args.inject_alpha
    met_t, _    = warp_family_a(met, zpt, alpha_inj)        # injected target
    h_target    = histo(met_t)
    print(f"\nInjection: Family A with α={alpha_inj:+.2f}")
    print(f"  Baseline total:  {h_baseline.sum():.0f}")
    print(f"  Target total:    {h_target.sum():.0f}")

    # ------------------------------------------------------------------
    # Fit both families to h_target
    # ------------------------------------------------------------------
    print("\nGrid-searching Family A (alpha)…")
    alpha_fit, chi2_a = fit_family_a(met, zpt, h_target)
    met_a, wa         = warp_family_a(met, zpt, alpha_fit)
    h_a               = histo(met_a, wa)
    print(f"  Best alpha = {alpha_fit:+.3f}   χ² = {chi2_a:.2f}")

    print("Grid-searching Family B (beta)…")
    beta_fit, chi2_b  = fit_family_b(met, njets, h_target)
    met_b, wb         = warp_family_b(met, njets, beta_fit)
    h_b               = histo(met_b, wb)
    print(f"  Best beta  = {beta_fit:+.3f}   χ² = {chi2_b:.2f}")

    if has_family_c:
        print("Grid-searching Family C (gamma)…")
        gamma_fit, chi2_c = fit_family_c(met, met_phi, delta_x, delta_y, h_target)
        met_c, wc         = warp_family_c(met, met_phi, delta_x, delta_y, gamma_fit)
        h_c               = histo(met_c, wc)
        dof_incl_c        = (len(MET_BINS) - 1) - 1
        p_value_c         = float(chi2_dist.sf(chi2_c, dof_incl_c))
        print(f"  Best gamma = {gamma_fit:+.3f}   χ² = {chi2_c:.2f}   p = {p_value_c:.4f}")
    else:
        gamma_fit = chi2_c = p_value_c = None
        met_c = wc = h_c = None

    ratio = chi2_b / max(chi2_a, 1e-9)

    # Formal p-value for Family B on the inclusive histogram.
    # dof = (number of bins) - 1 fitted parameter (beta).
    # MET_BINS has 13 edges → 12 histogram bins → dof = 11.
    dof_incl  = (len(MET_BINS) - 1) - 1
    p_value_b = float(chi2_dist.sf(chi2_b, dof_incl))

    print(f"\n  χ²(B) / χ²(A) = {ratio:.2f}")
    print(f"  Inclusive p-value  χ²(B)={chi2_b:.2f}, dof={dof_incl}, p={p_value_b:.4f}")
    if p_value_b > 0.05:
        print("  → Inclusive MET provides no statistical basis to reject Family B.")
        print("    The inclusive projection is non-identifiable at this working point.")
    else:
        print("  → Inclusive MET already rejects Family B (p < 0.05).")

    # ------------------------------------------------------------------
    # Per-stratum χ²
    # ------------------------------------------------------------------
    met_target_arr, _ = warp_family_a(met, zpt, alpha_inj)

    zpt_stratum_chi2 = {}
    print("\nPer-stratum χ² (pT(Z) bins):")
    print(f"  {'Stratum':<22}  {'N':>8}  {'χ²(A)':>8}  {'χ²(B)':>8}  {'Δχ²':>10}")
    print(f"  {'-'*22}  {'-'*8}  {'-'*8}  {'-'*8}  {'-'*10}")
    for lo, hi, lab in zip(ZPT_EDGES[:-1], ZPT_EDGES[1:], ZPT_LABELS):
        mask = (zpt >= lo) & (zpt < hi)
        n = int(mask.sum())
        if n < 5:
            continue
        c_a, c_b = compute_stratum_chi2(met_target_arr, met_a, wa, met_b, wb, mask)
        delta = c_b - c_a
        print(f"  {lab:<22}  {n:>8,}  {c_a:>8.2f}  {c_b:>8.2f}  {delta:>10.1f}")
        zpt_stratum_chi2[lab] = {"n": n, "chi2_A": round(c_a, 4), "chi2_B": round(c_b, 4),
                                  "delta_chi2": round(delta, 4)}

    nj_stratum_chi2 = {}
    nj_defs = [("0-jet", njets == 0), ("1-jet", njets == 1), ("geq2jet", njets >= 2)]
    print("\nPer-stratum χ² (jet-multiplicity bins):")
    print(f"  {'Stratum':<22}  {'N':>8}  {'χ²(A)':>8}  {'χ²(B)':>8}  {'Δχ²':>10}")
    print(f"  {'-'*22}  {'-'*8}  {'-'*8}  {'-'*8}  {'-'*10}")
    for lab, mask in nj_defs:
        n = int(mask.sum())
        if n < 5:
            continue
        c_a, c_b = compute_stratum_chi2(met_target_arr, met_a, wa, met_b, wb, mask)
        delta = c_b - c_a
        print(f"  {lab:<22}  {n:>8,}  {c_a:>8.2f}  {c_b:>8.2f}  {delta:>10.1f}")
        nj_stratum_chi2[lab] = {"n": n, "chi2_A": round(c_a, 4), "chi2_B": round(c_b, 4),
                                 "delta_chi2": round(delta, 4)}

    # Family C per-stratum chi2 (if available)
    if has_family_c:
        print("\nPer-stratum χ² for Family C (pT(Z) bins):")
        print(f"  {'Stratum':<22}  {'N':>8}  {'χ²(C)':>8}  {'Δχ²(C-A)':>12}")
        print(f"  {'-'*22}  {'-'*8}  {'-'*8}  {'-'*12}")
        for lo, hi, lab in zip(ZPT_EDGES[:-1], ZPT_EDGES[1:], ZPT_LABELS):
            mask = (zpt >= lo) & (zpt < hi)
            n = int(mask.sum())
            if n < 5:
                continue
            h_t = histo(met_target_arr[mask])
            c_a_s = chisq(h_t, histo(met_a[mask], wa[mask]))
            c_c_s = chisq(h_t, histo(met_c[mask], wc[mask]))
            delta_c = c_c_s - c_a_s
            print(f"  {lab:<22}  {n:>8,}  {c_c_s:>8.2f}  {delta_c:>12.1f}")
            zpt_stratum_chi2[lab]["chi2_C"] = round(c_c_s, 4)
            zpt_stratum_chi2[lab]["delta_chi2_C"] = round(delta_c, 4)

        print("\nPer-stratum χ² for Family C (jet-multiplicity bins):")
        print(f"  {'Stratum':<22}  {'N':>8}  {'χ²(C)':>8}  {'Δχ²(C-A)':>12}")
        print(f"  {'-'*22}  {'-'*8}  {'-'*8}  {'-'*12}")
        for lab, mask in nj_defs:
            n = int(mask.sum())
            if n < 5:
                continue
            h_t = histo(met_target_arr[mask])
            c_a_s = chisq(h_t, histo(met_a[mask], wa[mask]))
            c_c_s = chisq(h_t, histo(met_c[mask], wc[mask]))
            delta_c = c_c_s - c_a_s
            print(f"  {lab:<22}  {n:>8,}  {c_c_s:>8.2f}  {delta_c:>12.1f}")
            nj_stratum_chi2[lab]["chi2_C"] = round(c_c_s, 4)
            nj_stratum_chi2[lab]["delta_chi2_C"] = round(delta_c, 4)

    # ------------------------------------------------------------------
    # Figures
    # ------------------------------------------------------------------
    FIGURES_INCL.mkdir(parents=True, exist_ok=True)
    FIGURES_STRAT.mkdir(parents=True, exist_ok=True)
    OUT_DIR.mkdir(parents=True, exist_ok=True)

    plot_inclusive(
        h_baseline, h_target, h_a, h_b,
        alpha_inj, alpha_fit, beta_fit, chi2_a, chi2_b,
        FIGURES_INCL / "met_familyA_vs_familyB.png",
        h_c=h_c, gamma_fit=gamma_fit, chi2_c=chi2_c,
    )
    plot_stratified_zpt(
        met, zpt, njets,
        alpha_inj, alpha_fit, beta_fit,
        FIGURES_STRAT / "met_familyA_vs_familyB_by_zpt.png",
        met_phi=met_phi, delta_x=delta_x, delta_y=delta_y, gamma_fit=gamma_fit,
    )
    plot_stratified_njets(
        met, zpt, njets,
        alpha_inj, alpha_fit, beta_fit,
        FIGURES_STRAT / "met_familyA_vs_familyB_by_njets.png",
        met_phi=met_phi, delta_x=delta_x, delta_y=delta_y, gamma_fit=gamma_fit,
    )

    # ------------------------------------------------------------------
    # JSON
    # ------------------------------------------------------------------
    summary = {
        "input_file":        args.file,
        "events_processed":  n_before,
        "events_selected":   int(len(met)),
        "injection":         {"family": "A", "alpha": alpha_inj},
        "fit": {
            "family_a": {"alpha": float(alpha_fit), "chi2_inclusive": chi2_a},
            "family_b": {"beta":  float(beta_fit),  "chi2_inclusive": chi2_b,
                         "dof_inclusive": dof_incl,
                         "p_value_inclusive": round(p_value_b, 6)},
            **({"family_c": {"gamma": float(gamma_fit), "chi2_inclusive": chi2_c,
                             "dof_inclusive": dof_incl_c,
                             "p_value_inclusive": round(p_value_c, 6)}}
               if has_family_c else {}),
        },
        "chi2_ratio_B_over_A": round(ratio, 3),
        "stratum_chi2": {
            "zpt":   zpt_stratum_chi2,
            "njets": nj_stratum_chi2,
        },
        "strata_counts": {
            "zpt": {
                lab: int(((zpt >= lo) & (zpt < hi)).sum())
                for lo, hi, lab in zip(ZPT_EDGES[:-1], ZPT_EDGES[1:], ZPT_LABELS)
            },
            "njets": {
                "0-jet":   int((njets == 0).sum()),
                "1-jet":   int((njets == 1).sum()),
                "geq2jet": int((njets >= 2).sum()),
            },
        },
    }
    # Binning robustness (optional)
    if args.report_binning_robustness:
        print("\nInclusive binning robustness:")
        robustness = compute_binning_robustness(
            met_target_arr, met, zpt, njets, alpha_fit, beta_fit
        )
        print(f"  {'Binning':<10}  {'n_bins':>6}  {'dof':>4}  "
              f"{'χ²(A)':>8}  {'χ²(B)':>8}  {'p(B)':>8}")
        print(f"  {'-'*10}  {'-'*6}  {'-'*4}  {'-'*8}  {'-'*8}  {'-'*8}")
        for name, r in robustness.items():
            print(f"  {name:<10}  {r['n_bins']:>6}  {r['dof']:>4}  "
                  f"{r['chi2_A']:>8.3f}  {r['chi2_B']:>8.3f}  {r['p_value_B']:>8.4f}")
        summary["binning_robustness"] = robustness

    out_json = OUT_DIR / "minimal_fit_summary.json"
    with open(out_json, "w") as fh:
        json.dump(summary, fh, indent=2)
    print(f"Wrote: {out_json}")

    print("\nDone.")


if __name__ == "__main__":
    main()
