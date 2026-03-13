#!/usr/bin/env python3
"""
null_check.py — empirical null distribution for the inclusive chi² statistic

Takes the target MET array (Family A at the injection working point), splits it
randomly into two halves 1000 times, computes the same chi² statistic between
the two halves, and reports where Family B's inclusive chi² = 3.13 sits within
that null distribution.

Usage
-----
  python code/analysis/null_check.py --file data/raw/test.root --inject-alpha 0.05
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

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
sys.path.insert(0, str(Path(__file__).resolve().parent))

from phase1_baseline import build_zmumu_sample
from fit_minimal_nuisance_families import MET_BINS, warp_family_a, histo, chisq

N_TRIALS = 1000
RNG_SEED = 42


def null_distribution(met_target, n_trials=N_TRIALS, rng=None):
    """
    Empirical null: bootstrap N events with replacement from met_target,
    histogram the bootstrap sample, and compute chi² against the original
    target histogram. Repeat n_trials times.

    This compares two N-event histograms (bootstrap vs. original), exactly
    matching the scale of the Family B vs. target comparison, in which Family
    B's weighted histogram also sums to N events.

    Returns array of chi² values.
    """
    if rng is None:
        rng = np.random.default_rng(RNG_SEED)

    n = len(met_target)
    h_target = histo(met_target)
    chi2_vals = np.empty(n_trials)

    for i in range(n_trials):
        idx_boot = rng.choice(n, size=n, replace=True)
        h_boot   = histo(met_target[idx_boot])
        chi2_vals[i] = chisq(h_target, h_boot)

    return chi2_vals


def main():
    parser = argparse.ArgumentParser(
        description="Empirical null distribution for inclusive chi²"
    )
    parser.add_argument("--file",         required=True)
    parser.add_argument("--max-events",   type=int, default=None)
    parser.add_argument("--inject-alpha", type=float, default=0.05)
    parser.add_argument("--n-trials",     type=int, default=N_TRIALS)
    args = parser.parse_args()

    # ------------------------------------------------------------------
    # Load events and build target MET array
    # ------------------------------------------------------------------
    import uproot, awkward as ak, vector
    vector.register_awkward()

    print(f"Opening: {args.file}")
    f     = uproot.open(args.file)
    avail = set(f["Events"].keys())

    id_b  = next(b for b in ["Muon_tightId", "Muon_mediumId", "Muon_looseId"]
                 if b in avail)
    iso_b = next(b for b in ["Muon_pfRelIso04_all", "Muon_pfRelIso03_all",
                              "Muon_miniPFRelIso_all"] if b in avail)

    branches = [
        "Muon_pt", "Muon_eta", "Muon_phi", "Muon_charge", "Muon_mass",
        "nMuon", id_b, iso_b, "MET_pt", "Jet_pt", "Jet_eta", "Jet_phi",
    ]
    events = f["Events"].arrays(branches, entry_stop=args.max_events, library="ak")
    f.close()

    if id_b  != "Muon_tightId":        events["Muon_tightId"]        = events[id_b]
    if iso_b != "Muon_pfRelIso04_all": events["Muon_pfRelIso04_all"] = events[iso_b]

    sample, n_before = build_zmumu_sample(events)
    met   = sample["met_pt"].astype(float)
    zpt   = sample["dimuon_pt"].astype(float)
    print(f"  Selected {len(met):,} / {n_before:,} events.")

    # Build target distribution
    alpha_inj  = args.inject_alpha
    met_target, _ = warp_family_a(met, zpt, alpha_inj)
    print(f"  Injection α = {alpha_inj:+.2f}   N(target) = {len(met_target):,}")

    # ------------------------------------------------------------------
    # Run null distribution
    # ------------------------------------------------------------------
    print(f"\nRunning {args.n_trials} random half-sample splits…")
    rng  = np.random.default_rng(RNG_SEED)
    null = null_distribution(met_target, n_trials=args.n_trials, rng=rng)

    # Observed Family B chi²
    chi2_b = 3.13   # from main analysis at alpha=0.05

    pct_above = 100.0 * (null >= chi2_b).mean()
    pct_rank  = 100.0 * (null < chi2_b).mean()    # percentile of chi2_b in null

    print(f"\nNull distribution summary ({args.n_trials} trials):")
    print(f"  min    = {null.min():.3f}")
    print(f"  median = {np.median(null):.3f}")
    print(f"  mean   = {null.mean():.3f}")
    print(f"  max    = {null.max():.3f}")
    print(f"  90th pct = {np.percentile(null, 90):.3f}")
    print(f"  95th pct = {np.percentile(null, 95):.3f}")
    print(f"\nFamily B chi²(B) = {chi2_b:.2f}")
    print(f"  Percentile in null: {pct_rank:.1f}th  ({pct_above:.1f}% of null >= chi2_b)")

    if pct_rank < 50:
        print("  → chi²(B) is BELOW the median of the null distribution.")
    elif pct_rank < 75:
        print("  → chi²(B) is between the median and 75th percentile of the null.")
    else:
        print("  → chi²(B) is in the upper tail of the null distribution.")

    # ------------------------------------------------------------------
    # Figure
    # ------------------------------------------------------------------
    Path("figures/inclusive_fit").mkdir(parents=True, exist_ok=True)
    fig, ax = plt.subplots(figsize=(7, 4.5))
    ax.hist(null, bins=40, color="#1f77b4", alpha=0.75,
            label=f"Null: random half-sample splits (N={args.n_trials})")
    ax.axvline(chi2_b, color="#d62728", lw=2.0, ls="--",
               label=f"χ²(B) = {chi2_b:.2f}  (Family B inclusive fit)")
    ax.axvline(np.median(null), color="gray", lw=1.2, ls=":",
               label=f"Null median = {np.median(null):.2f}")
    ax.set_xlabel("χ²  (same metric as main analysis)", fontsize=11)
    ax.set_ylabel("Count", fontsize=11)
    ax.set_title(
        f"Empirical null distribution: random half-sample χ²\n"
        f"(injection α = {alpha_inj:+.2f}, {args.n_trials} trials)",
        fontsize=11,
    )
    ax.legend(fontsize=9)
    fig.tight_layout()
    out = Path("figures/inclusive_fit/null_chi2_distribution.png")
    fig.savefig(str(out), dpi=130, bbox_inches="tight")
    plt.close(fig)
    print(f"\nWrote: {out}")

    # ------------------------------------------------------------------
    # JSON
    # ------------------------------------------------------------------
    Path("data/processed/phase1").mkdir(parents=True, exist_ok=True)
    result = {
        "injection_alpha":   alpha_inj,
        "n_trials":          args.n_trials,
        "chi2_family_b":     chi2_b,
        "null_min":          round(float(null.min()), 4),
        "null_median":       round(float(np.median(null)), 4),
        "null_mean":         round(float(null.mean()), 4),
        "null_max":          round(float(null.max()), 4),
        "null_p90":          round(float(np.percentile(null, 90)), 4),
        "null_p95":          round(float(np.percentile(null, 95)), 4),
        "chi2_b_percentile": round(pct_rank, 2),
        "frac_null_above":   round(pct_above / 100.0, 4),
    }
    out_json = Path("data/processed/phase1/null_check.json")
    with open(out_json, "w") as fh:
        json.dump(result, fh, indent=2)
    print(f"Wrote: {out_json}")
    print("\nDone.")


if __name__ == "__main__":
    main()
