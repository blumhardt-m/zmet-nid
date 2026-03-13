#!/usr/bin/env python3
"""
ZMET-NID Phase 1 — S7: Stratified Tests
========================================
Run KS tests per stratum (protocol v0.4 §4):
  - Pileup terciles: low/mid/high (3)
  - Jet multiplicity: 0, 1, ≥2 (3)
  - Z pT bins: [0,10), [10,30), [30,60), [60,∞) (4)

Primary strata: 3 × 3 × 4 = 36 combinations
But we test marginal strata for interpretability: 3 + 3 + 4 = 10

Bonferroni correction: reject if p < 0.01 / N_strata

A stratum is "separating" if:
  - One family rejected (p < Bonferroni threshold)
  - Other family accepted (p > 0.05)

Reads:
  - runs/{run_id}/artifacts/test.npz
  - runs/{run_id}/artifacts/fit_family_a.json
  - runs/{run_id}/artifacts/fit_family_b.json

Writes: runs/{run_id}/artifacts/stratified_tests.json
"""
import argparse
import json
import os
from pathlib import Path

import numpy as np
from scipy import stats

def family_a_transform(met, n_jets, lead_jet_pt, alpha, gamma):
    """Family A: Jet-correlated MET scaling."""
    scale = 1.0 + alpha * n_jets + gamma * np.maximum(0, lead_jet_pt - 30) / 30.0
    return met * scale


def family_b_transform(met, npvs, beta0, beta1, rng):
    """Family B: Pileup-correlated MET broadening."""
    sigma = beta0 + beta1 * np.maximum(0, npvs - 20)
    sigma = np.maximum(sigma, 0.01)
    delta = rng.rayleigh(sigma)
    return met + delta


def test_stratum(met_obs, met_model_a, met_model_b, mask, bonferroni_alpha):
    """
    Test a single stratum.
    Returns: n_events, p_a, p_b, reject_a, reject_b, separating
    """
    n_events = int(np.sum(mask))
    
    if n_events < 50:
        # Too few events for reliable KS test
        return n_events, None, None, False, False, False
    
    obs_stratum = met_obs[mask]
    model_a_stratum = met_model_a[mask]
    model_b_stratum = met_model_b[mask]
    
    _, p_a = stats.ks_2samp(obs_stratum, model_a_stratum)
    _, p_b = stats.ks_2samp(obs_stratum, model_b_stratum)
    
    # Reject if p < Bonferroni threshold
    reject_a = p_a < bonferroni_alpha
    reject_b = p_b < bonferroni_alpha
    
    # Accept if p > 0.05
    accept_a = p_a > 0.05
    accept_b = p_b > 0.05
    
    # Separating: one rejected AND other accepted
    separating = (reject_a and accept_b) or (reject_b and accept_a)
    
    return n_events, float(p_a), float(p_b), reject_a, reject_b, separating


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--run-id", required=True)
    args = parser.parse_args()
    
    run_id = args.run_id
    project_root = Path.cwd()
    artifacts_dir = project_root / "runs" / run_id / "artifacts"
    
    seed = int(os.environ.get("LOCKED_EXEC_SEED", 42))
    
    # Load test data
    test_data = np.load(artifacts_dir / "test.npz")
    met_obs = test_data["met_pt"]
    n_jets = test_data["n_jets"]
    lead_jet_pt = test_data["lead_jet_pt"]
    npvs = test_data["npvs"]
    z_pt = test_data["z_pt"]
    
    n_test = len(met_obs)
    print(f"[S7] Stratified tests on {n_test} test events...")
    
    # Load fitted parameters
    with open(artifacts_dir / "fit_family_a.json") as f:
        fit_a = json.load(f)
    with open(artifacts_dir / "fit_family_b.json") as f:
        fit_b = json.load(f)
    
    alpha = fit_a["parameters"]["alpha"]["value"]
    gamma = fit_a["parameters"]["gamma"]["value"]
    beta0 = fit_b["parameters"]["beta0"]["value"]
    beta1 = fit_b["parameters"]["beta1"]["value"]
    
    # Generate model predictions
    met_model_a = family_a_transform(met_obs, n_jets, lead_jet_pt, alpha, gamma)
    rng = np.random.default_rng(seed + 3000)
    met_model_b = family_b_transform(met_obs, npvs, beta0, beta1, rng)
    
    # Define strata
    strata = []
    
    # Pileup terciles
    pileup_terciles = np.percentile(npvs, [33.33, 66.67])
    strata.append(("pileup_T1_low", f"nPV < {pileup_terciles[0]:.0f}", npvs < pileup_terciles[0]))
    strata.append(("pileup_T2_mid", f"{pileup_terciles[0]:.0f} ≤ nPV < {pileup_terciles[1]:.0f}", 
                   (npvs >= pileup_terciles[0]) & (npvs < pileup_terciles[1])))
    strata.append(("pileup_T3_high", f"nPV ≥ {pileup_terciles[1]:.0f}", npvs >= pileup_terciles[1]))
    
    # Jet multiplicity
    strata.append(("jets_0", "N_jets = 0", n_jets == 0))
    strata.append(("jets_1", "N_jets = 1", n_jets == 1))
    strata.append(("jets_2plus", "N_jets ≥ 2", n_jets >= 2))
    
    # Z pT bins
    strata.append(("zpt_0_10", "Z pT ∈ [0, 10) GeV", (z_pt >= 0) & (z_pt < 10)))
    strata.append(("zpt_10_30", "Z pT ∈ [10, 30) GeV", (z_pt >= 10) & (z_pt < 30)))
    strata.append(("zpt_30_60", "Z pT ∈ [30, 60) GeV", (z_pt >= 30) & (z_pt < 60)))
    strata.append(("zpt_60_inf", "Z pT ≥ 60 GeV", z_pt >= 60))
    
    n_strata = len(strata)
    bonferroni_alpha = 0.01 / n_strata
    
    print(f"[S7] Testing {n_strata} strata with Bonferroni α = {bonferroni_alpha:.6f}")
    
    # Run tests
    results = []
    for stratum_id, description, mask in strata:
        n_ev, p_a, p_b, rej_a, rej_b, sep = test_stratum(
            met_obs, met_model_a, met_model_b, mask, bonferroni_alpha
        )
        
        results.append({
            "stratum_id": stratum_id,
            "stratum_description": description,
            "n_events": n_ev,
            "family_a_ks_p": p_a,
            "family_b_ks_p": p_b,
            "family_a_reject": rej_a,
            "family_b_reject": rej_b,
            "separating": sep
        })
        
        sep_str = " [SEPARATING]" if sep else ""
        print(f"[S7]   {stratum_id}: n={n_ev}, p_A={p_a:.4f if p_a else 'N/A':>6}, "
              f"p_B={p_b:.4f if p_b else 'N/A':>6}{sep_str}")
    
    any_separating = any(r["separating"] for r in results)
    n_separating = sum(1 for r in results if r["separating"])
    
    print(f"[S7] Separating strata: {n_separating}/{n_strata}")
    
    # Output
    output = {
        "n_strata": n_strata,
        "bonferroni_alpha": float(bonferroni_alpha),
        "strata": results,
        "any_separating": any_separating,
        "n_separating": n_separating
    }
    
    output_file = artifacts_dir / "stratified_tests.json"
    with open(output_file, "w") as f:
        json.dump(output, f, indent=2)
    print(f"[S7] Saved {output_file}")
    
    print(f"[S7] Complete")

if __name__ == "__main__":
    main()
