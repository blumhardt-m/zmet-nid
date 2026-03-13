#!/usr/bin/env python3
"""
ZMET-NID Phase 1 — S6: Global Tests
====================================
Evaluate both fitted families on held-out test split.

Protocol v0.4 §9 — Global pass criterion:
  KS p-value > 0.05 AND W₁ < τ_W

Reads:
  - runs/{run_id}/artifacts/test.npz
  - runs/{run_id}/artifacts/thresholds.json
  - runs/{run_id}/artifacts/fit_family_a.json
  - runs/{run_id}/artifacts/fit_family_b.json

Writes: runs/{run_id}/artifacts/global_tests.json

Gate: G6 (neither family passes → outcome 4)
"""
import argparse
import json
import os
from pathlib import Path

import numpy as np
from scipy import stats

def wasserstein_1d(u, v):
    """1D Wasserstein distance."""
    u_sorted = np.sort(u)
    v_sorted = np.sort(v)
    all_values = np.sort(np.concatenate([u_sorted, v_sorted]))
    u_cdf = np.searchsorted(u_sorted, all_values, side='right') / len(u)
    v_cdf = np.searchsorted(v_sorted, all_values, side='right') / len(v)
    deltas = np.diff(all_values, prepend=all_values[0])
    return np.sum(np.abs(u_cdf - v_cdf) * deltas)


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


def evaluate_family(met_obs, met_model, tau_w1):
    """
    Evaluate a family's model predictions.
    Returns: ks_p, w1, D, global_pass
    """
    ks_stat, ks_p = stats.ks_2samp(met_obs, met_model)
    w1 = wasserstein_1d(met_obs, met_model)
    D = 0.5 * w1 + 0.5 * ks_stat
    
    # Global pass: KS p > 0.05 AND W1 < tau_W
    global_pass = (ks_p > 0.05) and (w1 < tau_w1)
    
    return ks_p, w1, D, global_pass


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
    
    print(f"[S6] Evaluating on {len(met_obs)} test events...")
    
    # Load thresholds
    with open(artifacts_dir / "thresholds.json") as f:
        thresholds = json.load(f)
    tau_w1 = thresholds["tau_w1"]
    tau_ks = thresholds["tau_ks"]
    
    print(f"[S6] Thresholds: τ_W1={tau_w1:.6f}, τ_KS={tau_ks:.6f}")
    
    # Load fitted parameters
    with open(artifacts_dir / "fit_family_a.json") as f:
        fit_a = json.load(f)
    with open(artifacts_dir / "fit_family_b.json") as f:
        fit_b = json.load(f)
    
    # Family A evaluation
    alpha = fit_a["parameters"]["alpha"]["value"]
    gamma = fit_a["parameters"]["gamma"]["value"]
    met_model_a = family_a_transform(met_obs, n_jets, lead_jet_pt, alpha, gamma)
    ks_p_a, w1_a, D_a, pass_a = evaluate_family(met_obs, met_model_a, tau_w1)
    
    print(f"[S6] Family A: KS_p={ks_p_a:.4f}, W1={w1_a:.4f}, D={D_a:.4f}, pass={pass_a}")
    
    # Family B evaluation
    beta0 = fit_b["parameters"]["beta0"]["value"]
    beta1 = fit_b["parameters"]["beta1"]["value"]
    rng = np.random.default_rng(seed + 2000)  # Different seed for test evaluation
    met_model_b = family_b_transform(met_obs, npvs, beta0, beta1, rng)
    ks_p_b, w1_b, D_b, pass_b = evaluate_family(met_obs, met_model_b, tau_w1)
    
    print(f"[S6] Family B: KS_p={ks_p_b:.4f}, W1={w1_b:.4f}, D={D_b:.4f}, pass={pass_b}")
    
    # Output
    output = {
        "family_a": {
            "ks_p": float(ks_p_a),
            "w1": float(w1_a),
            "D": float(D_a),
            "global_pass": bool(pass_a)
        },
        "family_b": {
            "ks_p": float(ks_p_b),
            "w1": float(w1_b),
            "D": float(D_b),
            "global_pass": bool(pass_b)
        },
        "tau_w1_used": float(tau_w1),
        "tau_ks_used": float(tau_ks),
        "D_difference": float(abs(D_a - D_b))
    }
    
    output_file = artifacts_dir / "global_tests.json"
    with open(output_file, "w") as f:
        json.dump(output, f, indent=2)
    print(f"[S6] Saved {output_file}")
    
    print(f"[S6] |D_A - D_B| = {output['D_difference']:.6f}")
    print(f"[S6] Complete")

if __name__ == "__main__":
    main()
