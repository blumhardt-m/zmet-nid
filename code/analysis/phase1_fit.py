#!/usr/bin/env python3
"""
ZMET-NID Phase 1 — S4/S5: Nuisance Family Fitting
==================================================
Protocol v0.4 §5:

Family A (Jet-correlated MET scaling):
  MET' = MET × (1 + α × N_jets + γ × max(0, pT_lead - 30)/30)
  Parameters: θ_A = {α, γ}
  Bounds: α ∈ [0, 0.3], γ ∈ [0, 0.3]

Family B (Pileup-correlated MET broadening):
  MET' = MET + |Δ|, where |Δ| ~ Rayleigh(σ)
  σ = β₀ + β₁ × max(0, nPV - 20)
  Parameters: θ_B = {β₀, β₁}
  Bounds: β₀ ∈ [0, 20], β₁ ∈ [0, 2]

Fit metric (protocol v0.4 §7):
  D = 0.5 × W₁ + 0.5 × KS

Fit on training split only.

Reads: runs/{run_id}/artifacts/train.npz
Writes: runs/{run_id}/artifacts/fit_family_{a,b}.json

Gates: G4 (Family A), G5 (Family B)
"""
import argparse
import json
import os
from pathlib import Path

import numpy as np
from scipy import optimize, stats

def wasserstein_1d(u, v):
    """1D Wasserstein distance."""
    u_sorted = np.sort(u)
    v_sorted = np.sort(v)
    all_values = np.sort(np.concatenate([u_sorted, v_sorted]))
    u_cdf = np.searchsorted(u_sorted, all_values, side='right') / len(u)
    v_cdf = np.searchsorted(v_sorted, all_values, side='right') / len(v)
    deltas = np.diff(all_values, prepend=all_values[0])
    return np.sum(np.abs(u_cdf - v_cdf) * deltas)


def fit_metric(observed, model):
    """Combined fit metric: D = 0.5 × W₁ + 0.5 × KS."""
    w1 = wasserstein_1d(observed, model)
    ks_stat, _ = stats.ks_2samp(observed, model)
    return 0.5 * w1 + 0.5 * ks_stat, w1, ks_stat


def family_a_transform(met, n_jets, lead_jet_pt, alpha, gamma):
    """
    Family A: Jet-correlated MET scaling.
    MET' = MET × (1 + α × N_jets + γ × max(0, pT_lead - 30)/30)
    """
    jet_term = alpha * n_jets
    lead_term = gamma * np.maximum(0, lead_jet_pt - 30) / 30.0
    scale = 1.0 + jet_term + lead_term
    return met * scale


def family_b_transform(met, npvs, beta0, beta1, rng):
    """
    Family B: Pileup-correlated MET broadening.
    MET' = MET + |Δ|, |Δ| ~ Rayleigh(σ)
    σ = β₀ + β₁ × max(0, nPV - 20)
    """
    sigma = beta0 + beta1 * np.maximum(0, npvs - 20)
    # Ensure sigma is positive
    sigma = np.maximum(sigma, 0.01)
    delta = rng.rayleigh(sigma)
    return met + delta


def fit_family_a(met_obs, n_jets, lead_jet_pt):
    """Fit Family A parameters by minimizing D."""
    
    def objective(params):
        alpha, gamma = params
        met_model = family_a_transform(met_obs, n_jets, lead_jet_pt, alpha, gamma)
        D, _, _ = fit_metric(met_obs, met_model)
        return D
    
    # Initial guess
    x0 = [0.05, 0.05]
    bounds = [(0, 0.3), (0, 0.3)]
    
    # Minimize
    result = optimize.minimize(
        objective, x0, method='L-BFGS-B', bounds=bounds,
        options={'maxiter': 500}
    )
    
    return result


def fit_family_b(met_obs, npvs, rng, seed_offset=0):
    """Fit Family B parameters by minimizing D."""
    
    # Use deterministic seed for Rayleigh draws during optimization
    base_seed = int(os.environ.get("LOCKED_EXEC_SEED", 42))
    
    def objective(params):
        beta0, beta1 = params
        # Fixed random state for reproducibility within optimization
        rng_fixed = np.random.default_rng(base_seed + seed_offset)
        met_model = family_b_transform(met_obs, npvs, beta0, beta1, rng_fixed)
        D, _, _ = fit_metric(met_obs, met_model)
        return D
    
    # Initial guess
    x0 = [5.0, 0.5]
    bounds = [(0, 20), (0, 2)]
    
    # Minimize
    result = optimize.minimize(
        objective, x0, method='L-BFGS-B', bounds=bounds,
        options={'maxiter': 500}
    )
    
    return result


def bootstrap_ci(met_obs, n_jets, lead_jet_pt, npvs, family, best_params, rng, n_boot=50):
    """Bootstrap confidence intervals for parameters."""
    n = len(met_obs)
    base_seed = int(os.environ.get("LOCKED_EXEC_SEED", 42))
    
    param_samples = []
    
    for i in range(n_boot):
        # Resample
        idx = rng.choice(n, size=n, replace=True)
        met_boot = met_obs[idx]
        
        if family == "A":
            n_jets_boot = n_jets[idx]
            lead_boot = lead_jet_pt[idx]
            
            def obj(params):
                met_model = family_a_transform(met_boot, n_jets_boot, lead_boot, params[0], params[1])
                D, _, _ = fit_metric(met_boot, met_model)
                return D
            
            bounds = [(0, 0.3), (0, 0.3)]
            res = optimize.minimize(obj, best_params, method='L-BFGS-B', bounds=bounds)
            param_samples.append(res.x)
            
        else:  # Family B
            npvs_boot = npvs[idx]
            
            def obj(params):
                rng_b = np.random.default_rng(base_seed + 1000 + i)
                met_model = family_b_transform(met_boot, npvs_boot, params[0], params[1], rng_b)
                D, _, _ = fit_metric(met_boot, met_model)
                return D
            
            bounds = [(0, 20), (0, 2)]
            res = optimize.minimize(obj, best_params, method='L-BFGS-B', bounds=bounds)
            param_samples.append(res.x)
    
    param_samples = np.array(param_samples)
    
    # 95% CI
    ci_lo = np.percentile(param_samples, 2.5, axis=0)
    ci_hi = np.percentile(param_samples, 97.5, axis=0)
    
    return ci_lo, ci_hi


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--run-id", required=True)
    parser.add_argument("--family", required=True, choices=["A", "B"])
    args = parser.parse_args()
    
    run_id = args.run_id
    family = args.family
    
    project_root = Path.cwd()
    artifacts_dir = project_root / "runs" / run_id / "artifacts"
    
    seed = int(os.environ.get("LOCKED_EXEC_SEED", 42))
    rng = np.random.default_rng(seed)
    
    # Load training data
    train_data = np.load(artifacts_dir / "train.npz")
    met_obs = train_data["met_pt"]
    n_jets = train_data["n_jets"]
    lead_jet_pt = train_data["lead_jet_pt"]
    npvs = train_data["npvs"]
    
    print(f"[S{4 if family == 'A' else 5}] Fitting Family {family} on {len(met_obs)} training events...")
    
    if family == "A":
        # Fit Family A
        result = fit_family_a(met_obs, n_jets, lead_jet_pt)
        alpha_opt, gamma_opt = result.x
        
        print(f"[S4] Optimal: α={alpha_opt:.6f}, γ={gamma_opt:.6f}")
        
        # Compute final metrics
        met_model = family_a_transform(met_obs, n_jets, lead_jet_pt, alpha_opt, gamma_opt)
        D, w1, ks_stat = fit_metric(met_obs, met_model)
        
        # Bootstrap CI
        ci_lo, ci_hi = bootstrap_ci(met_obs, n_jets, lead_jet_pt, npvs, "A", 
                                     [alpha_opt, gamma_opt], rng)
        
        # Check significance (CI excludes zero)
        alpha_sig = (ci_lo[0] > 1e-6) or (ci_hi[0] < -1e-6)
        gamma_sig = (ci_lo[1] > 1e-6) or (ci_hi[1] < -1e-6)
        any_sig = alpha_sig or gamma_sig
        
        output = {
            "family": "A",
            "converged": bool(result.success),
            "any_param_significant": any_sig,
            "parameters": {
                "alpha": {
                    "value": float(alpha_opt),
                    "ci_lo": float(ci_lo[0]),
                    "ci_hi": float(ci_hi[0])
                },
                "gamma": {
                    "value": float(gamma_opt),
                    "ci_lo": float(ci_lo[1]),
                    "ci_hi": float(ci_hi[1])
                }
            },
            "fit_metric_D": float(D),
            "w1": float(w1),
            "ks_statistic": float(ks_stat),
            "n_function_evals": int(result.nfev)
        }
        
        output_file = artifacts_dir / "fit_family_a.json"
        
    else:  # Family B
        # Fit Family B
        result = fit_family_b(met_obs, npvs, rng)
        beta0_opt, beta1_opt = result.x
        
        print(f"[S5] Optimal: β₀={beta0_opt:.6f}, β₁={beta1_opt:.6f}")
        
        # Compute final metrics
        rng_final = np.random.default_rng(seed + 999)
        met_model = family_b_transform(met_obs, npvs, beta0_opt, beta1_opt, rng_final)
        D, w1, ks_stat = fit_metric(met_obs, met_model)
        
        # Bootstrap CI
        ci_lo, ci_hi = bootstrap_ci(met_obs, n_jets, lead_jet_pt, npvs, "B",
                                     [beta0_opt, beta1_opt], rng)
        
        # Check significance
        beta0_sig = (ci_lo[0] > 1e-6) or (ci_hi[0] < -1e-6)
        beta1_sig = (ci_lo[1] > 1e-6) or (ci_hi[1] < -1e-6)
        any_sig = beta0_sig or beta1_sig
        
        output = {
            "family": "B",
            "converged": bool(result.success),
            "any_param_significant": any_sig,
            "parameters": {
                "beta0": {
                    "value": float(beta0_opt),
                    "ci_lo": float(ci_lo[0]),
                    "ci_hi": float(ci_hi[0])
                },
                "beta1": {
                    "value": float(beta1_opt),
                    "ci_lo": float(ci_lo[1]),
                    "ci_hi": float(ci_hi[1])
                }
            },
            "fit_metric_D": float(D),
            "w1": float(w1),
            "ks_statistic": float(ks_stat),
            "n_function_evals": int(result.nfev)
        }
        
        output_file = artifacts_dir / "fit_family_b.json"
    
    with open(output_file, "w") as f:
        json.dump(output, f, indent=2)
    print(f"[S{4 if family == 'A' else 5}] Saved {output_file}")
    
    print(f"[S{4 if family == 'A' else 5}] D={D:.6f}, converged={result.success}, any_sig={any_sig}")
    print(f"[S{4 if family == 'A' else 5}] Complete")

if __name__ == "__main__":
    main()
