#!/usr/bin/env python3
"""
ZMET-NID Phase 1 — S3: Bootstrap Thresholds
============================================
Protocol v0.4 §8:
  1. Draw 100 bootstrap resamples of observed MET (pilot data, first 100k)
  2. Compute pairwise W₁ and KS distances
  3. Set τ_W = median(W₁), τ_KS = median(KS)

Thresholds are frozen before test-split unblinding.

Reads: runs/{run_id}/artifacts/train.npz
Writes: runs/{run_id}/artifacts/thresholds.json

Gate: G3 (thresholds non-positive)
"""
import argparse
import json
import os
from pathlib import Path

import numpy as np
from scipy import stats

def wasserstein_1d(u, v):
    """1D Wasserstein distance between empirical distributions."""
    u_sorted = np.sort(u)
    v_sorted = np.sort(v)
    
    # Combine and sort all values
    all_values = np.concatenate([u_sorted, v_sorted])
    all_values = np.sort(all_values)
    
    # Compute CDFs at each point
    u_cdf = np.searchsorted(u_sorted, all_values, side='right') / len(u)
    v_cdf = np.searchsorted(v_sorted, all_values, side='right') / len(v)
    
    # Integrate absolute difference
    deltas = np.diff(all_values, prepend=all_values[0])
    return np.sum(np.abs(u_cdf - v_cdf) * deltas)

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--run-id", required=True)
    parser.add_argument("--n-bootstrap", type=int, default=100)
    parser.add_argument("--pilot-n", type=int, default=100000)
    args = parser.parse_args()
    
    run_id = args.run_id
    n_bootstrap = args.n_bootstrap
    pilot_n = args.pilot_n
    
    project_root = Path.cwd()
    artifacts_dir = project_root / "runs" / run_id / "artifacts"
    
    # Reproducibility: use LOCKED_EXEC_SEED
    seed = int(os.environ.get("LOCKED_EXEC_SEED", 42))
    rng = np.random.default_rng(seed)
    
    # Load training data
    train_data = np.load(artifacts_dir / "train.npz")
    met_pt = train_data["met_pt"]
    
    # Use pilot subset (first pilot_n events)
    n_available = len(met_pt)
    actual_pilot_n = min(pilot_n, n_available)
    met_pilot = met_pt[:actual_pilot_n]
    
    print(f"[S3] Computing thresholds on {actual_pilot_n} events with {n_bootstrap} bootstraps...")
    
    # Generate bootstrap resamples
    resamples = []
    for _ in range(n_bootstrap):
        idx = rng.choice(actual_pilot_n, size=actual_pilot_n, replace=True)
        resamples.append(met_pilot[idx])
    
    # Compute pairwise distances
    w1_distances = []
    ks_distances = []
    
    n_pairs = 0
    for i in range(n_bootstrap):
        for j in range(i + 1, n_bootstrap):
            w1 = wasserstein_1d(resamples[i], resamples[j])
            ks_stat, _ = stats.ks_2samp(resamples[i], resamples[j])
            w1_distances.append(w1)
            ks_distances.append(ks_stat)
            n_pairs += 1
    
    w1_arr = np.array(w1_distances)
    ks_arr = np.array(ks_distances)
    
    print(f"[S3] Computed {n_pairs} pairwise distances")
    
    # Set thresholds as medians
    tau_w1 = float(np.median(w1_arr))
    tau_ks = float(np.median(ks_arr))
    
    print(f"[S3] τ_W1 = {tau_w1:.6f}")
    print(f"[S3] τ_KS = {tau_ks:.6f}")
    
    # Output
    thresholds = {
        "tau_w1": tau_w1,
        "tau_ks": tau_ks,
        "n_bootstrap": n_bootstrap,
        "pilot_n": actual_pilot_n,
        "w1_distribution": {
            "mean": float(np.mean(w1_arr)),
            "std": float(np.std(w1_arr)),
            "p05": float(np.percentile(w1_arr, 5)),
            "p95": float(np.percentile(w1_arr, 95))
        },
        "ks_distribution": {
            "mean": float(np.mean(ks_arr)),
            "std": float(np.std(ks_arr)),
            "p05": float(np.percentile(ks_arr, 5)),
            "p95": float(np.percentile(ks_arr, 95))
        }
    }
    
    output_file = artifacts_dir / "thresholds.json"
    with open(output_file, "w") as f:
        json.dump(thresholds, f, indent=2)
    print(f"[S3] Saved {output_file}")
    
    print(f"[S3] Complete")

if __name__ == "__main__":
    main()
