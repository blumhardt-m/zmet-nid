#!/usr/bin/env python3
"""
ZMET-NID Phase 1 — S2: Train/Test Split
========================================
Protocol v0.4 §6: Split by lumisection parity
  - Even lumisections → train
  - Odd lumisections → test

Reads: runs/{run_id}/artifacts/selected_events.npz
Writes:
  - runs/{run_id}/artifacts/train.npz
  - runs/{run_id}/artifacts/test.npz
  - runs/{run_id}/meta/split_metrics.json
"""
import argparse
import json
from pathlib import Path

import numpy as np

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--run-id", required=True)
    args = parser.parse_args()
    
    run_id = args.run_id
    project_root = Path.cwd()
    
    artifacts_dir = project_root / "runs" / run_id / "artifacts"
    meta_dir = project_root / "runs" / run_id / "meta"
    
    # Load selected events
    input_file = artifacts_dir / "selected_events.npz"
    data = np.load(input_file)
    
    lumi_block = data["lumi_block"]
    n_total = len(lumi_block)
    
    # Split by parity
    train_mask = (lumi_block % 2 == 0)
    test_mask = (lumi_block % 2 == 1)
    
    n_train = int(np.sum(train_mask))
    n_test = int(np.sum(test_mask))
    
    print(f"[S2] Total: {n_total}, Train: {n_train}, Test: {n_test}")
    
    # Save train split
    train_file = artifacts_dir / "train.npz"
    np.savez_compressed(
        train_file,
        **{key: data[key][train_mask] for key in data.files}
    )
    print(f"[S2] Saved {train_file}")
    
    # Save test split
    test_file = artifacts_dir / "test.npz"
    np.savez_compressed(
        test_file,
        **{key: data[key][test_mask] for key in data.files}
    )
    print(f"[S2] Saved {test_file}")
    
    # Write split metrics
    metrics = {
        "n_total": n_total,
        "n_train": n_train,
        "n_test": n_test,
        "train_fraction": float(n_train / n_total) if n_total > 0 else 0.0,
        "test_fraction": float(n_test / n_total) if n_total > 0 else 0.0,
        "split_method": "lumisection_parity_even_train"
    }
    
    metrics_file = meta_dir / "split_metrics.json"
    with open(metrics_file, "w") as f:
        json.dump(metrics, f, indent=2)
    print(f"[S2] Saved {metrics_file}")
    
    print(f"[S2] Complete")

if __name__ == "__main__":
    main()
