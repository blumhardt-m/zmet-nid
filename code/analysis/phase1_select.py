#!/usr/bin/env python3
"""
ZMET-NID Phase 1 — S1: Event Selection
=======================================
Applies frozen selection cuts from protocol v0.4 §2:
  - 2 opposite-sign muons
  - pT(μ1) ≥ 25 GeV, pT(μ2) ≥ 15 GeV
  - |η(μ)| ≤ 2.4
  - Tight ID, PFIso < 0.15
  - M(μμ) ∈ [80, 100] GeV

Reads: data/raw/test.root
Writes:
  - runs/{run_id}/artifacts/selected_events.npz
  - runs/{run_id}/meta/selection_metrics.json

Gates: G1 (n_selected < 50k), G2 (Z peak outside [88,94] GeV)
"""
import argparse
import json
import os
import sys
from pathlib import Path

import numpy as np

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--run-id", required=True)
    args = parser.parse_args()
    
    run_id = args.run_id
    project_root = Path.cwd()
    
    # Output paths
    artifacts_dir = project_root / "runs" / run_id / "artifacts"
    meta_dir = project_root / "runs" / run_id / "meta"
    artifacts_dir.mkdir(parents=True, exist_ok=True)
    meta_dir.mkdir(parents=True, exist_ok=True)
    
    # Input
    input_file = project_root / "data" / "raw" / "test.root"
    if not input_file.exists():
        print(f"ERROR: Input file not found: {input_file}", file=sys.stderr)
        sys.exit(1)
    
    # Load NanoAOD
    try:
        import uproot
        import awkward as ak
    except ImportError as e:
        print(f"ERROR: Missing dependency: {e}", file=sys.stderr)
        sys.exit(1)
    
    print(f"[S1] Loading {input_file}...")
    tree = uproot.open(f"{input_file}:Events")
    
    # Read required branches
    branches = [
        "Muon_pt", "Muon_eta", "Muon_phi", "Muon_mass", "Muon_charge",
        "Muon_tightId", "Muon_pfRelIso04_all",
        "MET_pt", "MET_phi",
        "Jet_pt", "Jet_eta",
        "PV_npvs",
        "luminosityBlock"
    ]
    
    events = tree.arrays(branches, library="ak")
    n_raw = len(events)
    print(f"[S1] Loaded {n_raw} events")
    
    cut_flow = {"raw": n_raw}
    
    # Cut 1: Exactly 2 muons
    n_muons = ak.num(events.Muon_pt)
    mask = (n_muons == 2)
    events = events[mask]
    cut_flow["2_muons"] = len(events)
    print(f"[S1] After 2 muons: {len(events)}")
    
    # Cut 2: Opposite sign
    charge_sum = ak.sum(events.Muon_charge, axis=1)
    mask = (charge_sum == 0)
    events = events[mask]
    cut_flow["opposite_sign"] = len(events)
    print(f"[S1] After opposite sign: {len(events)}")
    
    # Cut 3: pT cuts (leading ≥ 25, subleading ≥ 15)
    pt_sorted = ak.sort(events.Muon_pt, axis=1, ascending=False)
    pt_lead = pt_sorted[:, 0]
    pt_sublead = pt_sorted[:, 1]
    mask = (pt_lead >= 25) & (pt_sublead >= 15)
    events = events[mask]
    cut_flow["pt_cuts"] = len(events)
    print(f"[S1] After pT cuts: {len(events)}")
    
    # Cut 4: |η| ≤ 2.4
    eta_max = ak.max(np.abs(events.Muon_eta), axis=1)
    mask = (eta_max <= 2.4)
    events = events[mask]
    cut_flow["eta_cut"] = len(events)
    print(f"[S1] After η cut: {len(events)}")
    
    # Cut 5: Tight ID (both muons)
    tight_both = ak.all(events.Muon_tightId, axis=1)
    mask = tight_both
    events = events[mask]
    cut_flow["tight_id"] = len(events)
    print(f"[S1] After tight ID: {len(events)}")
    
    # Cut 6: Isolation (PFIso < 0.15 for both)
    iso_pass = ak.all(events.Muon_pfRelIso04_all < 0.15, axis=1)
    mask = iso_pass
    events = events[mask]
    cut_flow["isolation"] = len(events)
    print(f"[S1] After isolation: {len(events)}")
    
    # Cut 7: M(μμ) ∈ [80, 100] GeV
    # Compute invariant mass
    mu1_pt = events.Muon_pt[:, 0]
    mu2_pt = events.Muon_pt[:, 1]
    mu1_eta = events.Muon_eta[:, 0]
    mu2_eta = events.Muon_eta[:, 1]
    mu1_phi = events.Muon_phi[:, 0]
    mu2_phi = events.Muon_phi[:, 1]
    mu1_mass = events.Muon_mass[:, 0]
    mu2_mass = events.Muon_mass[:, 1]
    
    # 4-vector components
    mu1_px = mu1_pt * np.cos(mu1_phi)
    mu1_py = mu1_pt * np.sin(mu1_phi)
    mu1_pz = mu1_pt * np.sinh(mu1_eta)
    mu1_E = np.sqrt(mu1_px**2 + mu1_py**2 + mu1_pz**2 + mu1_mass**2)
    
    mu2_px = mu2_pt * np.cos(mu2_phi)
    mu2_py = mu2_pt * np.sin(mu2_phi)
    mu2_pz = mu2_pt * np.sinh(mu2_eta)
    mu2_E = np.sqrt(mu2_px**2 + mu2_py**2 + mu2_pz**2 + mu2_mass**2)
    
    # Dimuon invariant mass
    mumu_px = mu1_px + mu2_px
    mumu_py = mu1_py + mu2_py
    mumu_pz = mu1_pz + mu2_pz
    mumu_E = mu1_E + mu2_E
    mumu_mass = np.sqrt(mumu_E**2 - mumu_px**2 - mumu_py**2 - mumu_pz**2)
    
    mask = (mumu_mass >= 80) & (mumu_mass <= 100)
    events = events[mask]
    mumu_mass = mumu_mass[mask]
    cut_flow["mass_window"] = len(events)
    print(f"[S1] After mass window: {len(events)}")
    
    n_selected = len(events)
    
    # Compute Z peak (mode of mass distribution)
    mass_np = ak.to_numpy(mumu_mass)
    hist, bin_edges = np.histogram(mass_np, bins=100, range=(80, 100))
    peak_bin = np.argmax(hist)
    peak_mass_gev = (bin_edges[peak_bin] + bin_edges[peak_bin + 1]) / 2
    print(f"[S1] Z peak at {peak_mass_gev:.2f} GeV")
    
    # Extract arrays for saving
    # Recompute for selected events
    mu1_pt_sel = ak.to_numpy(events.Muon_pt[:, 0])
    mu2_pt_sel = ak.to_numpy(events.Muon_pt[:, 1])
    mu1_eta_sel = ak.to_numpy(events.Muon_eta[:, 0])
    mu2_eta_sel = ak.to_numpy(events.Muon_eta[:, 1])
    
    met_pt = ak.to_numpy(events.MET_pt)
    met_phi = ak.to_numpy(events.MET_phi)
    
    # Jet counting (pT > 30 GeV, |η| < 2.4)
    jet_pt = events.Jet_pt
    jet_eta = events.Jet_eta
    good_jets = (jet_pt > 30) & (np.abs(jet_eta) < 2.4)
    n_jets = ak.sum(good_jets, axis=1)
    n_jets_np = ak.to_numpy(n_jets)
    
    # Leading jet pT (0 if no jets)
    lead_jet_pt = ak.fill_none(ak.firsts(jet_pt[good_jets]), 0)
    lead_jet_pt_np = ak.to_numpy(lead_jet_pt)
    
    npvs = ak.to_numpy(events.PV_npvs)
    lumi_block = ak.to_numpy(events.luminosityBlock)
    
    # Z pT
    z_pt = np.sqrt((mu1_pt_sel * np.cos(ak.to_numpy(events.Muon_phi[:, 0])) + 
                    mu2_pt_sel * np.cos(ak.to_numpy(events.Muon_phi[:, 1])))**2 +
                   (mu1_pt_sel * np.sin(ak.to_numpy(events.Muon_phi[:, 0])) + 
                    mu2_pt_sel * np.sin(ak.to_numpy(events.Muon_phi[:, 1])))**2)
    
    # Save selected events
    output_file = artifacts_dir / "selected_events.npz"
    np.savez_compressed(
        output_file,
        met_pt=met_pt,
        met_phi=met_phi,
        n_jets=n_jets_np,
        lead_jet_pt=lead_jet_pt_np,
        npvs=npvs,
        mumu_mass=mass_np,
        z_pt=z_pt,
        lumi_block=lumi_block,
        mu1_pt=mu1_pt_sel,
        mu2_pt=mu2_pt_sel,
        mu1_eta=mu1_eta_sel,
        mu2_eta=mu2_eta_sel
    )
    print(f"[S1] Saved {output_file}")
    
    # Write selection metrics
    metrics = {
        "n_raw": int(n_raw),
        "n_selected": int(n_selected),
        "peak_mass_gev": float(peak_mass_gev),
        "selection_efficiency": float(n_selected / n_raw) if n_raw > 0 else 0.0,
        "cut_flow": {k: int(v) for k, v in cut_flow.items()}
    }
    
    metrics_file = meta_dir / "selection_metrics.json"
    with open(metrics_file, "w") as f:
        json.dump(metrics, f, indent=2)
    print(f"[S1] Saved {metrics_file}")
    
    print(f"[S1] Complete: {n_selected} events selected ({100*n_selected/n_raw:.1f}%)")
    
if __name__ == "__main__":
    main()
