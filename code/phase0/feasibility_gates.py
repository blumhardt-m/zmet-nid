#!/usr/bin/env python3
"""
ZMET-NID Phase 0 Feasibility Gates

Hard gates:
  G1: Can we open the NanoAOD file?
  G2: Are required branches present?

Soft gate:
  G3: Do we see a Z->mumu peak? (sanity check only)

This script is deliberately conservative and diagnostic.
"""

from pathlib import Path
import sys

# -------------------------------------------------------------------
# Configuration
# -------------------------------------------------------------------

FILE_URL_PATH = Path("data/manifests/phase0_file_url.txt")

REQUIRED_BRANCHES = [
    "Muon_pt", "Muon_eta", "Muon_phi", "Muon_charge", "Muon_mass",
    "MET_pt", "MET_phi",
    "Jet_pt", "Jet_eta", "Jet_phi",
    "PV_npvs",
]

# -------------------------------------------------------------------
# Gate 0: Cheap HTTP preflight (diagnostic)
# -------------------------------------------------------------------

def preflight_url(url: str) -> bool:
    if not url.startswith("http"):
        print("[G0 SKIP] Non-HTTP URL (likely local file).")
        return True

    try:
        import fsspec
        with fsspec.open(url, "rb") as f:
            chunk = f.read(1024)
        print(f"[G0 PASS] HTTP preflight ok (read {len(chunk)} bytes).")
        return True
    except Exception as e:
        print(f"[G0 FAIL] HTTP preflight failed: {e}")
        return False

# -------------------------------------------------------------------
# Gate 1: Connectivity / uproot open
# -------------------------------------------------------------------

def gate1_connectivity(url: str):
    import uproot

    try:
        f = uproot.open(url)
        tree = f["Events"]
        print(f"[G1 PASS] File opened. Events entries = {tree.num_entries:,}")
        return tree
    except Exception as e:
        print(f"[G1 FAIL] Cannot open file with uproot: {e}")
        return None

# -------------------------------------------------------------------
# Gate 2: Required branches
# -------------------------------------------------------------------

def gate2_branches(tree) -> bool:
    if tree is None:
        print("[G2 SKIP] No tree available.")
        return False

    available = set(tree.keys())
    missing = [b for b in REQUIRED_BRANCHES if b not in available]

    if not missing:
        print(f"[G2 PASS] All {len(REQUIRED_BRANCHES)} required branches present.")
        return True

    print(f"[G2 FAIL] Missing branches: {missing}")
    for m in missing:
        prefix = m.split("_")[0]
        candidates = sorted(k for k in available if k.startswith(prefix))[:10]
        if candidates:
            print(f"  Candidates for {m}: {candidates}")
    return False

# -------------------------------------------------------------------
# Gate 3: Z peak sanity check (soft)
# -------------------------------------------------------------------

def gate3_sanity(tree, n_events=100_000):
    import awkward as ak
    import numpy as np

    if tree is None:
        print("[G3 SKIP] No tree available.")
        return

    try:
        arr = tree.arrays(REQUIRED_BRANCHES, entry_stop=n_events)
    except Exception as e:
        print(f"[G3 FAIL] Cannot read arrays: {e}")
        return

    # Muon selection
    mu_pt  = arr["Muon_pt"]
    mu_eta = arr["Muon_eta"]
    mu_phi = arr["Muon_phi"]
    mu_q   = arr["Muon_charge"]
    mu_m   = arr["Muon_mass"]

    mask = (mu_pt > 15) & (np.abs(mu_eta) < 2.4)
    mu_pt, mu_eta, mu_phi, mu_q, mu_m = (
        mu_pt[mask], mu_eta[mask], mu_phi[mask], mu_q[mask], mu_m[mask]
    )

    n_mu = ak.num(mu_pt)
    qsum = ak.sum(mu_q, axis=1)
    dimu = (n_mu == 2) & (qsum == 0)

    if ak.sum(dimu) < 100:
        print(f"[G3 WARN] Only {ak.sum(dimu)} dimuon candidates. Check selection.")
        return

    mu_pt  = mu_pt[dimu]
    mu_eta = mu_eta[dimu]
    mu_phi = mu_phi[dimu]
    mu_m   = mu_m[dimu]

    # Invariant mass
    px1 = mu_pt[:,0] * np.cos(mu_phi[:,0])
    py1 = mu_pt[:,0] * np.sin(mu_phi[:,0])
    pz1 = mu_pt[:,0] * np.sinh(mu_eta[:,0])
    e1  = np.sqrt(px1**2 + py1**2 + pz1**2 + mu_m[:,0]**2)

    px2 = mu_pt[:,1] * np.cos(mu_phi[:,1])
    py2 = mu_pt[:,1] * np.sin(mu_phi[:,1])
    pz2 = mu_pt[:,1] * np.sinh(mu_eta[:,1])
    e2  = np.sqrt(px2**2 + py2**2 + pz2**2 + mu_m[:,1]**2)

    m2 = (e1 + e2)**2 - (px1 + px2)**2 - (py1 + py2)**2 - (pz1 + pz2)**2
    mll = np.sqrt(np.maximum(ak.to_numpy(m2), 0))

    zwin = (mll > 80) & (mll < 100)
    nz   = zwin.sum()
    med  = np.median(mll[zwin]) if nz > 0 else None

    if nz > 1000 and med is not None and 88 < med < 94:
        print(f"[G3 PASS] Z peak visible: N(80–100)={nz}, median={med:.1f} GeV")
    else:
        print(f"[G3 WARN] Z peak unclear: N(80–100)={nz}, median={med}")
        print("          Sanity warning only; not a hard failure.")

# -------------------------------------------------------------------
# Main
# -------------------------------------------------------------------

def main():
    print("=" * 60)
    print("ZMET-NID Phase 0 Feasibility Gates")
    print("=" * 60)

    if not FILE_URL_PATH.exists():
        print(f"[FATAL] Missing {FILE_URL_PATH}. Generate manifest first.")
        sys.exit(1)

    url = FILE_URL_PATH.read_text().strip()
    print(f"[INFO] FILE_URL = {url}")

    preflight_url(url)

    tree = gate1_connectivity(url)
    branches_ok = gate2_branches(tree)
    gate3_sanity(tree)

    print("=" * 60)
    print("SUMMARY")
    print("=" * 60)

    if tree is None:
        print("RESULT: FAIL (G1 connectivity). Fix data access.")
        sys.exit(1)

    if not branches_ok:
        print("RESULT: FAIL (G2 branches). Assess MiniAOD or branch mapping.")
        sys.exit(1)

    print("RESULT: PASS. Proceed to Phase 0b (plots + thresholds).")
    sys.exit(0)

if __name__ == "__main__":
    main()
