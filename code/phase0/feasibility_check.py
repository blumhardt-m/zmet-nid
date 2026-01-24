#!/usr/bin/env python3
"""
ZMET-NID Phase 0 Feasibility Check
===================================

Gates:
  G1: File connectivity
  G2: Required branch availability  
  G3: Z peak sanity check

Usage:
  python feasibility_check.py [--local FILE]

If --local not provided, attempts XRootD access to CMS Open Data.
"""

import sys
import argparse
from pathlib import Path

# Configuration
REQUIRED_BRANCHES = [
    # Muons
    'Muon_pt', 'Muon_eta', 'Muon_phi', 'Muon_charge', 'Muon_mass',
    'Muon_tightId', 'Muon_pfRelIso04_all',
    'nMuon',
    # MET
    'MET_pt', 'MET_phi',
    # Jets
    'Jet_pt', 'Jet_eta', 'Jet_phi',
    'nJet',
    # Pileup
    'PV_npvs',
    # Event info
    'run', 'luminosityBlock', 'event',
]

# Alternative branch names to check
BRANCH_ALTERNATIVES = {
    'Muon_tightId': ['Muon_mediumId', 'Muon_looseId'],
    'Muon_pfRelIso04_all': ['Muon_pfRelIso03_all', 'Muon_miniPFRelIso_all'],
}

# Test file (DoubleMuon 2016 NanoAOD)
TEST_FILE_XROOTD = "root://eospublic.cern.ch//eos/opendata/cms/Run2016H/DoubleMuon/NANOAOD/UL2016_MiniAODv2_NanoAODv9-v1/2510000/127C2975-1B1C-A046-AABF-62B77E757A86.root"


def check_uproot_available():
    """Check if uproot is installed."""
    try:
        import uproot
        return True, uproot.__version__
    except ImportError:
        return False, None


def gate_g1_connectivity(filepath: str) -> dict:
    """G1: Test file connectivity."""
    import uproot
    
    result = {
        'gate': 'G1',
        'name': 'File Connectivity',
        'status': 'UNKNOWN',
        'details': {}
    }
    
    try:
        f = uproot.open(filepath)
        result['status'] = 'PASS'
        result['details']['keys'] = list(f.keys())[:10]
        result['details']['file_opened'] = True
        f.close()
    except Exception as e:
        result['status'] = 'FAIL'
        result['details']['error'] = str(e)
        
    return result


def gate_g2_branches(filepath: str) -> dict:
    """G2: Check required branch availability."""
    import uproot
    
    result = {
        'gate': 'G2',
        'name': 'Branch Availability',
        'status': 'UNKNOWN',
        'details': {
            'found': [],
            'missing': [],
            'alternatives_used': {}
        }
    }
    
    try:
        f = uproot.open(filepath)
        events = f['Events']
        available_branches = set(events.keys())
        
        for branch in REQUIRED_BRANCHES:
            if branch in available_branches:
                result['details']['found'].append(branch)
            else:
                alts = BRANCH_ALTERNATIVES.get(branch, [])
                found_alt = None
                for alt in alts:
                    if alt in available_branches:
                        found_alt = alt
                        break
                
                if found_alt:
                    result['details']['found'].append(found_alt)
                    result['details']['alternatives_used'][branch] = found_alt
                else:
                    result['details']['missing'].append(branch)
        
        f.close()
        
        if len(result['details']['missing']) == 0:
            result['status'] = 'PASS'
        else:
            result['status'] = 'FAIL'
            
    except Exception as e:
        result['status'] = 'FAIL'
        result['details']['error'] = str(e)
        
    return result


def gate_g3_zpeak(filepath: str, max_events: int = 100000) -> dict:
    """G3: Z peak sanity check."""
    import uproot
    import numpy as np
    
    result = {
        'gate': 'G3',
        'name': 'Z Peak Sanity',
        'status': 'UNKNOWN',
        'details': {}
    }
    
    try:
        f = uproot.open(filepath)
        events = f['Events']
        
        arrays = events.arrays(
            ['Muon_pt', 'Muon_eta', 'Muon_phi', 'Muon_charge', 'Muon_mass', 'nMuon'],
            entry_stop=max_events,
            library='np'
        )
        
        masses = []
        n_candidates = 0
        
        for i in range(len(arrays['nMuon'])):
            if arrays['nMuon'][i] >= 2:
                pt1, pt2 = arrays['Muon_pt'][i][0], arrays['Muon_pt'][i][1]
                eta1, eta2 = arrays['Muon_eta'][i][0], arrays['Muon_eta'][i][1]
                phi1, phi2 = arrays['Muon_phi'][i][0], arrays['Muon_phi'][i][1]
                q1, q2 = arrays['Muon_charge'][i][0], arrays['Muon_charge'][i][1]
                m1, m2 = arrays['Muon_mass'][i][0], arrays['Muon_mass'][i][1]
                
                if q1 * q2 < 0:
                    px1 = pt1 * np.cos(phi1)
                    py1 = pt1 * np.sin(phi1)
                    pz1 = pt1 * np.sinh(eta1)
                    e1 = np.sqrt(px1**2 + py1**2 + pz1**2 + m1**2)
                    
                    px2 = pt2 * np.cos(phi2)
                    py2 = pt2 * np.sin(phi2)
                    pz2 = pt2 * np.sinh(eta2)
                    e2 = np.sqrt(px2**2 + py2**2 + pz2**2 + m2**2)
                    
                    mll_sq = (e1 + e2)**2 - (px1 + px2)**2 - (py1 + py2)**2 - (pz1 + pz2)**2
                    if mll_sq > 0:
                        mll = np.sqrt(mll_sq)
                        if 60 < mll < 120:
                            masses.append(mll)
                            n_candidates += 1
        
        f.close()
        
        if len(masses) > 100:
            masses = np.array(masses)
            peak_region = masses[(masses > 85) & (masses < 97)]
            
            if len(peak_region) > 50:
                peak_mass = np.median(peak_region)
                result['details']['peak_mass'] = float(peak_mass)
                result['details']['n_candidates'] = n_candidates
                result['details']['events_processed'] = len(arrays['nMuon'])
                
                if 88 < peak_mass < 94:
                    result['status'] = 'PASS'
                else:
                    result['status'] = 'FAIL'
                    result['details']['reason'] = f'Peak at {peak_mass:.1f} GeV, expected ~91 GeV'
            else:
                result['status'] = 'FAIL'
                result['details']['reason'] = 'Insufficient candidates in peak region'
        else:
            result['status'] = 'FAIL'
            result['details']['reason'] = f'Only {len(masses)} Z candidates found'
            
    except Exception as e:
        result['status'] = 'FAIL'
        result['details']['error'] = str(e)
        
    return result


def print_result(result: dict):
    """Pretty-print a gate result."""
    status_symbol = '✓' if result['status'] == 'PASS' else '✗'
    print(f"\n{'='*60}")
    print(f"Gate {result['gate']}: {result['name']}")
    print(f"Status: [{status_symbol}] {result['status']}")
    print(f"{'='*60}")
    
    for key, value in result['details'].items():
        if isinstance(value, list) and len(value) > 5:
            print(f"  {key}: [{len(value)} items]")
            for item in value[:3]:
                print(f"    - {item}")
            print(f"    ... and {len(value)-3} more")
        elif isinstance(value, dict):
            print(f"  {key}:")
            for k, v in value.items():
                print(f"    {k}: {v}")
        else:
            print(f"  {key}: {value}")


def main():
    parser = argparse.ArgumentParser(description='ZMET-NID Phase 0 Feasibility Check')
    parser.add_argument('--local', type=str, help='Path to local NanoAOD file')
    args = parser.parse_args()
    
    print("Checking dependencies...")
    uproot_ok, uproot_ver = check_uproot_available()
    if not uproot_ok:
        print("ERROR: uproot not installed. Run: pip install uproot awkward")
        sys.exit(1)
    print(f"  uproot version: {uproot_ver}")
    
    if args.local:
        filepath = args.local
        print(f"\nUsing local file: {filepath}")
    else:
        filepath = TEST_FILE_XROOTD
        print(f"\nUsing XRootD: {filepath}")
        print("(This may take a moment...)")
    
    results = []
    
    print("\n" + "="*60)
    print("Running Gate G1: File Connectivity...")
    g1 = gate_g1_connectivity(filepath)
    print_result(g1)
    results.append(g1)
    
    if g1['status'] != 'PASS':
        print("\n[ABORT] G1 failed - cannot proceed without file access")
        sys.exit(1)
    
    print("\nRunning Gate G2: Branch Availability...")
    g2 = gate_g2_branches(filepath)
    print_result(g2)
    results.append(g2)
    
    if g2['status'] != 'PASS':
        print("\n[WARNING] G2 failed - missing required branches")
        print("Consider checking MiniAOD or adjusting analysis scope")
    
    print("\nRunning Gate G3: Z Peak Sanity Check...")
    g3 = gate_g3_zpeak(filepath)
    print_result(g3)
    results.append(g3)
    
    print("\n" + "="*60)
    print("PHASE 0 SUMMARY")
    print("="*60)
    
    all_pass = all(r['status'] == 'PASS' for r in results)
    
    for r in results:
        symbol = '✓' if r['status'] == 'PASS' else '✗'
        print(f"  [{symbol}] {r['gate']}: {r['name']} - {r['status']}")
    
    print()
    if all_pass:
        print("[PROCEED] All gates passed. Ready for Phase 1 execution.")
    else:
        print("[HOLD] Some gates failed. Review before proceeding.")
    
    return 0 if all_pass else 1


if __name__ == '__main__':
    sys.exit(main())
