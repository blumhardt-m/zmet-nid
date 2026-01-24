# ZMET-NID: MET Non-Identifiability Study

**Title:** Limits of missing transverse momentum interpretation with public LHC data: a Z→μ⁺μ⁻ identifiability study using CMS Open Data

**Status:** Protocol frozen, manuscript 80% complete, awaiting Phase 0 execution

## Project Summary

This study quantifies the degree to which MET tails in Z→μμ events are non-identifiable under reduced public-tier observability (NanoAOD), even when standard global validation checks pass.

**Core question:** Can multiple incompatible nuisance explanations fit the same MET distribution while passing standard validation?

## Directory Structure

```
ZMET_NID/
├── manuscript/          # Paper sections (JINST target)
├── protocol/            # Frozen protocol v0.4
├── code/
│   ├── phase0/          # Feasibility gate scripts
│   ├── analysis/        # Main analysis code
│   └── utils/           # Helper functions
├── data/
│   ├── raw/             # Downloaded NanoAOD files
│   └── processed/       # Processed outputs
├── figures/             # Generated plots
├── notes/               # Working notes
└── literature/          # Reference papers
```

## Execution Checklist

### Phase 0 (Feasibility Gates)
- [ ] G1: File connectivity test
- [ ] G2: Branch availability audit
- [ ] G3: Z peak sanity check

### Phase 1 (Analysis)
- [ ] Event selection implementation
- [ ] Nuisance family implementation
- [ ] Train/test split
- [ ] Fitting procedure
- [ ] Stratified tests

### Phase 2 (Writing)
- [ ] Fill Section 7 with results
- [ ] Write Section 3 (selection details)
- [ ] Write Sections 8-9 (discussion/conclusions)
- [ ] Final review

## Key Files

- `protocol/PROTOCOL_v0.4.md` - Frozen analysis protocol
- `manuscript/00_abstract.md` - Publication-ready abstract
- `code/phase0/feasibility_check.py` - Day-1 feasibility script

## Kill Criteria

Abort if:
1. Required NanoAOD branches absent
2. Cannot reproduce Z peak (not at 91±2 GeV)
3. Event yield < 50k after selection
4. Runtime > 2 weeks on available hardware

## Target Venue

Journal of Instrumentation (JINST)
