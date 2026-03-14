# ZMET-NID: MET Non-Identifiability Study

This repository contains the reproducible analysis supporting the paper:

"Nuisance Parameter Non-Identifiability in Compressed Collider Observables."

The project demonstrates how compression of high-dimensional collider events into
reduced observables (such as missing transverse energy) can produce degeneracy between
distinct nuisance mechanisms in inference pipelines.

## Repository Structure

    analysis/        Analysis scripts
    data/            Data manifests and SHA-256 checksums (no raw data)
    manuscript/      Paper draft and figures
    protocol/        Locked analysis protocol
    outputs/         Result files and plots

## Analysis Summary

Z→μ⁺μ⁻ candidate events are selected from CMS Open Data (DoubleMuon 2016, NanoAOD
tier, 13 TeV). Two phenomenological nuisance families — topology-mixture reweighting
(Family B) and recoil-response scale distortion (Family A) — are fitted to an injected
target MET distribution. The inclusive MET projection is found to be insufficient to
discriminate between the two families at a physically motivated working point.
Stratified projections conditioned on pT(Z) and jet multiplicity restore discrimination,
yielding max Δχ² = 3,415.

A third nuisance family (Family C, derived from CMS unclustered-energy systematic
branches) is tested in the boundary regime of inclusive discriminability (χ² = 18.28,
p = 0.075), demonstrating that identifiability is graded rather than binary.

## Protocol Governance

The analysis protocol is frozen in `protocol/PROTOCOL.locked.yaml`. All analysis steps
were executed against this locked protocol. SHA-256 checksums for input data files and
the protocol document itself are recorded in `data/manifests/`.

## Key Files

- `protocol/PROTOCOL.locked.yaml` — Frozen analysis protocol
- `manuscript/ZMET_NID_paper.pdf` — Compiled manuscript
- `manuscript/00_abstract.md` through `manuscript/12_references.md` — Paper sections
- `code/analysis/` — Analysis scripts (phase1_baseline.py, fit_minimal_nuisance_families.py, etc.)
- `data/manifests/` — SHA-256 checksums and data source records

## Target Venue

Journal of Instrumentation (JINST)
