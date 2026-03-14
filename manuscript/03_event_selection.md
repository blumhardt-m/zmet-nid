# 3. Event Selection

## 3.1 Input observables

The selection and baseline analysis read the following NanoAOD branches from the `Events` tree:

- **Muon kinematics and identity:** `Muon_pt`, `Muon_eta`, `Muon_phi`, `Muon_mass`, `Muon_charge`
- **Muon quality:** `Muon_tightId` (boolean), `Muon_pfRelIso04_all` (PF relative isolation, cone Δ*R* = 0.4)
- **Missing transverse momentum:** `MET_pt`
- **Jets:** `Jet_pt`, `Jet_eta`, `Jet_phi`
- **Event count:** `nMuon`

No lower-level detector quantities are used. The analysis operates entirely at the NanoAOD reconstructed-object tier. The muon identification branch is resolved at runtime in order of preference: `Muon_tightId`, `Muon_mediumId`, `Muon_looseId`; the isolation branch is resolved in order: `Muon_pfRelIso04_all`, `Muon_pfRelIso03_all`, `Muon_miniPFRelIso_all`. For the dataset used in this analysis, `Muon_tightId` and `Muon_pfRelIso04_all` are present and are used throughout.

## 3.2 Muon selection and dimuon construction

Per-muon quality cuts are applied before any event-level selection. A muon is retained if it satisfies all of the following requirements simultaneously:

- *p*_T > 15 GeV
- |η| < 2.4
- `Muon_tightId` = 1 (CMS Tight muon identification)
- `Muon_pfRelIso04_all` < 0.15 (PF relative isolation, Δ*R* = 0.4 cone)

Both muons are required to satisfy p*_T* > 15 GeV and |η| < 2.4. In addition, the leading muon in the selected pair must satisfy p*_T* ≥ 25 GeV. All muons passing the per-muon cuts are retained as candidates. No multiplicity requirement is applied at this stage. Four-momentum vectors are constructed from (`Muon_pt`, `Muon_eta`, `Muon_phi`, `Muon_mass`) using the `vector` library Momentum4D representation.

Dimuon pairs are formed by taking all unique two-muon combinations within each event and retaining only opposite-sign pairs. For events containing multiple quality-selected muons, all opposite-sign combinations are evaluated. The invariant mass of each OS pair is computed from the four-vector sum.

If multiple opposite-sign muon pairs satisfy the baseline selection, the Z candidate is defined as the pair whose invariant mass is closest to the nominal Z boson mass (91.1876 GeV). This "best pair" criterion is applied event by event without iterating over alternative selections. The leading muon p*_T* of the best pair — defined as the higher of the two muon transverse momenta — must satisfy p*_T*(μ_1) ≥ 25 GeV. Events in which no OS pair exists, or in which the best-pair leading muon fails this threshold, are rejected.

## 3.3 Z-candidate definition

The Z-candidate mass is taken from the four-vector sum of the best OS pair selected in Section 3.2. A mass window requirement of 80 GeV ≤ *m*(μ⁺μ⁻) ≤ 100 GeV is applied as a hard selection cut. Events whose best-pair mass falls outside this window are excluded from the analysis sample. This window defines the Z→μ⁺μ⁻ signal region used in all downstream inclusive and stratified comparisons.

The dimuon transverse momentum p*_T*(Z) is taken as the magnitude of the transverse component of the four-vector sum of the two selected muons. It serves as the primary recoil variable for the pT(Z) stratification in Section 6.2.

A separate histogram-mode estimator in the narrower window 85–97 GeV is used in the Phase 0 feasibility check to confirm the Z peak position without fitting assumptions. This narrower window is a diagnostic convention and does not define the analysis sample.

## 3.4 Jet cleaning and jet-multiplicity definition

Hadronic activity is quantified by counting jets reconstructed by the CMS particle-flow algorithm. A jet is included in the multiplicity count if it satisfies:

- p*_T*(jet) > 30 GeV
- Δ*R*(jet, μ_1) > 0.4 **and** Δ*R*(jet, μ_2) > 0.4

where μ_1 and μ_2 are the two muons forming the selected Z candidate, and Δ*R* = √(Δη² + Δφ²). Jets that fail either cleaning condition are removed from the count before applying the p*_T* threshold.

This cleaning step prevents the particle-flow jet algorithm from counting the selected muons as jets. Any jet within a cone of Δ*R* < 0.4 around either Z-candidate muon is removed, ensuring that the jet multiplicity reflects genuine hadronic recoil activity rather than jet-muon object overlap.

Events are assigned to one of three exclusive jet-multiplicity categories:

- **0-jet:** no jets with p*_T* > 30 GeV survive the cleaning requirement
- **1-jet:** exactly one such jet survives
- **≥2-jet:** two or more such jets survive

These three categories are used for the jet-multiplicity stratification in Section 6.3.

## 3.5 Analysis sample used in Section 6

The selection described in Sections 3.2–3.4 was applied to the full CMS Open Data DoubleMuon 2016 NanoAOD dataset. A total of 2,147,195 events were processed. After applying all muon quality cuts, the opposite-sign requirement, the leading muon p*_T* threshold of 25 GeV, and the Z-mass window of 80–100 GeV, 233,524 events were retained. The overall selection efficiency is approximately 10.9%.

This sample of 233,524 Z→μ⁺μ⁻ candidates constitutes the analysis sample for all comparisons in Section 6. The same selected events are used for the inclusive MET comparison (Section 6.1), the pT(Z)-stratified comparisons (Section 6.2), the jet-multiplicity-stratified comparisons (Section 6.3), and the MET binning robustness checks (Section 6.4). No additional event-level requirements are applied within the results section; the stratification into sub-samples uses the dimuon p*_T* and jet-multiplicity variables already computed during the selection. This selected sample is used for all inclusive and stratified comparisons reported in Section 6.
