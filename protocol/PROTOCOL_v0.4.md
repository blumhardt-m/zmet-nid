# ZMET-NID Protocol v0.4 (FROZEN)

**Frozen:** 2026-01-24
**Do not modify without explicit version increment**

## 1. Data Scope

- **Format:** CMS 2016 NanoAOD
- **Primary dataset:** DoubleMuon
- **MC:** DY→μμ + jets
- **Validated runs:** CMS JSON with checksums
- **Expected yield:** O(10⁶) Z→μμ candidates pre-selection
- **Minimum floor:** 50k selected events

## 2. Event Selection (FROZEN)

```
- Exactly 2 opposite-sign muons
- pT(μ1) ≥ 25 GeV, pT(μ2) ≥ 15 GeV
- |η(μ)| ≤ 2.4
- Muon ID: Tight (or best available)
- Isolation: PFIso < 0.15
- M(μμ) ∈ [80, 100] GeV
- Standard event filters per CMS guidance
```

## 3. Primary Observable

`MET_pt` (magnitude of missing transverse momentum)

## 4. Stratification Variables

| Variable | Binning |
|----------|---------|
| Pileup (PV_npvs) | Terciles |
| Jets (N_jets, pT > 30 GeV) | {0, 1, ≥2} |
| Z pT | [0-10), [10-30), [30-60), [60+) GeV |
| MET φ | 8 equal bins |

## 5. Nuisance Families (FROZEN)

### Family A: Jet-correlated MET scaling
```
MET' = MET × (1 + α × N_jets + γ × max(0, pT_lead - 30)/30)
Parameters: θ_A = {α, γ}
Bounds: α ∈ [0, 0.3], γ ∈ [0, 0.3]
```

### Family B: Pileup-correlated MET broadening
```
MET' = MET + |Δ|
where |Δ| ~ Rayleigh(β₀ + β₁ × max(0, nPV - 20))
Direction: uniform random in φ
Parameters: θ_B = {β₀, β₁}
Bounds: β₀ ∈ [0, 20], β₁ ∈ [0, 2]
```

## 6. Train/Test Split

- **Train:** Even lumisections
- **Test:** Odd lumisections

## 7. Fitting Metrics

```
D = 0.5 × W₁(MET histogram) + 0.5 × KS(MET CDF)
```

## 8. Threshold Determination (Preregistered)

On pilot data (first 100k events, training split):
1. Draw 100 bootstrap resamples of observed MET distribution
2. Compute pairwise W₁ and KS distances
3. Set τ_W = median(W₁), τ_KS = median(KS)
4. Freeze thresholds before test-split unblinding

## 9. Decision Criteria

**Global pass:** KS p > 0.05 AND W₁ < τ_W

**Non-identifiability declared if ALL:**
1. Both families pass global tests
2. |D(A) - D(B)| < 0.05 on test set
3. Each family has ≥1 parameter with 95% CI excluding zero

**Stratified separation:** Any structure test rejects one family (p < 0.01/24 Bonferroni) while accepting other (p > 0.05)

## 10. Outcomes

1. **Non-identifiable globally, separated structurally** → "Standard validation passes underdetermined inference"
2. **Non-identifiable globally AND structurally** → "Available observables cannot distinguish mechanisms"
3. **One family fails globally** → "Identifiable with these observables"
4. **Both fail globally** → "Neither model adequate; information bottleneck finding"

## 11. Kill Criteria

- Required branches absent in NanoAOD and MiniAOD
- Cannot reproduce Z peak (not at 91±2 GeV)
- Event yield < 50k after selection
- Runtime > 2 weeks on available hardware
- Fitting fails to converge for either family
- Results change qualitatively under reasonable binning variations
