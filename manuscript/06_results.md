# 6. Results

After the selection described in Section 3, a total of 233,524 Z→μ⁺μ⁻
candidate events remain from 2,147,195 processed events. All quantitative
comparisons in this section use this selected sample.

The three nuisance families examined here occupy distinct positions along an
identifiability spectrum: Family A represents the injected truth, Family B
exhibits deep inclusive degeneracy, and Family C lies near the boundary of
inclusive discriminability. Figure 5 provides a summary of this landscape at the
end of this section.

The results are organised as follows. Section 6.1 presents the inclusive
missing transverse momentum comparison, which establishes the working point
for the identifiability demonstration. Section 6.2 examines the same
comparison after conditioning on pT(Z), the primary recoil variable.
Section 6.3 presents the jet-multiplicity stratification. Section 6.4
documents the robustness of the inclusive result to MET binning. Section 6.5
presents a physically grounded nuisance variation derived from CMS
unclustered-energy systematic branches.

## 6.1 Inclusive MET degeneracy

The analysis sample consists of 233,524 Z→μ⁺μ⁻ candidates selected from
CMS Open Data (DoubleMuon, 2016). Both nuisance families were fitted
independently to the inclusive MET histogram. Family A (recoil-response
distortion) was used to construct the target distribution at the working
point α = 0.05, representing a 5% MET scale distortion growing linearly
with min(pT(Z), 100 GeV)/100. Family B (topology-mixture reweighting) was
then fitted to the same target by grid search over β ∈ [−0.40, +0.40].

The best-fit parameters and inclusive fit quality are summarised in Table 1.

**Table 1.** Inclusive MET fit quality. The target distribution was constructed
by applying Family A at the working point α = 0.05. All families were fitted
independently to this target using a χ² statistic with (obs + 1) denominator.
Degrees of freedom are (number of histogram bins) − 1 fitted parameter.
Family A, which generated the target, trivially achieves χ² = 0.

| Model    | Best parameter | χ²    | dof | p-value |
|----------|----------------|-------|-----|---------|
| Family A | α = +0.050     |  0.00 | —   | —       |
| Family B | β = +0.290     |  3.13 | 11  | 0.989   |
| Family C | γ = +0.200     | 18.28 | 11  | 0.075   |

Family B achieves χ²(B) = 3.13 on 11 degrees of freedom, corresponding to
p = 0.989. The inclusive projection provides no statistical basis to reject
the topology-mixture model at this working point. The result is shown
visually in Fig. 1, which presents the inclusive MET distributions for the
target, Family A, and Family B on a logarithmic scale, together with the
ratio of each model to the target in the lower panel. The three distributions
are visually indistinguishable across the bulk of the spectrum (MET < 100 GeV);
a small excess from Family B appears in the sparse high-MET tail above
100 GeV, which accounts for the residual χ² = 3.13.

These results demonstrate that the inclusive MET distribution, despite
containing 233,524 selected events, admits a topology-mixture explanation
with goodness of fit statistically equivalent to the injected recoil model.

As an empirical noise-floor check, the inclusive χ² statistic was recomputed
on 1,000 bootstrap resamplings of the target MET distribution using the same
12-bin histogram and the same χ² metric. Each bootstrap sample draws 233,524
events with replacement from the target and is compared against the target
histogram directly. The resulting null distribution has a median of χ² = 10.6
and spans the range 2.2–29.8. The χ²(B) = 3.13 obtained for the
topology-mixture family falls at the 1.2nd percentile of this distribution
(Fig. 2), demonstrating that the inclusive agreement of Family B is not merely
within the formal rejection threshold but is substantially better than the
typical level of sampling variability intrinsic to the target histogram itself.
The topology-mixture model, with its smooth optimised-weight histogram, provides
a fit quality well inside the noise floor of the target distribution.

![Inclusive MET distributions for the target (black), Family A (dashed red), Family B (dotted blue), and Family C (dash-dot green) on a logarithmic scale. The lower panel shows the ratio of each model to the target. Family B achieves χ²(B) = 3.13 (dof = 11, p = 0.989); Family C achieves χ²(C) = 18.28 (dof = 11, p = 0.075). All families are fitted to the inclusive histogram only. The sample comprises 233,524 Z→μ⁺μ⁻ candidates from CMS Open Data (DoubleMuon 2016).](figures/inclusive_fit/met_familyA_vs_familyB.png){width=85%}

![Empirical null distribution of the inclusive χ² statistic from 1,000 bootstrap resamplings of the target MET histogram (blue). Each bootstrap sample draws 233,524 events with replacement from the target. The dashed red line marks χ²(B) = 3.13, the inclusive fit for Family B. The null median is χ² = 10.6 (dotted gray). Family B falls at the 1.2nd percentile of the null, confirming that the inclusive agreement is well within the sampling noise floor of the target distribution.](figures/inclusive_fit/null_chi2_distribution.png){width=75%}

## 6.2 Recoil-stratified MET distributions

To probe whether the inclusive degeneracy conceals differences in the
internal event structure, both models were evaluated—without refitting—in
four pT(Z) strata following the protocol binning: [0, 10), [10, 30), [30, 60),
and ≥60 GeV. The separation between the two nuisance families is quantified
using Δχ² = χ²(B) − χ²(A), evaluated within each stratum while keeping the
nuisance parameters fixed at the values obtained from the inclusive fit.
Because Family A generated the target distribution, χ²(A) ≈ 0 in each
stratum and Δχ² ≈ χ²(B).

**Table 2.** Per-stratum discrimination by recoil activity. Δχ² = χ²(B) − χ²(A)
is evaluated in each pT(Z) bin using the best-fit parameters from the inclusive
fit only. Both families were fitted to the inclusive distribution; the
stratified comparison uses the inclusive best-fit parameters without refitting.

| pT(Z) bin (GeV)  | N      | Δχ²   |
|------------------|--------|-------|
| < 10             | 92,608 |   392 |
| 10 ≤ pT(Z) < 30  | 89,387 |    95 |
| 30 ≤ pT(Z) < 60  | 34,312 |   513 |
| ≥ 60             | 17,217 | 1,432 |

The overall trend is that Δχ² grows with recoil activity, with the largest
separation in the highest pT(Z) bin (Δχ² = 1,432) and the smallest in the
10–30 GeV bin (Δχ² = 95). This pattern is consistent with the design of
Family A, which produces MET distortions proportional to pT(Z): events with
higher pT(Z) are shifted more, and Family B—which does not scale MET with
pT(Z)—cannot reproduce that pattern when the distributions are conditioned
on recoil.

The Δχ² in the lowest pT(Z) bin (< 10 GeV, Δχ² = 392) is larger than in
the 10–30 GeV bin. At very low pT(Z) the recoil effect is weak, so both
families agree on the MET shape; however, the topology-mixture reweighting
perturbs the event composition in this stratum even when the MET-scale
effect is negligible, producing a larger discrepancy than in the bin where
the two effects partially cancel.

The per-stratum MET distributions are shown in Fig. 3. In the three lower
pT(Z) bins, Family A and Family B are visually close to the target. In the
highest pT(Z) bin, Family B is systematically elevated relative to the
target at low MET values and slightly suppressed at higher MET, while
Family A tracks the target throughout.

![MET distributions in four pT(Z) strata: < 10 GeV (N = 92,608), 10–30 GeV (N = 89,387), 30–60 GeV (N = 34,312), and ≥60 GeV (N = 17,217). Target (black), Family A (dashed red), Family B (dotted blue), and Family C (dash-dot green) are shown at inclusive best-fit parameters without per-stratum refitting. The per-stratum discrimination Δχ²(B) = χ²(B) − χ²(A) grows with recoil activity, reaching Δχ² = 1,432 in the highest bin. Family C produces Δχ²(C) = 58 in the highest pT(Z) bin.](figures/stratified_tests/met_familyA_vs_familyB_by_zpt.png){width=100%}

## 6.3 Jet-multiplicity stratification

As a complementary test, both models were evaluated after conditioning on
the number of jets with pT > 30 GeV, cleaned against the selected Z muons
using ΔR > 0.4. The three jet categories follow the protocol binning:
0-jet, 1-jet, and ≥2-jet. Results are shown in Table 3 and Fig. 4.

**Table 3.** Per-stratum discrimination by jet multiplicity. Δχ² = χ²(B) − χ²(A)
evaluated using best-fit parameters from the inclusive fit.

| Jet bin | N       | Δχ²   |
|---------|---------|-------|
| 0-jet   | 171,385 | 1,405 |
| 1-jet   |  44,706 | 1,421 |
| ≥2-jet  |  17,433 | 3,415 |

Because the topology-mixture model modifies the relative weights of jet
multiplicity categories globally, conditioning on jet count directly tests
whether the model preserves the observed event population across topologies.

The jet-multiplicity stratification tests population consistency across event
topologies rather than the MET–recoil relationship. Family B assigns global
event weights of 1, (1 + β), and (1 + 2β) to the 0-jet, 1-jet, and ≥2-jet
categories respectively, normalised to preserve the inclusive event count.
At best-fit β = +0.290, the effective per-event weights after normalisation
are approximately 0.91, 1.17, and 1.44 for the three categories. When the
distribution is conditioned on jet multiplicity, this reweighting manifests
as an inflated event count in the ≥2-jet stratum relative to the target, and
a correspondingly depleted count in the 0-jet stratum.

Family A, by contrast, applies a continuous MET scale factor that does not
alter the relative population of each jet stratum. The conditional event
counts are therefore preserved by Family A and violated by Family B.
The ≥2-jet bin produces the largest per-stratum Δχ² (3,415), reflecting the
44% count inflation relative to the target. The discrimination in the 0-jet
and 1-jet bins (Δχ² ≈ 1,400 each) reflects the corresponding count deficit
and surplus in those categories.

The per-event Δχ² (Δχ² / N) is 0.0082, 0.032, and 0.196 for the 0-jet,
1-jet, and ≥2-jet bins respectively, confirming that the ≥2-jet stratum
provides the strongest population-level discrimination per event.

![MET distributions in three jet-multiplicity categories: 0-jet (N = 171,385), 1-jet (N = 44,706), and ≥2-jet (N = 17,433). Target (black), Family A (dashed red), Family B (coloured dotted), and Family C (dash-dot green) are evaluated at inclusive best-fit parameters without refitting. Family B inflates the ≥2-jet event count by approximately 44% relative to the target, yielding Δχ²(B) = 3,415. Family C produces Δχ²(C) ≤ 26 across all jet bins. Family A achieves χ²(A) ≈ 0 throughout.](figures/stratified_tests/met_familyA_vs_familyB_by_njets.png){width=90%}

Together, Tables 2 and 3 show that both the recoil-stratified and
jet-stratified projections expose clear discrepancies in Family B that are
absent from the inclusive distribution. The two stratifications probe
different latent assumptions: pT(Z) conditioning tests the MET–recoil
relationship, while jet-multiplicity conditioning tests the population model.

## 6.4 Robustness to MET binning

The inclusive degeneracy demonstrated in Section 6.1 was evaluated under
two alternative MET binning schemes in addition to the baseline (12 bins):
a coarser binning with 9 bins merging adjacent low-MET intervals, and a
tail-emphasised binning with modified boundaries above 50 GeV. In all three
cases, χ²(A) remained zero and the p-value for Family B exceeded 0.96.
The inclusive degeneracy at the working point α = 0.05 persisted under
modest changes in MET binning, including coarser and tail-emphasised bin
choices, with p(B) > 0.96 across all tested binning schemes. Detailed
binning robustness results are reported in Supplementary Table S1.

Taken together, these tests demonstrate that while the inclusive MET
distribution admits multiple nuisance explanations at modest distortion
amplitudes, conditional projections restore discriminability by exposing
violations of either recoil–response structure or event-population
consistency.

## 6.5 Physically grounded public-tier nuisance variation

To test whether the inclusive degeneracy demonstrated above extends beyond the
phenomenological nuisance families introduced in Section 4, a third nuisance family
is constructed based on detector-level systematic variations exposed directly in
the CMS NanoAOD format.

NanoAOD provides event-level shifts corresponding to the CMS unclustered-energy
systematic variation through the branches `MET_MetUnclustEnUpDeltaX` and
`MET_MetUnclustEnUpDeltaY`. These quantities represent the change in the
reconstructed missing transverse momentum vector when the unclustered-energy
component of the event is shifted upward within the CMS reconstruction framework.

Family C — unclustered-energy vector perturbation — is defined by applying a scaled
version of this shift to the nominal MET vector:

    MET_C = sqrt((MET_x + γ ΔX)² + (MET_y + γ ΔY)²)

where MET_x = MET_pt × cos(φ), MET_y = MET_pt × sin(φ), and (ΔX, ΔY) are the
NanoAOD unclustered-energy variation components. The parameter γ scales the
magnitude of the systematic shift, with γ = 1 corresponding to the full
unclustered-energy-up variation provided by CMS. The parameter search range is
γ ∈ [−3.0, +3.0] in steps of 0.05.

Fitting Family C to the inclusive target distribution produced by the
recoil-response injection (α = 0.05) yields a best-fit value γ = +0.200,
corresponding to approximately 20% of the full CMS unclustered-energy shift.
The resulting inclusive fit quality is χ²(C) = 18.28 (dof = 11, p = 0.075).

This value lies above the nominal rejection threshold (p = 0.05), indicating that
the unclustered-energy perturbation is not rejected by the inclusive projection at
this working point. However, the fit quality is substantially worse than that
obtained for the topology-mixture family (Family B), which achieves χ² = 3.13 and
p = 0.989.

Conditional projections again expose discrepancies between the nuisance family and
the injected recoil model. The recoil-stratified comparison produces Δχ²(C) of
{17, 5, 19, 58} across the four pT(Z) strata.

The contrast with the topology-mixture family is particularly pronounced in the
highest recoil bin. In the pT(Z) ≥ 60 GeV stratum, the topology-mixture model
produces Δχ²(B) = 1,432, whereas the unclustered-energy perturbation yields
Δχ²(C) = 58. The ratio of these values (≈ 25) quantifies how much more
efficiently the recoil stratification exposes the topology-mixture mechanism
relative to the detector-level perturbation.

The jet-multiplicity stratification yields Δχ²(C) of {11, 26, 24} for the 0-jet,
1-jet, and ≥2-jet categories respectively.

The smaller stratified discrepancies observed for Family C do not indicate that it
provides a superior explanation of the target distribution. Rather, they reflect a
difference in failure mode. The unclustered-energy perturbation modifies the MET
vector while preserving the event-topology composition of the sample. As a result,
jet-multiplicity conditioning does not directly expose a population inconsistency,
in contrast to the topology-mixture family whose global reweighting violates the
relative event counts across jet categories. The discriminating power of a given
stratification variable therefore depends on which latent assumption of the nuisance
model it probes.

These values are substantially smaller than those obtained for the
topology-mixture family (Section 6.3), reflecting the fact that the
unclustered-energy perturbation modifies the MET vector without altering the
relative population of jet multiplicity categories. Consequently, the
jet-stratified diagnostics expose Family C less strongly than Family B, which
explicitly distorts the event-topology mixture.

Figure 1 shows the inclusive MET comparison including Family C alongside Families A
and B. The unclustered-energy perturbation closely tracks the target distribution
over most of the spectrum, but exhibits mild deviations in the intermediate MET
region that account for the larger inclusive χ².

Taken together, these results indicate that physically grounded perturbations
derived from CMS systematic variation branches can approach inclusive acceptability
while remaining distinguishable through conditional projections. The strength of
the degeneracy therefore depends on the underlying nuisance mechanism:
topology-mixture distortions produce deep inclusive degeneracy, whereas
detector-level perturbations such as unclustered-energy shifts exhibit weaker,
partially recoverable degeneracy.

**Table 4.** Summary of identifiability across all three nuisance families at
the working point α = 0.05. Inclusive p-value and maximum per-stratum Δχ² are
reported for each family. A higher inclusive p-value indicates deeper degeneracy
with the injected target (better inclusive fit); a higher max Δχ² indicates
stronger conditional exposure.

| Family | Mechanism                 | Inclusive χ² | p-value | max Δχ² |
|--------|---------------------------|--------------|---------|---------|
| A      | Injected truth            | 0.00         | —       | 0       |
| B      | Topology-mixture          | 3.13         | 0.989   | 3,415   |
| C      | Unclustered-energy shift  | 18.28        | 0.075   | 58      |

![Identifiability landscape of the three nuisance families. The horizontal axis shows the inclusive p-value for goodness-of-fit to the injected target: higher values indicate better inclusive fit and therefore deeper degeneracy risk. The vertical axis shows the maximum conditional discrepancy across all recoil and jet stratifications, displayed on a symlog scale (linear below Δχ² = 100, logarithmic above). The dotted vertical line marks the p = 0.05 inclusive rejection threshold. Family B occupies a regime of deep inclusive degeneracy (p = 0.989) with catastrophic conditional failure (Δχ² = 3,415). Family C lies near the boundary of inclusive discriminability (p = 0.075) with modest conditional discrepancies (Δχ² = 58). The injected recoil model (Family A) anchors the upper-right corner at p ≈ 1 and Δχ² = 0. Shaded regions illustrate practical regimes of nuisance-model identifiability observed in this analysis: inclusive degeneracy (Δχ² < 10), marginal identifiability (10–100), and clear separation (Δχ² > 100). The dashed horizontal line indicates the χ² rejection threshold corresponding to p = 0.05 for 11 degrees of freedom. Here Δχ² denotes the difference between the stratified χ² values of the competing nuisance families evaluated at the inclusive best-fit parameters (Section 5.3).](figures/inclusive_fit/identifiability_landscape.png){width=80%}
