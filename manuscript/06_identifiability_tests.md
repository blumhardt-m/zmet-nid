# 6. Identifiability Test Design

This section describes the statistical framework used to test whether the two
nuisance families defined in Section 4 can be distinguished using the inclusive
MET distribution, and how conditional projections are used to recover
discriminating information absent from the inclusive view. The analysis uses the
full selected sample of 233,524 Z→μ⁺μ⁻ events described in Section 3; no
train/test split is applied.

## 6.1 Injection Framework

A target distribution is constructed by applying Family A (recoil-response
distortion) to the observed MET values at a fixed injection parameter α₀:

    MET_target(i) = MET(i) × (1 + α₀ × clip(pT(Z)_i, 0, 100) / 100)

This target represents the observed distribution as it would appear under a
known, controlled nuisance distortion. Family A trivially recovers the target
because it generated it; the nontrivial question is whether Family B
(topology-mixture reweighting) can also reproduce the same target on the inclusive
MET histogram, despite operating through a mechanistically different mechanism.

The injection framework has two advantages over a direct comparison of models to
data. First, it provides a ground truth: the injected distortion amplitude α₀ is
known exactly, so the quality of Family A's recovery is a calibration check.
Second, it isolates the identifiability question from the absolute normalization
of the data: both families are evaluated against the same target histogram, so
differences in their fit quality reflect differences in shape-matching ability,
not overall normalisation.

## 6.2 Inclusive Fit Criterion

Both families are fitted independently to the target histogram by grid search
over their respective parameters. The fit metric is:

    χ²(model, target) = Σ_k (h_k^model − h_k^target)² / (h_k^target + 1)

where the sum runs over all MET histogram bins k, h_k denotes the bin count, and
the (h_k^target + 1) denominator provides numerical stability in bins with low
target occupancy. This form is analogous to a Pearson χ² but uses the target bin
count rather than the model prediction in the denominator, avoiding sensitivity to
model parameters in the uncertainty estimate.

The MET histogram uses 12 variable-width bins spanning 0–200 GeV, with boundaries
at [0, 5, 10, 15, 20, 30, 40, 50, 60, 80, 100, 150, 200] GeV. Finer binning at
low MET (where most events lie) provides sensitivity to shape changes in the bulk
of the distribution; coarser binning above 50 GeV reflects the sparse high-MET
tail.

Family A is searched over α ∈ [−0.30, +0.30] in steps of 0.01. Family B is
searched over β ∈ [−0.40, +0.40] in steps of 0.01. The best-fit parameter for
each family is the value minimising χ² over the grid.

A formal p-value for Family B is computed from the chi-squared survival function
with dof = (number of bins) − 1 = 11 degrees of freedom, where the subtraction of
one accounts for the single fitted parameter β. Family A, which generated the
target distribution, achieves χ²(A) = 0 trivially and is not assigned a p-value.
The working criterion for inclusive non-identifiability is p(B) > 0.05: if Family B
cannot be rejected at the 5% level on the inclusive projection, the two families
are operationally non-identifiable in that view.

## 6.3 Conditional Evaluation Without Refitting

After the inclusive fit, both families are evaluated in conditional subsamples
defined by pT(Z) and jet multiplicity. Crucially, the nuisance parameters are
held fixed at the values obtained from the inclusive fit; no per-stratum
refitting is performed.

This design is deliberate. Per-stratum refitting would allow each family to
optimise against the conditional distribution independently, potentially masking
genuine differences in their internal structure. By holding parameters fixed, the
stratified evaluation tests whether the inclusive best-fit parameters also produce
acceptable descriptions of the conditional distributions,
without any additional degrees of freedom.

The per-stratum discrimination is quantified by:

    Δχ²(stratum) = χ²(B, stratum) − χ²(A, stratum)

Because χ²(A) ≈ 0 in every stratum (Family A generated the target distribution at
the same parameter values used in the stratified evaluation), Δχ² ≈ χ²(B) in each
stratum. A large Δχ² indicates that Family B's best-fit parameter, which was
sufficient for the inclusive distribution, fails to reproduce the target when
conditioned on the stratification variable.

The pT(Z) stratification uses four bins: [0, 10), [10, 30), [30, 60), and
≥60 GeV. The jet-multiplicity stratification uses three categories: 0-jet, 1-jet,
and ≥2-jet, with the same jet definition as in the event selection (Section 3.4).

The conditional evaluation deliberately fixes nuisance parameters at the values
obtained from the inclusive fit rather than allowing refitting within each stratum.
This design converts the stratified comparison into a consistency test rather than
a parameter-estimation procedure. Allowing parameters to refit independently in
each stratum would permit each nuisance family to absorb conditional structure
separately, potentially masking mechanistic inconsistencies between the model and
the injected target. The fixed-parameter evaluation therefore functions analogously to a held-out
diagnostic: the question being tested is whether the inclusive best-fit parameters
remain consistent with the conditional distributions without additional degrees of
freedom. The stratified evaluation therefore constitutes a held-out consistency
test rather than a parameter-estimation step. Allowing the nuisance parameters to
float independently in each conditional subset would introduce additional degrees of
freedom and would test the flexibility of the parameterisation rather than the
internal consistency of the mechanism inferred from the inclusive distribution.

## 6.4 Working-Point Selection

The injection amplitude α₀ determines the size of the distortion injected into
the target distribution and therefore controls how far the target departs from
the nominal MET distribution. At large α₀, the distortion is large enough that
even the inclusive distribution distinguishes the two families; at small α₀, the
distortion is undetectable and neither family can be constrained.

To identify the practical degeneracy regime, the analysis was run at four injection
amplitudes: α₀ ∈ {0.02, 0.05, 0.10, 0.15}. The inclusive p-value for Family B
as a function of injection amplitude provides a threshold map:

- At α₀ = 0.02, both families fit the target with p(B) near unity.
- At α₀ = 0.05, Family B achieves χ²(B) = 3.13, dof = 11, p = 0.989.
- At α₀ = 0.10, Family B begins to fail the inclusive criterion (χ²(B) ≈ 41).
- At α₀ = 0.15, Family B is clearly rejected on the inclusive distribution.

The working point α₀ = 0.05 is selected for the main demonstration because it
lies firmly within the non-identifiable regime: at this amplitude, a 5% maximum
MET scale distortion in the highest-recoil events is invisible to the inclusive
projection despite a sample of more than 2 × 10⁵ events, yet conditional
projections recover clear discrimination. This choice represents a physically
plausible distortion amplitude — small enough to pass inclusive validation, large
enough to have interpretive consequences — rather than a worst-case or best-case
scenario.
