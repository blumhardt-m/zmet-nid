# 4. Nuisance Model Families

## 4.1 Motivation for phenomenological nuisance families

The nuisance models considered here are intentionally low-dimensional
phenomenological parameterisations. Their purpose is not to reproduce
collaboration-grade detector systematic models, but to represent distinct classes
of physical explanations for distortions in the missing transverse momentum
distribution that an external analyst might entertain when working with CMS Open
Data. The analysis evaluates whether such incompatible mechanisms can produce
indistinguishable agreement with the inclusive observable when only reduced-tier
information is available.

Two families are considered, each parameterised by a single continuous scalar.
Low dimensionality is deliberate: excessive model flexibility would trivially
explain any distributional feature, making non-identifiability uninformative. A
single-parameter family imposes a strong constraint on the allowed distortion
shape; degeneracy under such a constraint is therefore a meaningful result.

The two families are mutually incompatible in mechanism: one alters the MET
magnitude of individual events, the other alters the population weights across
event topologies. Agreement on the inclusive distribution is therefore not
tautological but reflects a genuine degeneracy of the marginal projection.

## 4.2 Fitting strategy

The identifiability test is implemented as a controlled injection. Family A is
applied to the nominal MET distribution at a known working point α₀, producing a
target distribution that represents a specific, calibrated departure from the
nominal. Family B is then independently fitted to the same target by grid search
over its parameter β. Fit quality is quantified using the χ² statistic:

    χ² = Σ_k (h_k^model − h_k^target)² / (h_k^target + 1)

where h_k denotes the bin count in histogram bin k, and the (obs + 1) denominator
provides stability in bins with low target occupancy. Degrees of freedom for
Family B are n_bins − 1 fitted parameter.

Inclusive fit quality determines the working point. For a given α₀, if Family B
achieves an acceptable p-value on the inclusive χ², the two families are
non-identifiable on the inclusive projection. Stratified comparisons are then
applied at the inclusive best-fit parameters, without per-stratum refitting, to
test whether conditional projections restore discrimination.

The simplified χ² goodness-of-fit criterion should be understood as a diagnostic
tool for comparing nuisance mechanisms rather than as a full profile-likelihood
inference in the standard LHC sense [12]. The goal is not optimal parameter
estimation but the detection of inclusive degeneracy and its resolution through
conditional projections.

## 4.3 Family A — Recoil-Response Distortion

Family A represents the hypothesis that MET tails arise from a miscalibrated
hadronic recoil response. In Z→μ⁺μ⁻ events, the Z boson recoils against hadronic
activity, and the magnitude of the MET reflects the vector imbalance of that
recoil. A fractional miscalibration in the hadronic response would shift the
reconstructed MET proportionally to the recoil magnitude, which is directly
approximated by pT(Z), the transverse momentum of the dimuon system.

The model modifies the event-level MET as:

    MET'_i = MET_i × (1 + α × clip(pT(Z)_i, 0, 100) / 100)

where pT(Z)_i is the dimuon transverse momentum in GeV, clip(·, 0, 100) caps the
effective recoil scale at 100 GeV to prevent unbounded scaling in the high-recoil
tail, and α is the distortion amplitude. Positive α inflates MET; negative α
deflates it.

The parameter search range is α ∈ [−0.30, +0.30] in steps of 0.01. The injection
working point used throughout this analysis is α₀ = 0.05, corresponding to a 5%
maximum MET scale distortion in the highest-recoil events, scaling linearly with
pT(Z) below 100 GeV. Family A leaves the relative population of jet-multiplicity
categories unchanged.

## 4.4 Family B — Topology-Mixture Reweighting

Family B represents the hypothesis that MET tails arise from a shift in the
event-topology composition of the sample rather than from a miscalibration of the
MET scale within any given topology. Such distortions may arise from mismodelling
of jet production rates, residual efficiency corrections that vary across jet
categories, or other effects that alter the effective mixture of zero-jet, one-jet,
and multi-jet events contributing to the inclusive distribution.

The model applies per-event weights according to jet multiplicity:

| Jet category | Weight |
|---|---|
| 0 jets | 1 |
| 1 jet  | 1 + β |
| ≥2 jets | 1 + 2β |

Jet multiplicity is counted using the same cleaned jet definition as in the event
selection (pT > 30 GeV, ΔR > 0.4 from both Z-candidate muons; Section 3.4).
The weights are normalised so that the total weighted event count equals the
unweighted count, preserving the inclusive sample size.

The parameter β controls the degree to which jet-rich topologies are upweighted
relative to zero-jet events. The search range is β ∈ [−0.40, +0.40] in steps of
0.01. At β > 0, events with more jets are upweighted, shifting the inclusive MET
distribution toward higher values because multi-jet events have systematically
larger MET than zero-jet events.

Crucially, Family B does not alter the MET value of any individual event. It
changes only the effective population of each jet-multiplicity stratum.

The nuisance families are chosen to represent qualitatively distinct physical mechanisms: topology mixture modifies the composition of event kinematics entering the MET distribution, while recoil-scale perturbations model detector-level response distortions. Demonstrating degeneracy across these mechanism classes indicates that the indistinguishability arises from observable compression rather than from trivial functional parameterizations. The nuisance families correspond to perturbations arising at different stages of the event reconstruction chain, representing distinct physical mechanisms rather than alternative parameterizations of a single model.

## 4.5 Mechanistic Contrast

The two families modify the inclusive MET distribution through fundamentally
different mechanisms:

| | Family A | Family B |
|---|---|---|
| **Mechanism** | Recoil-response distortion | Topology-mixture reweighting |
| **Parameter** | α | β |
| **Effect on event-level MET** | Scales MET with pT(Z) | None |
| **Effect on jet-category populations** | None | Upweights multi-jet events |
| **Injection working point** | α = 0.05 | β = +0.290 (best fit) |

Family A is sensitive to the recoil–MET relationship and is exposed by
conditioning on pT(Z). Family B is sensitive to population consistency across jet
categories and is exposed by conditioning on jet multiplicity. Neither diagnostic
is implicit in the inclusive MET distribution, which is what makes the inclusive
degeneracy possible.

## 4.6 Scope of Nuisance Parameterisations

The nuisance models used in this study are intentionally simplified representations
of broader classes of detector and modelling effects. The goal is not to reproduce
collaboration-grade systematic uncertainties, which require detailed detector
simulation and calibration studies, but to construct phenomenological families that
capture distinct mechanisms by which the inclusive MET distribution may be
distorted. The analysis therefore tests identifiability at the level of mechanism
rather than at the level of detailed detector modelling.

If two mechanistically distinct nuisance explanations cannot be distinguished in
the inclusive projection even in this simplified setting, the ambiguity will
persist — and likely worsen — for more realistic high-dimensional systematic
models. The simplified parameterisations therefore provide a conservative lower
bound on the identifiability challenge: demonstrating degeneracy here implies that
more complex nuisance families face at least as severe an identifiability
constraint. Because the nuisance families are restricted to a single parameter,
the degeneracy demonstrated in the inclusive projection represents a conservative
lower bound on the identifiability challenge. Increasing the dimensionality of
the nuisance parameterisation without introducing additional observables generally
enlarges the space of distributions compatible with the data rather than reducing it.
