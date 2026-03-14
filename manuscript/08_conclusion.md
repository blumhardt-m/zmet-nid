# 8. Conclusions

A total of 233,524 Z→μ⁺μ⁻ candidate events were selected from CMS Open Data
and used to compare two physically distinct nuisance families on the inclusive
missing transverse momentum distribution. The topology-mixture family (Family B),
fitted by grid search to the same target as the injected recoil-response family
(Family A), achieves χ²(B) = 3.13 on 11 degrees of freedom (p = 0.989). The
inclusive MET distribution does not reject the topology-mixture nuisance family
at the working point α = 0.05, despite more than 2 × 10⁵ selected events.

Stratified analyses applied without refitting recover the discriminating
information absent from the inclusive projection. Conditioning on pT(Z) yields
Δχ² = χ²(B) − χ²(A) = 1,432 in the highest pT(Z) bin (pT(Z) ≥ 60 GeV),
with values growing monotonically with recoil activity. Conditioning on jet
multiplicity yields Δχ² = 3,415 in the ≥2-jet bin, reflecting the failure of
the topology-mixture model to preserve the conditional event population across
jet categories. The two stratifications probe different latent assumptions of
the nuisance model — the recoil–MET relationship and population consistency
across event topologies — and each exposes discrepancies absent from the
inclusive view.

A supplementary test using the CMS unclustered-energy systematic variation shows
that physically grounded perturbations can also approach inclusive degeneracy,
though with weaker conditional discrepancies than the topology-mixture family.
The comparison of three distinct nuisance families demonstrates that
reduced-tier identifiability is not an all-or-nothing property but a spectrum
determined by how the underlying mechanism projects onto the observable space.

The general principle illustrated is that inclusive validation constrains the
marginal distribution P(MET), while stratified diagnostics probe conditional
distributions P(MET | S), where S is a physics-motivated stratification
variable. When analyses are performed on compressed data formats such as
NanoAOD, the reduced dimensionality of the observable space increases the
probability that P(MET) is degenerate across distinct nuisance models.
Conditional projections provide a minimal mechanism to restore partial
identifiability using the variables already available at the NanoAOD tier,
without requiring access to lower-level detector information.

Identifiability diagnostics of the type demonstrated here should accompany
inclusive validation tests in analyses that rely on inclusive nuisance
estimation from compressed public datasets. The approach extends beyond the
MET observable and the specific families studied: any inclusive comparison
that constrains only P(X) for some observable X may be subject to similar
degeneracies, and a targeted set of conditional projections P(X | S) provides
a practical diagnostic framework. The present analysis is illustrative rather
than exhaustive; a complete treatment would require a broader family of
nuisance models and a systematic study of stratification variables. The
central finding — that inclusive MET in Z→μ⁺μ⁻ events admits multiple
nuisance explanations at physically plausible distortion amplitudes, while
stratified projections distinguish them — demonstrates the interpretability
risk in reduced-observability analyses and points toward a straightforward
mitigation strategy. The use of real CMS Open Data events ensures that the
observed degeneracy arises from properties of reduced collider observables
themselves rather than artifacts of a particular simulation framework.

All analysis scripts and the frozen protocol used in this study are available
in the accompanying repository, tagged at the version corresponding to this
submission.
