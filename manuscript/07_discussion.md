# 7. Discussion

## 7.1 Inclusive degeneracy and reduced observability

As demonstrated in Section 6, the inclusive MET distribution admits both nuisance families despite the large event sample. The topology-mixture model achieves χ²(B) = 3.13 on 11 degrees of freedom (p = 0.989), a fit quality statistically indistinguishable from the injected recoil-response model that generated the target distribution.

This result indicates a practical non-identifiability of nuisance explanations in the inclusive projection: the observed data are consistent with more than one generative model, and no hypothesis test on the inclusive projection can select among them. The non-identifiability here is not a consequence of insufficient data — 233,524 events is a large sample by typical analysis standards — but of insufficient observable dimensionality relative to the complexity of the nuisance space being probed.

The explanation lies in the nature of the inclusive projection. Both Family A and Family B alter the single-observable distribution P(MET) through mechanistically different operations. Family A applies a continuous MET scale factor that grows with pT(Z), shifting the MET distribution rightward in a recoil-correlated manner. Family B applies global per-event weights that shift the event-topology mixture, altering the effective composition of the sample without altering the MET scale within any given topology. At a sufficiently small distortion amplitude, the integrated shape changes produced by these two mechanisms can be made to agree to within statistical noise on the marginal distribution. The existence of a best-fit β that achieves this agreement — at a value of β = +0.290, well inside the search range — confirms that the degeneracy is not accidental but reflects the limited discriminating power of the inclusive one-dimensional projection.

This generalises a known challenge in MET analyses: inclusive validation distributes sensitivity across the entire observable space, potentially concentrating statistical power where signal-to-background from nuisance effects is weakest. The present study provides a quantitative illustration of this effect at a controlled and physically motivated working point. While the empirical example studied here involves missing transverse energy, the mechanism arises from observable compression and therefore applies to any collider analysis relying on reduced summary observables.

## 7.2 Stratified diagnostics as identifiability tests

The stratified analyses demonstrate that the non-identifiability of the inclusive distribution does not reflect a fundamental indistinguishability of the two nuisance models. Conditional projections P(MET | S) expose clear model differences that are invisible to the inclusive comparison P(MET).

The pT(Z) stratification tests the recoil–MET relationship. Family A is constructed so that its MET distortion is proportional to pT(Z): events with large pT(Z) receive large MET shifts, and events with pT(Z) ≈ 0 receive negligible shifts. When the analysis is conditioned on pT(Z), this design is directly tested — Family A must reproduce the target within each pT(Z) stratum using the same parameter fit to the inclusive distribution, and it does so with χ²(A) ≈ 0 in all bins. Family B, which does not scale its distortion with pT(Z), cannot reproduce the conditional distributions while maintaining the parameter values that were optimised for the inclusive projection. The result is Δχ² values ranging from 95 in the 10–30 GeV stratum to 1,432 in the ≥60 GeV stratum, a discrimination of three orders of magnitude relative to the inclusive χ²(B) = 3.13.

The jet-multiplicity stratification tests a different latent assumption: population consistency across event topologies. Family B assigns different event weights to the three jet categories, which globally alters the relative event counts within each category. When the distribution is conditioned on jet multiplicity, this population shift is directly exposed: Family B produces approximately 44% more events in the ≥2-jet stratum than the target, a count inflation that generates Δχ² = 3,415 in that bin. Family A, which does not alter jet-category populations, achieves χ²(A) ≈ 0 throughout.

The two stratifications therefore probe different latent assumptions of the alternative model: the pT(Z) test examines the recoil–MET relationship, while the jet-multiplicity test examines population consistency across event topologies. The pT(Z) test is sensitive to the functional dependence of MET on recoil activity; the jet test is sensitive to population consistency across topologies. A nuisance model that passes both conditioning tests would need to jointly preserve the recoil–MET relationship and the jet-category population fractions — a substantially more constrained requirement than agreement on the marginal MET distribution alone. The practical implication is that a two-dimensional stratification of this type can serve as a minimal identifiability test for inclusive nuisance estimates, without requiring access to additional observables beyond those already provided at the NanoAOD tier.

Considering the three nuisance families together reveals that identifiability in
reduced-observability analyses is graded rather than binary. The topology-mixture
family produces a deep inclusive degeneracy that is strongly exposed by
stratification, whereas the unclustered-energy perturbation occupies a boundary
regime in which the inclusive projection is only marginally discriminating and
conditional diagnostics produce moderate discrepancies. The efficiency with which
stratification restores discrimination depends on which latent structural assumption
of the nuisance family is violated.

## 7.3 Implications for analyses using compressed public datasets

Public analysis formats such as NanoAOD provide a reduced set of reconstructed observables and omit low-level detector quantities, effectively compressing the dimensionality of the observable event space. The CMS NanoAOD format in particular reduces per-event storage by retaining only high-level reconstructed quantities, discarding lower-level detector information. For external analyses working exclusively with public NanoAOD, the observable space available for nuisance characterisation is materially narrower than that available in a full AOD or MiniAOD analysis. This reduction has a direct implication for identifiability: as the dimensionality of the observable space decreases, the probability that P(X) is degenerate across distinct nuisance models increases for any single observable X.

The families studied here were chosen to be simple enough to fit and interpret with public data, but the degeneracy they exhibit is not a pathology of the specific choice. Any two nuisance models that produce compatible marginal distributions will be inclusive-indistinguishable on the compressed observable, regardless of their physical motivation. The analysis therefore points toward a systematic issue rather than an artifact: nuisance estimation from inclusive comparisons on compressed data may be subject to identifiability failures that remain invisible without explicit diagnostic checks.

Within collaboration analyses, Z→ll control samples are routinely used to constrain
recoil-related nuisances in precision measurements such as the W boson mass [15];
the question addressed here is how much of that discriminating power remains
available at the NanoAOD tier to an external analyst without full
systematic-variation access.

The conditional projections used here represent a minimal response to this issue, requiring only the stratification variables pT(Z) and jet multiplicity that are already present in NanoAOD. The approach requires no additional data access and imposes no significant computational overhead beyond the initial inclusive fit. It does require that the analyst commit to stratification variables in advance, before examining the stratified distributions, to avoid post-hoc selection of the stratification that best supports a preferred model — a standard requirement in confirmatory analysis that is made explicit in the frozen protocol accompanying this study.

## 7.4 Limitations and scope of the demonstration

Several limitations bound the conclusions that can be drawn from this study.

First, only two nuisance families are compared. The demonstration establishes that at least one topology-based explanation exists for the inclusive MET distribution at the working point α = 0.05, but it does not characterise the full space of inclusive-indistinguishable models. A broader family survey would be required to map the degeneracy manifold and understand which model features are and are not constrained by the inclusive projection.

Second, the analysis is restricted to a single decay channel (Z→μ⁺μ⁻) and a single observable (MET). The qualitative finding — that conditional projections restore discrimination absent from the marginal — is expected to hold more broadly, but the specific Δχ² values and the working point at which degeneracy arises will depend on the channel, the nuisance families, and the observable.

Third, the nuisance families were constructed to be analytically tractable, not to represent the most realistic models of detector effects. Family A's linear recoil response and Family B's jet-count reweighting are simplifications of more complex physical processes. A study using more realistic nuisance families would require a more sophisticated fitting procedure and potentially a larger range of stratification variables to expose the relevant failure modes.

Fourth, the inclusive degeneracy demonstrated here is established at a specific working point (α = 0.05, β = +0.290). The threshold sweep described in the supplementary material shows that inclusive discrimination is restored for α ≥ 0.10. The present result therefore characterises a regime of modest distortion amplitudes; whether similar degeneracies arise at other working points, or under different nuisance parameterisations, requires separate investigation.

These limitations notwithstanding, the central finding is robust: inclusive MET in Z→μ⁺μ⁻ events from CMS Open Data admits a topology-mixture explanation at a physically plausible distortion amplitude, and two targeted conditional projections suffice to expose the failure of that explanation. The framework demonstrated here — inclusive fit, conditional evaluation, Δχ² quantification — provides a practical template for identifiability diagnostics in compressed-data analyses. Within these limitations, the example provides a concrete demonstration that inclusive validation tests alone may be insufficient to guarantee interpretability in reduced-observability collider analyses.

A supplementary test using a physically grounded nuisance variation derived from
CMS-exposed unclustered-energy systematic branches further clarifies the structure
of the degeneracy. A scaled version of the unclustered-energy shift is not rejected
by the inclusive test at the 5% level (χ² = 18.28, 11 dof, p = 0.075), but
produces noticeably poorer agreement than the topology-mixture model and yields
only modest discrepancies under conditional stratification. This comparison
suggests that non-identifiability in reduced-tier analyses is mechanism-dependent.
Some nuisance families — particularly those altering event-topology mixtures —
can produce deep inclusive degeneracies, whereas detector-level perturbations that
preserve event topology remain more readily recoverable once conditional structure
is examined.

The degeneracies demonstrated here arise from the restricted observable space of
compressed public data tiers such as NanoAOD; they do not imply limitations in
collaboration-internal analyses, which typically exploit higher-dimensional
observables and systematic variation branches unavailable in the public format.
