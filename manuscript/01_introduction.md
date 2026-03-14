# 1. Introduction

Analyses performed on compressed or reduced-tier datasets face an intrinsic statistical challenge: distinct physical mechanisms can produce indistinguishable predictions for marginal observables when latent structure is integrated out. In such settings, agreement between model and data in a one-dimensional validation distribution does not guarantee that the underlying nuisance mechanism is uniquely identified. This problem becomes particularly relevant for external analyses using publicly released collider datasets, where the available observables are intentionally compressed relative to the full detector reconstruction. This work addresses a general inference problem in collider analyses: when high-dimensional event records are compressed into reduced observables, distinct nuisance mechanisms can become indistinguishable in standard validation distributions.

Collider analyses typically operate on reduced observables obtained from high-dimensional event records. Denoting the full reconstructed event record as X and the reduced observable as Y = R(X), the projection R discards information present in X that is not retained in Y. Because Y arises from such a projection, it is not generally a sufficient statistic for the underlying nuisance mechanisms: distinct values of a nuisance parameter θ can produce the same distribution P(Y | θ) even when P(X | θ) differs. Consequently, distinct physical mechanisms can produce indistinguishable inclusive distributions of a reduced observable even when their full-event descriptions differ.

Many collider measurements rely on reduced analysis tiers, where high-dimensional detector information is compressed into summary observables such as missing transverse momentum (MET). While this reduction enables efficient analysis and public data release, it can introduce an identifiability limitation: distinct physical or instrumental mechanisms may produce indistinguishable distributions once events are projected into a compressed observable space. In statistical terms, the reduction R(X) → Y can induce degeneracy where P(Y | M₁) ≈ P(Y | M₂) despite P(X | M₁) ≠ P(X | M₂). However, partial recovery of identifiability may still be possible when conditioning on structural variables that survive the reduction, such as hadronic recoil or jet multiplicity. This work investigates such compression-induced degeneracy in missing transverse momentum inference using CMS Open Data and demonstrates that mechanisms indistinguishable in inclusive MET distributions become separable once conditional structure is examined. Such effects are particularly relevant for analyses performed on compressed public data tiers, where access to full detector-level information is unavailable and inference must proceed using reduced observables alone.

In the broader statistics literature this situation corresponds to an identifiability problem, in which distinct latent mechanisms produce indistinguishable predictions once high-dimensional structure is compressed into low-dimensional summary observables. Modern collider analyses normally mitigate such degeneracies using simultaneous multi-bin profile likelihood fits across orthogonal control regions, but these techniques rely on access to the full observable space and collaboration-internal systematic variations that are not available in reduced-tier public datasets.

Missing transverse momentum (MET) provides a particularly useful observable for examining this issue. MET is a highly composite quantity, constructed from the vector sum of reconstructed particle-flow candidates and subject to corrections that depend on jet calibration, pileup mitigation, and detector response. Its high-value tail is often used as a diagnostic of detector performance. At the same time, the NanoAOD format used in CMS Open Data exposes MET primarily as a reconstructed end-product rather than as a quantity derived from constituent detector information. This makes MET an ideal setting for studying how reduced observability affects the identifiability of nuisance mechanisms. Within the CMS Collaboration, MET reconstruction and performance have been studied in detail using full detector information, dedicated calibration constants, and systematic uncertainty variations, making extensive use of control samples such as Z→ll and γ+jet events to characterize MET scale, resolution, and anomalous behavior. This study does not evaluate the CMS reconstruction pipeline itself; rather, it examines inference limits that arise when analyses operate on reduced public data products such as NanoAOD.

In parallel, CMS has released a substantial fraction of its LHC data to the public, enabling external researchers to perform analyses using well-defined reduced data tiers such as NanoAOD. These releases have supported a growing body of peer-reviewed work, including studies of jet substructure, dimuon spectra, and methodological analyses based on public data. However, the reduced-tier formats provided for open data necessarily compress detector-level information and expose high-level observables, such as MET, primarily as reconstructed end-products. While this design greatly improves accessibility, it also limits the range of detector and reconstruction effects that can be explicitly modeled or varied by external analysts.

As a consequence, external analyses using CMS Open Data face a distinct interpretive challenge: even when agreement between data and simulation is observed for standard validation distributions, it is not always clear whether the underlying causal mechanisms responsible for that agreement are uniquely determined by the available observables. In other words, agreement does not automatically imply identifiability. This distinction is particularly relevant for MET, where multiple detector-related effects can plausibly contribute to similar distributional features.

The possibility that distinct nuisance mechanisms can produce indistinguishable
predictions for marginal observables is well recognised in collider statistics.
Contemporary LHC analyses typically address such degeneracies using multidimensional
control regions and simultaneous profile-likelihood fits across correlated
observables, often supplemented by auxiliary calibration measurements [12]. More
recently, simulation-based inference and related likelihood-learning techniques have
explored similar identifiability questions in high-dimensional observable spaces
[13, 14]. These approaches implicitly rely on access to rich detector-level
information and extensive sets of systematic-variation branches available within
collaboration analysis frameworks. By contrast, external analyses using CMS Open
Data operate on compressed representations such as NanoAOD, where many of the
underlying degrees of freedom are no longer directly accessible. The present work
therefore addresses a complementary question: to what extent competing nuisance
explanations remain distinguishable when inference is restricted to the reduced
observable space available at the public data tier.

Reduced collider observables compress high-dimensional event structure into a limited set of summary statistics; this compression can erase the information required to distinguish physically distinct nuisance mechanisms, creating non-identifiability in inclusive validation tests.

The present study addresses this issue by focusing on Z→μ⁺μ⁻ events, which provide a well-established standard candle for MET studies. In these events, genuine MET is expected to be negligible, such that the observed MET distribution is dominated by instrumental effects. The analysis is restricted to publicly available NanoAOD-level observables, explicitly adopting the perspective of an external analyst operating within the constraints of CMS Open Data.

Importantly, the degeneracy demonstrated here is constructed using real collider events from CMS Open Data rather than purely simulated examples. This allows the analysis to expose inference structure that arises in practical collider datasets rather than in idealized Monte Carlo studies.

Rather than attempting to reproduce collaboration-grade MET systematics or to propose improvements to MET reconstruction, the objective is to assess the discriminating power of public-tier observables with respect to competing, physically plausible nuisance explanations for MET tails. Two phenomenological nuisance families are examined, representing distinct mechanisms by which the inclusive MET distribution may be distorted: a recoil-response scale distortion correlated with pT(Z), and a topology-mixture reweighting that alters the relative population of jet-multiplicity categories. The two nuisance families correspond to perturbations arising at different stages of the reconstruction pipeline, representing distinct physical mechanisms rather than alternative parameterizations of a single model. These mechanisms can produce similar distortions in the marginal MET distribution but imply different conditional structure when the data are stratified by recoil or by event topology.

Using a controlled injection framework and a χ² fit criterion, the analysis evaluates whether these nuisance families can be tuned to reproduce the same target MET distribution and whether the inclusive projection is sufficient to discriminate between them. The effectiveness of stratified diagnostics — conditioned on pT(Z) and jet multiplicity, both available in NanoAOD — is further examined for their ability to break the inclusive degeneracy.

By quantifying when and how such degeneracies persist, this work aims to delineate an inference boundary relevant to external analyses of CMS Open Data. The results are not intended to challenge the validity of CMS MET performance studies, but rather to clarify what kinds of mechanistic conclusions are and are not supported by reduced-tier public representations. More broadly, the approach provides a template for assessing identifiability limits of high-level observables in open collider datasets. The goal of this study is therefore not to propose an improved model of MET, but to map how different nuisance mechanisms project onto compressed observables and to identify simple diagnostics that restore mechanistic discrimination.

Figure \ref{fig:identifiability_landscape} summarizes the central argument of this paper. Distinct nuisance mechanisms can reproduce indistinguishable inclusive distributions of a reduced observable, while conditional diagnostics reveal statistically significant divergence.

\begin{figure}[t]
\centering
\includegraphics[width=0.8\linewidth]{figures/identifiability_landscape.pdf}
\caption{
Identifiability landscape for nuisance mechanism families studied in this work.
The horizontal axis shows inclusive goodness-of-fit to the injected target MET
distribution; higher values indicate stronger inclusive agreement and deeper
degeneracy risk. The vertical axis shows the maximum conditional discrepancy
across recoil and jet-multiplicity strata (symlog scale).
Configurations in the upper-right region pass inclusive validation while
producing large conditional discrepancies, illustrating how compression of
collider events into reduced observables can mask mechanistic differences
between nuisance models. Quantitative results are presented in Section 6.
}
\label{fig:identifiability_landscape}
\end{figure}
