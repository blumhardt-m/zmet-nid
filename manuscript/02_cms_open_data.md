# 2. CMS Open Data and Analysis Scope

The analysis presented in this work is based exclusively on publicly released CMS Open Data corresponding to proton–proton collisions recorded in 2016 at a center-of-mass energy of 13 TeV by the CMS detector. The CMS Open Data program provides external access to research-grade datasets, associated Monte Carlo simulations, and documentation intended to enable reproducible analyses outside the collaboration.

## 2.1 Open Data Tiers and Design Philosophy

CMS Open Data are released in multiple data tiers that differ in information content, file size, and intended use. For Run 2 data, the NanoAOD format represents the most compact tier designed for physics analysis, containing high-level reconstructed objects and derived quantities in a flat ntuple structure. NanoAOD is explicitly optimized for accessibility and computational efficiency, enabling analyses to be performed without reliance on the full CMS software stack or detector-level reconstruction code.

This design necessarily entails a loss of low-level information. In particular, NanoAOD does not provide access to individual particle-flow candidates, detailed reconstruction intermediates, or the full sequence of corrections and calibrations applied internally within CMS. As a result, certain detector effects—especially those associated with unclustered energy, rare reconstruction pathologies, or detailed pileup mitigation algorithms—cannot be independently recomputed or varied by external analysts.

The present work explicitly embraces these constraints. Rather than treating reduced observability as a limitation to be overcome, such compression is treated as a defining feature of the inference problem under study.

## 2.2 MET Representation in NanoAOD

In NanoAOD, missing transverse momentum is provided as a reconstructed object (`MET_pt`, `MET_phi`) that reflects the final CMS MET definition after standard corrections and pileup mitigation procedures. While jet kinematics, muon kinematics, and basic pileup proxies (e.g., number of reconstructed primary vertices) are available, the internal composition of MET—including contributions from unclustered energy—is not exposed.

Consequently, MET in NanoAOD should be understood as an end-product observable rather than a quantity that can be re-derived from constituent objects. This has important implications for causal interpretation: different detector- or reconstruction-level mechanisms may map onto similar MET distributions without being distinguishable at the public tier.

The analysis in this work does not attempt to reverse-engineer or reinterpret the CMS MET reconstruction. Instead, MET is treated as an observable whose distributional behavior is to be compared between data and simulation under phenomenological perturbations applied at the level of available observables.

## 2.3 Scope of Simulation Use

Monte Carlo simulation samples corresponding to Z+jets processes are used as a reference model for comparison with data. Consistent with the public-data context, only the simulation products provided through the CMS Open Data portal are used, without access to internal systematic variations or detector recalibration tools.

Nuisance models introduced in this work are therefore applied at the level of the simulated MET observable or via correlations with other available high-level quantities. These models are not intended to reproduce CMS systematic uncertainty procedures, but to represent distinct, physically plausible explanations that an external analyst could reasonably entertain given the accessible information.

## 2.4 Deliberate Scope Limitations

For clarity, the following scope restrictions apply:

- This work does not attempt to validate or supersede CMS MET performance studies.
- No detector-level recalibration or re-reconstruction is performed.
- No claims are made regarding optimal MET modeling or correction strategies.
- All conclusions are conditional on the information available at the NanoAOD tier.

Within this deliberately constrained scope, the analysis addresses a narrower but well-defined question: to what extent can distinct instrumental explanations for MET tails be discriminated using public-tier observables alone?
