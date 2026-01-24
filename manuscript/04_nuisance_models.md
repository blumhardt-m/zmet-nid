# 4. Nuisance Model Families

The objective of this section is to define a small set of competing nuisance model families that represent distinct, physically plausible explanations for missing transverse momentum (MET) tails, as could reasonably be entertained by an external analyst using CMS Open Data at the NanoAOD tier. These models are not intended to reproduce collaboration-grade systematic uncertainty procedures, nor to provide a complete description of detector effects. Rather, they are designed to probe identifiability: whether different causal explanations can be distinguished given the available observables.

## 4.1 Design Principles

All nuisance models introduced in this work satisfy the following constraints:

**Public-tier implementability:**
Models operate only on observables available in NanoAOD and do not require recomputation of MET from lower-level constituents.

**Low dimensionality:**
Each model is parameterized by a small number of free parameters, to avoid overfitting and to ensure that any observed degeneracy is not an artifact of excessive flexibility.

**Physical plausibility:**
Model dependencies are motivated by well-established correlations between MET behavior and event-level quantities such as jet activity or pileup.

**Mutual incompatibility:**
The model families encode qualitatively different mechanisms, such that agreement between them is not tautological.

## 4.2 Reference Observable

Let E_T^miss denote the reconstructed MET magnitude provided in NanoAOD. All nuisance models act multiplicatively on this quantity in simulation:

E_T,model^miss = E_T^miss × f(x; θ)

where x denotes a set of available event-level observables and θ denotes the nuisance parameters.

The multiplicative form is chosen to preserve the qualitative shape of the MET distribution while allowing controlled broadening of the tail.

## 4.3 Family A: Jet-Correlated MET Broadening

The first nuisance family represents the hypothesis that MET tails are primarily correlated with hadronic activity, reflecting jet energy mismeasurement, rare jet reconstruction pathologies, or residual jet-related effects not fully captured in simulation.

We parameterize this family as:

f_A(x; θ_A) = 1 + α·N_jet + γ·max(0, (p_T^lead - p_T^0)/p_T^0)

where:
- N_jet is the number of reconstructed jets above a fixed transverse momentum threshold
- p_T^lead is the transverse momentum of the leading jet
- p_T^0 is a fixed reference scale (chosen here as 30 GeV)
- θ_A = {α, γ} are free parameters

This form captures two intuitive effects: a global scaling with jet multiplicity and an additional contribution from particularly energetic leading jets. No attempt is made to model jet response at the constituent level, which is inaccessible in NanoAOD.

## 4.4 Family B: Pileup-Correlated MET Broadening

The second nuisance family represents the hypothesis that MET tails are dominated by pileup-related effects, such as imperfect subtraction of diffuse energy from additional proton–proton interactions.

This family is parameterized as:

f_B(x; θ_B) = 1 + β·max(0, (N_PV - N_PV^0)/N_PV^0)

where:
- N_PV is the number of reconstructed primary vertices
- N_PV^0 is a fixed reference value
- θ_B = {β} is the nuisance parameter

This model encodes a monotonic dependence of MET broadening on pileup activity without introducing additional shape degrees of freedom.

## 4.5 Scope and Interpretation

Both nuisance families are deliberately phenomenological. They do not correspond to official CMS systematic uncertainty variations, nor are they intended to provide a faithful detector simulation. Their role is instead comparative: to test whether distinct, physically motivated explanations can be tuned to reproduce observed MET distributions while implying different causal narratives.

If these families can be distinguished using public-tier observables, this indicates identifiability under the given representation. Conversely, if they remain degenerate under standard validation strategies, this indicates a limitation on mechanistic attribution imposed by reduced observability.
