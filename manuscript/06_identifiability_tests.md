# 6. Identifiability Tests and Metrics

This section defines the criteria by which competing nuisance model families are assessed for distinguishability using CMS Open Data observables. The goal is not to identify the "true" mechanism responsible for MET tails, but to determine whether the available data representation supports unique mechanistic attribution under reasonable validation strategies.

## 6.1 Operational Definition of Identifiability

Let D denote the observed data distribution of a target observable (here, the MET magnitude), and let M_i(θ_i) denote a nuisance model family with parameters θ_i.

We say that two nuisance families M_A and M_B are **non-identifiable** with respect to observable set O if:

1. There exist parameter settings θ_A* and θ_B* such that M_A(θ_A*) ≈ D and M_B(θ_B*) ≈ D according to pre-registered goodness-of-fit criteria; and

2. The parameter settings θ_A* and θ_B* correspond to distinct causal hypotheses (e.g., jet-correlated vs. pileup-correlated mechanisms); and

3. No test constructed from observables in O rejects one family while accepting the other at the specified sensitivity level.

Identifiability is therefore a property of the model–observable pair, not of the model in isolation.

## 6.2 Data Splitting and Fit Discipline

To avoid overfitting and implicit data reuse, all analyses follow a fixed split procedure:

- The dataset is partitioned into a **training split** and an **evaluation split**, using lumisection-parity partitioning.
- Nuisance parameters θ_i are fitted only on the training split.
- All identifiability tests are evaluated exclusively on the evaluation split.

This discipline ensures that apparent agreement is not driven by trivial parameter tuning on the same data being evaluated.

## 6.3 Primary Fit Quality Metrics

### 6.3.1 Wasserstein-1 Distance

The primary metric for distributional agreement is the first Wasserstein distance W_1, defined for one-dimensional distributions P and Q as:

W_1(P, Q) = ∫₀¹ |F_P⁻¹(u) - F_Q⁻¹(u)| du

where F⁻¹ denotes the quantile function.

The Wasserstein distance is sensitive to differences across the full distribution, including the tail region, and is less dominated by localized fluctuations than pointwise metrics.

### 6.3.2 Kolmogorov–Smirnov Statistic

As a complementary metric, we compute the Kolmogorov–Smirnov (KS) statistic:

D_KS = sup_x |F_P(x) - F_Q(x)|

The KS statistic provides sensitivity to localized shape differences and is widely used in high-energy physics validation contexts.

### 6.3.3 Threshold Determination (Noise Floor)

Rather than fixing arbitrary thresholds, acceptance criteria are determined via a preregistered bootstrap procedure:

1. From the training split, draw N_boot bootstrap resamples of the observed MET distribution.
2. Compute pairwise W_1 and KS distances between bootstrap resamples.
3. Define: τ_W = median(W_1^boot), τ_KS = median(D_KS^boot)

A model is said to achieve distributional agreement if its evaluation-split distance to data satisfies:
W_1 ≤ τ_W and D_KS ≤ τ_KS

This procedure defines agreement relative to the intrinsic sampling variability of the data.

## 6.4 Structure-Sensitive Diagnostics

Global agreement in a single marginal distribution does not guarantee mechanistic distinguishability. To probe latent structure, we apply a set of pre-declared stratified tests, each evaluated independently.

Let S_k denote a stratification variable (e.g., jet multiplicity, pileup proxy, Z-boson transverse momentum, or MET azimuthal angle). For each bin b ∈ S_k, we compute:
W_1^(b), D_KS^(b)
between model predictions and data.

A nuisance family is considered **rejected by structure tests** if it fails the agreement criteria in any pre-declared stratification with sufficient statistics (minimum bin yield requirement).

## 6.5 Identifiability Outcomes

For a given pair of nuisance families (M_A, M_B), we define three mutually exclusive outcomes:

**Identifiable:**
Exactly one family satisfies global and structure-level agreement.

**Non-identifiable:**
Multiple families satisfy global agreement and no structure test rejects any of them.

**Unexplainable:**
No family satisfies global agreement, indicating missing mechanisms or insufficient observability.

These outcomes are evaluated per dataset and per split, with uncertainty estimated via resampling where applicable.

## 6.6 Interpretation Boundaries

It is emphasized that:

- Non-identifiability does not imply that the underlying detector effects are unknowable in principle.
- It implies only that, under the chosen data representation and observable set, mechanistic attribution is underdetermined.

The results therefore characterize limitations of inference under public-tier observability, rather than deficiencies in detector modeling or reconstruction.
