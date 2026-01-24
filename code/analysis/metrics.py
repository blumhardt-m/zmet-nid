#!/usr/bin/env python3
"""
ZMET-NID Fit Quality Metrics
============================

Implements W1, KS, and bootstrap threshold determination.
"""

import numpy as np
from scipy import stats

def wasserstein_1d(x, y):
    """Compute 1D Wasserstein distance."""
    return stats.wasserstein_distance(x, y)

def ks_statistic(x, y):
    """Compute Kolmogorov-Smirnov statistic."""
    stat, _ = stats.ks_2samp(x, y)
    return stat

def bootstrap_thresholds(data, n_boot=100, seed=42):
    """
    Determine noise-floor thresholds via bootstrap.
    
    Returns:
        tau_w: W1 threshold
        tau_ks: KS threshold
    """
    rng = np.random.default_rng(seed)
    n = len(data)
    
    w1_distances = []
    ks_distances = []
    
    for _ in range(n_boot):
        idx1 = rng.choice(n, size=n, replace=True)
        idx2 = rng.choice(n, size=n, replace=True)
        
        sample1 = data[idx1]
        sample2 = data[idx2]
        
        w1_distances.append(wasserstein_1d(sample1, sample2))
        ks_distances.append(ks_statistic(sample1, sample2))
    
    tau_w = np.median(w1_distances)
    tau_ks = np.median(ks_distances)
    
    return tau_w, tau_ks
