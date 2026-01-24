#!/usr/bin/env python3
"""
ZMET-NID Nuisance Model Families
================================

Family A: Jet-correlated MET scaling
Family B: Pileup-correlated MET broadening
"""

import numpy as np

class FamilyA:
    """Jet-correlated MET broadening."""
    
    def __init__(self, alpha=0.0, gamma=0.0, pt0=30.0):
        self.alpha = alpha
        self.gamma = gamma
        self.pt0 = pt0
    
    def scale_factor(self, n_jets, pt_lead):
        """Compute MET scale factor."""
        jet_term = self.alpha * n_jets
        lead_term = self.gamma * np.maximum(0, (pt_lead - self.pt0) / self.pt0)
        return 1.0 + jet_term + lead_term
    
    def apply(self, met, n_jets, pt_lead):
        """Apply scaling to MET."""
        return met * self.scale_factor(n_jets, pt_lead)


class FamilyB:
    """Pileup-correlated MET broadening."""
    
    def __init__(self, beta=0.0, npv0=20):
        self.beta = beta
        self.npv0 = npv0
    
    def scale_factor(self, npv):
        """Compute MET scale factor."""
        return 1.0 + self.beta * np.maximum(0, (npv - self.npv0) / self.npv0)
    
    def apply(self, met, npv):
        """Apply scaling to MET."""
        return met * self.scale_factor(npv)
