#!/usr/bin/env python3
"""
make_landscape_figure.py — identifiability landscape scatter plot

Plots inclusive p-value vs. max per-stratum delta-chi² for the three nuisance
families, with horizontal interpretation bands and symlog y-axis.

Output: figures/inclusive_fit/identifiability_landscape.png
"""
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from pathlib import Path

Path("figures/inclusive_fit").mkdir(parents=True, exist_ok=True)

# Data: inclusive p-value and max per-stratum delta-chi² at alpha=0.05
# columns: label, p-value, max-Dchi2, color, marker, size, edgecolor, edgewidth
families = [
    ("A  (injected truth)",     1.000,  0.0,    "#333333", "*", 120, "white", 1.0),
    ("B  (topology-mixture)",   0.989,  3415.0, "#d62728", "o", 120, "black", 1.2),
    ("C  (unclustered energy)", 0.075,  58.0,   "#2ca02c", "^",  90, "white", 1.0),
]

fig, ax = plt.subplots(figsize=(8.5, 6.0))

# --- symlog y-axis: linear below 100, log above ---
ax.set_yscale("symlog", linthresh=100)
ax.set_ylim(-8, 5000)
ax.set_yticks([0, 10, 20, 50, 100, 200, 500, 1000, 3415])
ax.yaxis.set_major_formatter(matplotlib.ticker.ScalarFormatter())

# -----------------------------------------------------------------------
# Zone shading: two-axis classification (p-value × Δχ²)
#   p_thresh = 0.05  separates  "inclusive rejection" (left) from "passes" (right)
#   dchi2_thresh = 100  separates low and high conditional discrimination
#
#   GREEN  — p > 0.05, Δχ² < 100  : inclusive degeneracy zone
#            (passes inclusive, conditional also uninformative)
#   RED    — p > 0.05, Δχ² ≥ 100  : hidden failure zone
#            (passes inclusive, but conditional reveals mechanistic inconsistency)
#   YELLOW — narrow strip near p = 0.05 boundary : marginal identifiability
# -----------------------------------------------------------------------
pval_thresh  = 0.05
dchi2_thresh = 100

# Green: lower-right quadrant
ax.fill_between([pval_thresh, 1.08], -10, dchi2_thresh,
                color="#2ca02c", alpha=0.10, zorder=0)
# Red: upper-right quadrant
ax.fill_between([pval_thresh, 1.08], dchi2_thresh, 5000,
                color="#d62728", alpha=0.09, zorder=0)
# Yellow: marginal strip bracketing the p = 0.05 boundary
ax.axvspan(0.025, 0.15, color="#ffcc00", alpha=0.10, zorder=0)

# Thin boundary lines
ax.axhline(dchi2_thresh, color="gray", linewidth=0.6, alpha=0.4)

# Zone labels
ax.text(0.60,  4.5, "Inclusive degeneracy zone",
        fontsize=9, color="#2a7a2a", va="center", style="italic")
ax.text(0.60,  700, "Hidden failure zone",
        fontsize=9, color="#b02020", va="center", style="italic", fontweight="bold")
ax.text(0.065, 4.5, "Marginal\nidentifiability",
        fontsize=8, color="#996600", va="center", style="italic")

# --- Δχ² rejection reference line (p = 0.05, dof = 11, chi2 ≈ 19.7) ---
chi2_ref = 19.7
ax.axhline(chi2_ref, linestyle="--", linewidth=0.9, color="gray", alpha=0.65)
ax.text(1.05, chi2_ref * 1.08, "p = 0.05 rejection threshold",
        fontsize=8, color="gray", ha="right", va="bottom")

# --- vertical p-value threshold line (marks inclusive rejection boundary) ---
ax.axvline(pval_thresh, color="#996600", ls=":", lw=1.2, alpha=0.70,
           label="p = 0.05 inclusive rejection threshold")

# --- data points with label annotations ---
# (xytext offsets in data coordinates, chosen to avoid overlaps on symlog axis)
annotations = [
    (1.000,  0.0,    0.84, -3),     # A: label left and below
    (0.989,  3415.0, 0.60, 1800),   # B: label left and lower (log region)
    (0.075,  58.0,   0.16, 75),     # C: label right and above
]
for (label, pval, max_d, col, mk, sz, ec, elw), (px, py, lx, ly) in zip(families, annotations):
    ax.scatter(pval, max_d, s=sz, color=col, marker=mk, zorder=6,
               edgecolors=ec, linewidths=elw)
    ax.annotate(label, xy=(pval, max_d), xytext=(lx, ly),
                fontsize=10, ha="left", color=col,
                arrowprops=dict(arrowstyle="-", color=col, lw=0.8, alpha=0.7))

# --- annotation arrow for hidden failure zone (pointing to Family B) ---
ax.annotate(
    "Inclusive validation\nsuggests agreement,\nbut stratified diagnostics\nreveal mechanistic\ninconsistency",
    xy=(0.989, 3415.0),
    xytext=(0.55, 3000),
    fontsize=8,
    color="#8b0000",
    ha="right",
    va="center",
    arrowprops=dict(arrowstyle="->", color="#8b0000", lw=1.0),
    bbox=dict(boxstyle="round,pad=0.3", fc="white", ec="#d62728", alpha=0.8),
    zorder=7,
)

# --- axes and formatting ---
ax.set_xlabel(
    "Inclusive p-value  (higher = better inclusive fit = deeper degeneracy risk)",
    fontsize=11,
)
ax.set_ylabel("max Δχ²  across all stratifications", fontsize=11)
ax.set_title(
    "Identifiability landscape: inclusive fit quality vs. conditional discrimination\n"
    "(injection working point α = 0.05)",
    fontsize=10,
)
ax.set_xlim(-0.05, 1.08)
ax.grid(alpha=0.15)
ax.legend(fontsize=9, loc="upper left")

fig.tight_layout()
out = Path("figures/inclusive_fit/identifiability_landscape.png")
fig.savefig(str(out), dpi=130, bbox_inches="tight")
plt.close(fig)
print(f"Wrote: {out}")
