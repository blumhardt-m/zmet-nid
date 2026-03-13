#!/usr/bin/env python3
"""
make_landscape_figure.py — identifiability landscape scatter plot

Plots inclusive p-value vs. max per-stratum delta-chi² for the three nuisance
families, visualising the spectrum from deep degeneracy (B) through boundary
regime (C) to the injected truth (A).

Output: figures/inclusive_fit/identifiability_landscape.png
"""
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from pathlib import Path

Path("figures/inclusive_fit").mkdir(parents=True, exist_ok=True)

# Inclusive p-value and max per-stratum delta-chi² at alpha=0.05 working point
families = [
    ("A  (injected truth)",     1.000,  0.0,    "#333333", "s"),
    ("B  (topology-mixture)",   0.989,  3415.0, "#d62728", "^"),
    ("C  (unclustered energy)", 0.075,  58.0,   "#2ca02c", "o"),
]

fig, ax = plt.subplots(figsize=(8.5, 6.0))

thresh = 0.05   # p-value rejection boundary

# --- zone shading ---
# Red zone: high p-value (good inclusive fit) → deep degeneracy risk
ax.fill_betweenx([0, 4400], thresh, 1.05,
                 color="#d62728", alpha=0.08, label="_nolegend_")
# Yellow zone: near the inclusive rejection boundary
ax.fill_betweenx([0, 4400], thresh - 0.02, thresh + 0.02,
                 color="#ff7f0e", alpha=0.15, label="_nolegend_")
# Green zone: high conditional discrimination → stratification restores identifiability
ax.fill_between([0, 1.05], 800, 4400,
                color="#2ca02c", alpha=0.08, label="_nolegend_")

# Zone labels
ax.text(0.55, 120,  "inclusive\ndegeneracy",  fontsize=8, color="#d62728",
        ha="center", va="center", style="italic")
ax.text(0.05, 120,  "boundary\nregime",       fontsize=8, color="#ff7f0e",
        ha="center", va="center", style="italic")
ax.text(0.50, 3000, "conditional discrimination\nrestores identifiability",
        fontsize=8, color="#2ca02c", ha="center", va="center", style="italic")

# Rejection threshold vertical line
ax.axvline(thresh, color="gray", ls="--", lw=1.3, alpha=0.75,
           label=f"p = 0.05 inclusive rejection threshold")

# Data points
label_offsets = [(-0.04, 280), (-0.08, -480), (0.025, 220)]
for (label, pval, max_d, col, mk), (xo, yo) in zip(families, label_offsets):
    ax.scatter(pval, max_d, s=160, color=col, marker=mk, zorder=6, linewidths=1.5,
               edgecolors="white")
    ax.annotate(label, xy=(pval, max_d), xytext=(pval + xo, max_d + yo),
                fontsize=10, ha="left", color=col,
                arrowprops=dict(arrowstyle="-", color=col, lw=0.8, alpha=0.7))

ax.set_xlabel("Inclusive p-value  (higher = better inclusive fit = deeper degeneracy risk)", fontsize=11)
ax.set_ylabel("max Δχ²  across all stratifications", fontsize=11)
ax.set_title(
    "Identifiability landscape: inclusive fit quality vs. conditional discrimination\n"
    "(injection working point α = 0.05)",
    fontsize=10,
)
ax.set_xlim(-0.05, 1.08)
ax.set_ylim(-400, 4400)
ax.legend(fontsize=9, loc="upper left")

fig.tight_layout()
out = Path("figures/inclusive_fit/identifiability_landscape.png")
fig.savefig(str(out), dpi=130, bbox_inches="tight")
plt.close(fig)
print(f"Wrote: {out}")
