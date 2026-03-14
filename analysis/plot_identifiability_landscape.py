#!/usr/bin/env python3
"""
plot_identifiability_landscape.py — identifiability landscape from CSV

Reads outputs/identifiability_landscape.csv and produces:
  outputs/identifiability_landscape.png
  outputs/identifiability_landscape.pdf
  figures/identifiability_landscape.pdf  (for manuscript inclusion)

Zones:
  Rejected models   : p ≤ 0.05 (left of threshold)
  Inclusive degeneracy : p > 0.05 AND max Δχ² < dchi2_thresh
  Hidden divergence    : p > 0.05 AND max Δχ² ≥ dchi2_thresh
"""
import csv
import shutil
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.ticker
from pathlib import Path

# ── load data ───────────────────────────────────────────────────────────────
csv_path = Path("outputs/identifiability_landscape.csv")
families = []
with open(csv_path) as f:
    reader = csv.DictReader(f)
    for row in reader:
        families.append((
            row["family"],
            float(row["inclusive_p_value"]),
            float(row["max_conditional_delta_chi2"]),
        ))

# ── layout ──────────────────────────────────────────────────────────────────
marker_map  = {"A": "*", "B": "o", "C": "^"}
color_map   = {"A": "#333333", "B": "#d62728", "C": "#2ca02c"}
size_map    = {"A": 130, "B": 120, "C": 100}
label_map   = {
    "A": "A  (injected truth)",
    "B": "B  (topology-mixture)",
    "C": "C  (unclustered energy)",
}
ann_offset  = {
    "A": (0.84,  -3),
    "B": (0.58, 1800),
    "C": (0.15,   75),
}

pval_thresh  = 0.05
dchi2_thresh = 100
chi2_ref     = 19.7          # chi2 critical value for p=0.05, dof=11

fig, ax = plt.subplots(figsize=(8.5, 6.0))

# ── symlog y-axis ────────────────────────────────────────────────────────────
ax.set_yscale("symlog", linthresh=100)
ax.set_ylim(-8, 5000)
ax.set_yticks([0, 10, 20, 50, 100, 200, 500, 1000, 3415])
ax.yaxis.set_major_formatter(matplotlib.ticker.ScalarFormatter())

# ── zone shading ─────────────────────────────────────────────────────────────
# Rejected models: left of p = 0.05
ax.axvspan(-0.05, pval_thresh, color="#aaaaaa", alpha=0.12, zorder=0)
ax.text(0.01, 700, "Rejected\nmodels",
        fontsize=8, color="#555555", va="center", style="italic")

# Inclusive degeneracy: right + below dchi2_thresh
ax.fill_between([pval_thresh, 1.08], -10, dchi2_thresh,
                color="#2ca02c", alpha=0.10, zorder=0)
ax.text(0.60, 4.5, "Inclusive degeneracy",
        fontsize=9, color="#2a7a2a", va="center", style="italic")

# Hidden divergence: right + above dchi2_thresh
ax.fill_between([pval_thresh, 1.08], dchi2_thresh, 5000,
                color="#d62728", alpha=0.09, zorder=0)
ax.text(0.60, 700, "Hidden divergence",
        fontsize=9, color="#b02020", va="center", style="italic", fontweight="bold")

# ── reference lines ──────────────────────────────────────────────────────────
ax.axvline(pval_thresh, color="#996600", ls=":", lw=1.2, alpha=0.75,
           label="p = 0.05 inclusive threshold")
ax.axhline(chi2_ref, ls="--", lw=0.9, color="gray", alpha=0.65)
ax.axhline(dchi2_thresh, color="gray", lw=0.6, alpha=0.40)
ax.text(1.05, chi2_ref * 1.08, "p = 0.05 rejection threshold",
        fontsize=8, color="gray", ha="right", va="bottom")

# ── data points ──────────────────────────────────────────────────────────────
for fam, pval, max_d in families:
    ax.scatter(pval, max_d, s=size_map[fam], color=color_map[fam],
               marker=marker_map[fam], zorder=6,
               edgecolors="black" if fam == "B" else "white", linewidths=1.0)
    lx, ly = ann_offset[fam]
    ax.annotate(label_map[fam], xy=(pval, max_d), xytext=(lx, ly),
                fontsize=10, ha="left", color=color_map[fam],
                arrowprops=dict(arrowstyle="-", color=color_map[fam],
                                lw=0.8, alpha=0.7))

# Annotation for hidden divergence (Family B)
ax.annotate(
    "Inclusive validation\nsuggests agreement,\nbut stratified diagnostics\nreveal mechanistic\ninconsistency",
    xy=(0.989, 3415.0), xytext=(0.55, 3000),
    fontsize=8, color="#8b0000", ha="right", va="center",
    arrowprops=dict(arrowstyle="->", color="#8b0000", lw=1.0),
    bbox=dict(boxstyle="round,pad=0.3", fc="white", ec="#d62728", alpha=0.85),
    zorder=7,
)

# ── axes and formatting ───────────────────────────────────────────────────────
ax.set_xlabel(
    "Inclusive p-value  (higher = better inclusive fit = deeper degeneracy risk)",
    fontsize=11)
ax.set_ylabel("max Δχ²  across all stratifications", fontsize=11)
ax.set_title(
    "Identifiability landscape: inclusive fit quality vs. conditional discrimination\n"
    "(injection working point α = 0.05)",
    fontsize=10)
ax.set_xlim(-0.05, 1.08)
ax.grid(alpha=0.15)
ax.legend(fontsize=9, loc="upper left")
fig.tight_layout()

# ── save ─────────────────────────────────────────────────────────────────────
Path("outputs").mkdir(exist_ok=True)
for ext in ("png", "pdf"):
    out = Path(f"outputs/identifiability_landscape.{ext}")
    fig.savefig(str(out), dpi=130, bbox_inches="tight")
    print(f"Wrote: {out}")

# copy PDF to figures/ for manuscript inclusion
Path("figures").mkdir(exist_ok=True)
shutil.copy("outputs/identifiability_landscape.pdf",
            "figures/identifiability_landscape.pdf")
print("Copied: figures/identifiability_landscape.pdf")

plt.close(fig)
