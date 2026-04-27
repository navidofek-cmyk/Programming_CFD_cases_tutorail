#!/usr/bin/env python3
"""
Plot Cp distribution and compare with AGARD Cook et al. 1979 experimental data.
RAE 2822 Case 9: M=0.73, AoA=2.79°, Re=6.5e6
Our simulation:   M=0.729, AoA=2.31°, Re=6.5e6

Usage:
    python3 scripts/plot_cp.py [output/cp.dat]
"""
import sys
import numpy as np
import matplotlib.pyplot as plt
import os

# AGARD RAE 2822 Case 9 experimental data (Cook et al. 1979)
# Columns: x/c, Cp_upper, Cp_lower
AGARD_CASE9 = np.array([
    [0.000, 1.005, 1.005],
    [0.005, 0.069, 0.690],
    [0.015, -0.385, 0.629],
    [0.025, -0.582, 0.563],
    [0.050, -0.764, 0.456],
    [0.075, -0.844, 0.385],
    [0.100, -0.897, 0.333],
    [0.125, -0.931, 0.288],
    [0.150, -0.951, 0.249],
    [0.175, -0.970, 0.216],
    [0.200, -0.972, 0.188],
    [0.250, -0.990, 0.137],
    [0.300, -0.999, 0.091],
    [0.350, -1.005, 0.051],
    [0.400, -1.008, 0.016],
    [0.450, -1.017, -0.013],
    [0.500, -1.023, -0.036],
    [0.550, -1.007, -0.055],
    [0.575, -0.987, -0.062],
    [0.600, -0.920, -0.067],
    [0.625, -0.710, -0.070],
    [0.650, -0.432, -0.071],
    [0.675, -0.290, -0.070],
    [0.700, -0.223, -0.067],
    [0.725, -0.177, -0.062],
    [0.750, -0.148, -0.057],
    [0.775, -0.126, -0.051],
    [0.800, -0.110, -0.044],
    [0.825, -0.093, -0.037],
    [0.850, -0.079, -0.030],
    [0.875, -0.065, -0.023],
    [0.900, -0.053, -0.016],
    [0.925, -0.038, -0.009],
    [0.950, -0.025, -0.002],
    [0.975, -0.012,  0.004],
    [1.000,  0.000,  0.000],
])

def load_cp(path):
    """Load x/c, Cp from solver output."""
    data = np.loadtxt(path, comments='#')
    return data[:, 0], data[:, 1]

def split_surfaces(xc, cp):
    """
    Split cp.dat into upper and lower surfaces.
    cp.dat is written from i_TEl (lower TE, x≈1) going:
      i_TEl → i_LE  : lower surface, x decreasing 1→0
      i_LE  → i_TEu : upper surface, x increasing 0→1
    LE is the minimum-x point.
    """
    LE_idx = np.argmin(xc)
    # Lower surface (pressure side): first half, x goes 1→0 → reverse to 0→1
    x_lo  = xc[:LE_idx+1][::-1]
    cp_lo = cp[:LE_idx+1][::-1]
    # Upper surface (suction side): second half, x already goes 0→1
    x_up  = xc[LE_idx:]
    cp_up = cp[LE_idx:]
    return x_up, cp_up, x_lo, cp_lo

def main():
    cp_file = sys.argv[1] if len(sys.argv) > 1 else "output/cp.dat"
    if not os.path.exists(cp_file):
        print(f"File not found: {cp_file}")
        sys.exit(1)

    xc, cp = load_cp(cp_file)
    x_up, cp_up, x_lo, cp_lo = split_surfaces(xc, cp)

    fig, ax = plt.subplots(figsize=(10, 7))

    # Experimental data
    ax.plot(AGARD_CASE9[:, 0], AGARD_CASE9[:, 1],
            'ko', markersize=5, label='Exp. upper (Case 9)', zorder=5)
    ax.plot(AGARD_CASE9[:, 0], AGARD_CASE9[:, 2],
            'k^', markersize=5, label='Exp. lower (Case 9)', zorder=5)

    # CFD results
    ax.plot(x_up, cp_up, 'b-', linewidth=1.5, label='CFD upper (SA-RANS)')
    ax.plot(x_lo, cp_lo, 'r-', linewidth=1.5, label='CFD lower (SA-RANS)')

    ax.set_xlabel('x/c', fontsize=13)
    ax.set_ylabel(r'$C_p$', fontsize=13)
    ax.set_title('RAE 2822  M=0.729  AoA=2.31°  Re=6.5×10⁶  (SA-RANS)\n'
                 'vs AGARD Case 9  M=0.73  AoA=2.79°', fontsize=12)
    ax.invert_yaxis()
    ax.set_xlim(-0.02, 1.02)
    ax.legend(fontsize=11)
    ax.grid(True, alpha=0.3)
    ax.axhline(0, color='k', linewidth=0.5, linestyle='--')

    out = "output/cp_comparison.png"
    plt.tight_layout()
    plt.savefig(out, dpi=150)
    print(f"Saved {out}")
    plt.show()

if __name__ == "__main__":
    main()
