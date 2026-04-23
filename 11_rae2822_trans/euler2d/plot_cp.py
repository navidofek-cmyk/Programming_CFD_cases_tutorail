#!/usr/bin/env python3
"""Plot Cp distribution along the airfoil from output/cp.dat"""

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import sys, os

def read_cp(path):
    data = []
    with open(path) as f:
        for line in f:
            if line.startswith('#'): continue
            parts = line.split()
            if len(parts) >= 2:
                data.append((float(parts[0]), float(parts[1])))
    return np.array(data)

if __name__ == '__main__':
    cp_file = sys.argv[1] if len(sys.argv) > 1 else 'output/cp.dat'
    data = read_cp(cp_file)
    x, cp = data[:,0], data[:,1]

    # Split into lower (x decreasing from TE to LE) and upper (x increasing from LE to TE)
    # cp.dat loops from i_TEl to i_TEu-1 (lower TE→LE, upper LE→TE)
    n = len(x)
    mid = n // 2
    x_lo, cp_lo = x[:mid], cp[:mid]
    x_up, cp_up = x[mid:], cp[mid:]

    fig, ax = plt.subplots(figsize=(10, 6), facecolor='#0d1117')
    ax.set_facecolor('#0d1117')

    # RAE 2822 convention: Cp axis inverted (negative up)
    ax.plot(x_up, cp_up, color='#ef5350', lw=2.0, label='Upper surface')
    ax.plot(x_lo, cp_lo, color='#42a5f5', lw=2.0, label='Lower surface')

    # Mark shock region (approximate)
    shock_x = x_up[np.argmax(np.diff(cp_up)) + 1] if len(x_up) > 2 else 0.6
    ax.axvline(shock_x, color='yellow', lw=0.8, ls='--', alpha=0.6, label=f'Shock x/c≈{shock_x:.2f}')

    ax.invert_yaxis()
    ax.set_xlim(-0.05, 1.05)
    ax.set_xlabel('x/c', color='#ccc', fontsize=12)
    ax.set_ylabel('Cp  (–Cp upward)', color='#ccc', fontsize=12)
    ax.tick_params(colors='#aaa')
    for sp in ax.spines.values(): sp.set_color('#333')
    ax.set_title('RAE 2822  M=0.729  AoA=2.31°  (Euler, MUSCL)', color='white', fontsize=13)
    ax.legend(facecolor='#1a1a2e', edgecolor='#444', labelcolor='white', fontsize=10)
    ax.grid(color='#222', lw=0.5)

    # stats box
    cp_min = float(cp_up.min())
    cp_shock_jump = float(max(np.diff(cp_up)))
    txt = (f"Cp_min  (upper) = {cp_min:.3f}\n"
           f"ΔCp at shock    ≈ {cp_shock_jump:.3f}")
    ax.text(0.02, 0.05, txt, transform=ax.transAxes, va='bottom',
            fontsize=9, family='monospace', color='#e0e0e0',
            bbox=dict(boxstyle='round,pad=0.4', fc='#0d1117', ec='#4fc3f7', alpha=0.9))

    os.makedirs('output', exist_ok=True)
    out = 'output/cp.png'
    fig.tight_layout()
    fig.savefig(out, dpi=150, bbox_inches='tight', facecolor=fig.get_facecolor())
    plt.close(fig)
    print(f'Saved {out}')
