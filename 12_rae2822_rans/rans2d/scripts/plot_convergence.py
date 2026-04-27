#!/usr/bin/env python3
"""
Plot solver convergence history from output/convergence.dat.

Usage:
    python3 scripts/plot_convergence.py [output/convergence.dat]
"""
import sys
import numpy as np
import matplotlib.pyplot as plt
import os

def load_convergence(path):
    iters, res, res_norm = [], [], []
    warmup_end = None
    with open(path) as f:
        for line in f:
            if line.startswith('#'):
                if 'order switch' in line:
                    # mark where warmup ended
                    if iters:
                        warmup_end = iters[-1]
                continue
            parts = line.split()
            if len(parts) < 3:
                continue
            iters.append(int(parts[0]))
            res.append(float(parts[1]))
            res_norm.append(float(parts[2]))
    return (np.array(iters), np.array(res), np.array(res_norm), warmup_end)

def main():
    conv_file = sys.argv[1] if len(sys.argv) > 1 else "output/convergence.dat"
    if not os.path.exists(conv_file):
        print(f"File not found: {conv_file}")
        sys.exit(1)

    iters, res, res_norm, warmup_end = load_convergence(conv_file)

    fig, axes = plt.subplots(1, 2, figsize=(14, 5))

    for ax, y, ylabel in zip(
        axes,
        [res, res_norm],
        [r'$\|\mathbf{R}_\rho\|_2$  (absolute)', r'$\|\mathbf{R}_\rho\|/\|\mathbf{R}_0\|$  (relative)']
    ):
        ax.semilogy(iters, y, 'b-', linewidth=0.8)
        if warmup_end is not None:
            ax.axvline(warmup_end, color='orange', linewidth=1.5,
                       linestyle='--', label=f'Order switch (iter {warmup_end})')
            ax.legend(fontsize=10)
        ax.set_xlabel('Iteration', fontsize=12)
        ax.set_ylabel(ylabel, fontsize=12)
        ax.grid(True, which='both', alpha=0.3)

    axes[0].set_title('Absolute residual', fontsize=12)
    axes[1].set_title('Relative residual', fontsize=12)

    fig.suptitle('RAE 2822  SA-RANS convergence history  M=0.729  AoA=2.31°  Re=6.5×10⁶',
                 fontsize=13)
    plt.tight_layout()

    out = "output/convergence.png"
    plt.savefig(out, dpi=150)
    print(f"Saved {out}")
    plt.show()

if __name__ == "__main__":
    main()
