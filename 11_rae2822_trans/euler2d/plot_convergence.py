#!/usr/bin/env python3
"""Plot convergence history from output/convergence.dat"""

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import sys, os

if __name__ == '__main__':
    conv_file = sys.argv[1] if len(sys.argv) > 1 else 'output/convergence.dat'

    iters, res, res_norm = [], [], []
    switch_iter = None
    with open(conv_file) as f:
        for line in f:
            if line.startswith('#'):
                if 'switched' in line:
                    # extract iter number
                    parts = line.split('iter')
                    if len(parts) > 1:
                        try: switch_iter = int(parts[-1].strip())
                        except: pass
                continue
            parts = line.split()
            if len(parts) >= 3:
                iters.append(int(parts[0]))
                res.append(float(parts[1]))
                res_norm.append(float(parts[2]))

    iters = np.array(iters)
    res   = np.array(res)
    res_n = np.array(res_norm)

    fig, axes = plt.subplots(1, 2, figsize=(14, 5), facecolor='#0d1117')

    for ax, y, ylabel, title in zip(
            axes,
            [res, res_n],
            ['L2 density residual (abs)', 'Residual / Res₀'],
            ['Absolute residual', 'Normalized residual']):
        ax.set_facecolor('#0d1117')
        ax.semilogy(iters, y, color='#29b6f6', lw=1.2)
        if switch_iter:
            ax.axvline(switch_iter, color='#ff7043', lw=1.2, ls='--',
                       label=f'Order switch @ iter {switch_iter}')
            ax.legend(facecolor='#1a1a2e', edgecolor='#444', labelcolor='white', fontsize=9)
        ax.set_xlabel('Iteration', color='#ccc')
        ax.set_ylabel(ylabel, color='#ccc')
        ax.set_title(title, color='white')
        ax.tick_params(colors='#aaa')
        for sp in ax.spines.values(): sp.set_color('#333')
        ax.grid(color='#222', which='both', lw=0.5)

    fig.suptitle('RAE 2822 Euler MUSCL — Convergence', color='white', fontsize=13)
    fig.tight_layout()
    out = 'output/convergence.png'
    fig.savefig(out, dpi=150, bbox_inches='tight', facecolor=fig.get_facecolor())
    plt.close(fig)
    print(f'Saved {out}')
