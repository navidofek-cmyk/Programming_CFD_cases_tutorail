#!/usr/bin/env python3
"""Generate mesh visualisation images from output/mesh.dat"""

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import sys, os

# ---- read Tecplot mesh -------------------------------------------------------
def read_tecplot(path):
    nodes = []
    ni = nj = None
    with open(path) as f:
        for line in f:
            line = line.strip()
            if line.startswith('ZONE'):
                for tok in line.split():
                    if tok.startswith('I='):
                        ni = int(tok.split('=')[1].rstrip(','))
                    elif tok.startswith('J='):
                        nj = int(tok.split('=')[1].rstrip(','))
            elif line and not line.startswith(('TITLE','VARIABLES')):
                parts = line.split()
                if len(parts) == 2:
                    nodes.append((float(parts[0]), float(parts[1])))
    x = np.array([n[0] for n in nodes]).reshape(nj, ni)
    y = np.array([n[1] for n in nodes]).reshape(nj, ni)
    return x, y, ni, nj

# ---- mesh statistics ---------------------------------------------------------
def mesh_stats(x, y, ni, nj):
    # Cell areas (shoelace for quad)
    x0, x1 = x[:-1,:-1], x[:-1,1:]
    x2, x3 = x[1:,1:],   x[1:,:-1]
    y0, y1 = y[:-1,:-1], y[:-1,1:]
    y2, y3 = y[1:,1:],   y[1:,:-1]
    ac_x = x2-x0; ac_y = y2-y0
    bd_x = x3-x1; bd_y = y3-y1
    area = 0.5*np.abs(ac_x*bd_y - ac_y*bd_x)

    # First-cell normal height (j=0→1)
    dxj = x[1,:] - x[0,:]
    dyj = y[1,:] - y[0,:]
    h1 = np.sqrt(dxj**2 + dyj**2)

    # Streamwise spacing at wall (i-direction)
    dxi = np.diff(x[0,:])
    dyi = np.diff(y[0,:])
    ds_i = np.sqrt(dxi**2 + dyi**2)

    # Aspect ratio at wall (avoid /0)
    h1_cell = 0.5*(h1[:-1]+h1[1:])
    with np.errstate(divide='ignore', invalid='ignore'):
        ar = np.where(h1_cell > 1e-15, ds_i / h1_cell, 0.0)

    # Airfoil indices: guess from x-range
    on_airfoil = (x[0,:] >= 0.0) & (x[0,:] <= 1.0) & (np.abs(y[0,:]) < 0.15)
    idx_airfoil = np.where(on_airfoil)[0]

    stats = {
        'ni': ni, 'nj': nj,
        'n_cells': (ni-1)*(nj-1),
        'area_min': float(area[area > 0].min()),
        'area_max': float(area.max()),
        'h1_airfoil_mean': float(h1[idx_airfoil].mean()) if len(idx_airfoil) else 0.0,
        'h1_airfoil_max':  float(h1[idx_airfoil].max())  if len(idx_airfoil) else 0.0,
        'ar_max':          float(ar[np.isfinite(ar)].max()),
        'ar_mean_airfoil': float(ar[idx_airfoil[:-1]].mean()) if len(idx_airfoil) > 1 else 0.0,
    }
    return stats, area

# ---- draw a subset of grid lines ---------------------------------------------
def draw_grid(ax, x, y, ni, nj, skip_i=4, skip_j=2, **kw):
    for j in range(0, nj, skip_j):
        ax.plot(x[j,:], y[j,:], **kw)
    for i in range(0, ni, skip_i):
        ax.plot(x[:,i], y[:,i], **kw)

def stats_box(ax, txt, x=0.02, y=0.97, ha='left', color='#4fc3f7'):
    ax.text(x, y, txt, transform=ax.transAxes, va='top', ha=ha,
            fontsize=8, family='monospace', color='#e0e0e0',
            bbox=dict(boxstyle='round,pad=0.4', fc='#0d1117', ec=color, alpha=0.90))

def style_ax(ax, title):
    ax.set_aspect('equal')
    ax.set_xlabel('x/c', color='#ccc'); ax.set_ylabel('y/c', color='#ccc')
    ax.tick_params(colors='#aaa', labelsize=8)
    for sp in ax.spines.values(): sp.set_color('#333')
    ax.set_title(title, color='white', fontsize=12, pad=6)
    ax.set_facecolor('#0d1117')

# ---- plot 1: full mesh -------------------------------------------------------
def plot_full(x, y, ni, nj, stats, out):
    fig, ax = plt.subplots(figsize=(9, 9), facecolor='#0d1117')

    # Every 2nd j-line, every 4th i-line so individual cells are visible
    draw_grid(ax, x, y, ni, nj, skip_i=4, skip_j=2,
              color='#29b6f6', lw=0.4, alpha=0.7)
    ax.plot(x[0,:],  y[0,:],  color='#ff7043', lw=1.4, label='wall / wake cut', zorder=5)
    ax.plot(x[-1,:], y[-1,:], color='#66bb6a', lw=1.4, label='farfield', zorder=5)

    ax.set_xlim(-28, 28); ax.set_ylim(-28, 28)
    style_ax(ax, 'RAE 2822 C-grid — full domain')
    ax.legend(facecolor='#1a1a2e', edgecolor='#444', labelcolor='white',
              fontsize=9, loc='lower right')

    txt = (f"Nodes         {stats['ni']} × {stats['nj']}\n"
           f"Cells         {stats['n_cells']:,}\n"
           f"Area  min     {stats['area_min']:.2e} c²\n"
           f"Area  max     {stats['area_max']:.3f} c²\n"
           f"h₁ wall avg   {stats['h1_airfoil_mean']:.4f} c\n"
           f"AR  max       {stats['ar_max']:.0f}")
    stats_box(ax, txt)

    fig.tight_layout()
    fig.savefig(out, dpi=160, bbox_inches='tight', facecolor=fig.get_facecolor())
    plt.close(fig)
    print(f'Saved {out}')

# ---- plot 2: airfoil closeup ------------------------------------------------
def plot_closeup(x, y, ni, nj, stats, out):
    fig, axes = plt.subplots(1, 2, figsize=(15, 6), facecolor='#0d1117')

    # Left: medium closeup
    ax = axes[0]
    draw_grid(ax, x, y, ni, nj, skip_i=1, skip_j=1,
              color='#42a5f5', lw=0.45, alpha=0.65)
    ax.plot(x[0,:], y[0,:], color='#ff7043', lw=1.6, zorder=5)
    ax.set_xlim(-0.15, 1.2); ax.set_ylim(-0.35, 0.35)
    style_ax(ax, 'Airfoil region — all grid lines')
    txt = (f"h₁ avg  {stats['h1_airfoil_mean']:.5f} c\n"
           f"h₁ max  {stats['h1_airfoil_max']:.5f} c\n"
           f"AR avg  {stats['ar_mean_airfoil']:.1f}")
    stats_box(ax, txt, x=0.98, ha='right', color='#42a5f5')

    # Right: wall zoom (first ~10 % chord, first ~15 j-lines)
    ax2 = axes[1]
    jmax = 15
    draw_grid(ax2, x[:jmax+1,:], y[:jmax+1,:], ni, jmax+1, skip_i=1, skip_j=1,
              color='#66bb6a', lw=0.6, alpha=0.8)
    ax2.plot(x[0,:], y[0,:], color='#ff7043', lw=1.8, zorder=5)
    ax2.set_xlim(0.35, 0.75); ax2.set_ylim(-0.04, 0.12)
    style_ax(ax2, 'Near-wall detail (first 15 cells, x/c 0.35–0.75)')
    # annotate first cell height
    i_mid = ni // 2
    ax2.annotate('', xy=(x[1,i_mid], y[1,i_mid]), xytext=(x[0,i_mid], y[0,i_mid]),
                 arrowprops=dict(arrowstyle='<->', color='yellow', lw=1.5))
    h1_mid = np.sqrt((x[1,i_mid]-x[0,i_mid])**2 + (y[1,i_mid]-y[0,i_mid])**2)
    ax2.text(x[0,i_mid]+0.01, (y[0,i_mid]+y[1,i_mid])/2,
             f'h₁={h1_mid:.4f}c', color='yellow', fontsize=8)

    fig.suptitle('RAE 2822 mesh quality', color='white', fontsize=13, y=1.01)
    fig.tight_layout()
    fig.savefig(out, dpi=160, bbox_inches='tight', facecolor=fig.get_facecolor())
    plt.close(fig)
    print(f'Saved {out}')

# ---- main -------------------------------------------------------------------
if __name__ == '__main__':
    dat = sys.argv[1] if len(sys.argv) > 1 else 'output/mesh.dat'
    x, y, ni, nj = read_tecplot(dat)
    stats, area = mesh_stats(x, y, ni, nj)

    print("=== Mesh statistics ===")
    for k, v in stats.items():
        print(f"  {k:20s}: {v}")

    os.makedirs('output', exist_ok=True)
    plot_full    (x, y, ni, nj, stats, 'output/mesh_full.png')
    plot_closeup (x, y, ni, nj, stats, 'output/mesh_closeup.png')
