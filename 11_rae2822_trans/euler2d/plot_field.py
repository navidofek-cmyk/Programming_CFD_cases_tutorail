#!/usr/bin/env python3
"""Read field.vts and produce contour plots: Mach, velocity magnitude, density."""

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import sys, os, re

# ---- parse VTK StructuredGrid (.vts) ASCII -----------------------------------

def read_vts(path):
    with open(path) as f:
        text = f.read()

    # Grid dimensions from WholeExtent="0 NI 0 NJ 0 0"
    m = re.search(r'WholeExtent="0\s+(\d+)\s+0\s+(\d+)\s+0\s+0"', text)
    nci, ncj = int(m.group(1)), int(m.group(2))  # cell counts
    ni, nj   = nci + 1, ncj + 1                  # node counts

    def extract_array(name, ncomp=1):
        # find DataArray with Name="name"
        pat = rf'Name="{name}"[^>]*>(.*?)</DataArray>'
        m = re.search(pat, text, re.DOTALL)
        if not m:
            return None
        vals = np.fromstring(m.group(1), sep='\n' if '\n' in m.group(1) else ' ')
        if ncomp > 1:
            vals = vals.reshape(-1, ncomp)
        return vals

    # Node coordinates
    pts_raw = extract_array('', ncomp=3)   # Points section
    # Alternative: parse Points block directly
    pts_m = re.search(r'<Points>.*?<DataArray[^>]*>(.*?)</DataArray>',
                      text, re.DOTALL)
    pts = np.fromstring(pts_m.group(1), sep='\n').reshape(ni * nj, 3)
    X = pts[:, 0].reshape(nj, ni)
    Y = pts[:, 1].reshape(nj, ni)

    fields = {}
    for name in ('density', 'pressure', 'mach', 'cp', 'temperature', 'entropy'):
        arr = extract_array(name)
        if arr is not None:
            fields[name] = arr.reshape(ncj, nci)

    vel_m = re.search(r'Name="velocity"[^>]*>(.*?)</DataArray>', text, re.DOTALL)
    if vel_m:
        vel = np.fromstring(vel_m.group(1), sep='\n').reshape(ncj * nci, 3)
        fields['vel_u'] = vel[:, 0].reshape(ncj, nci)
        fields['vel_v'] = vel[:, 1].reshape(ncj, nci)
        fields['vel_mag'] = np.sqrt(vel[:, 0]**2 + vel[:, 1]**2).reshape(ncj, nci)

    return X, Y, ni, nj, nci, ncj, fields


# ---- plotting helpers --------------------------------------------------------

def style_ax(ax, title, cbar_label, cmap, data, X, Y,
             levels=None, vmin=None, vmax=None, clip_j=None):
    """Filled contour plot on the C-grid."""
    # Optionally clip j range (focus on near-airfoil region)
    jmax = clip_j if clip_j else -1
    Xp = X[:jmax+1, :] if clip_j else X
    Yp = Y[:jmax+1, :] if clip_j else Y
    Dp = data[:jmax, :] if clip_j else data

    # Cell-center coordinates (average of 4 nodes)
    xc = 0.25 * (Xp[:-1, :-1] + Xp[:-1, 1:] + Xp[1:, :-1] + Xp[1:, 1:])
    yc = 0.25 * (Yp[:-1, :-1] + Yp[:-1, 1:] + Yp[1:, :-1] + Yp[1:, 1:])

    if levels is not None:
        cf = ax.contourf(xc, yc, Dp, levels=levels, cmap=cmap, extend='both')
    else:
        cf = ax.pcolormesh(Xp, Yp, Dp, cmap=cmap, vmin=vmin, vmax=vmax,
                           shading='auto', rasterized=True)

    cb = plt.colorbar(cf, ax=ax, fraction=0.03, pad=0.02)
    cb.set_label(cbar_label, color='#ccc', fontsize=9)
    cb.ax.yaxis.set_tick_params(color='#aaa', labelsize=8)
    plt.setp(cb.ax.yaxis.get_ticklabels(), color='#aaa')

    ax.set_aspect('equal')
    ax.set_facecolor('#0d1117')
    ax.set_title(title, color='white', fontsize=11, pad=5)
    ax.set_xlabel('x/c', color='#ccc', fontsize=9)
    ax.set_ylabel('y/c', color='#ccc', fontsize=9)
    ax.tick_params(colors='#aaa', labelsize=8)
    for sp in ax.spines.values(): sp.set_color('#333')

    # Airfoil wall (j=0 line)
    ax.plot(Xp[0, :], Yp[0, :], color='#ff7043', lw=1.2, zorder=5)
    return cf


# ---- main -------------------------------------------------------------------

if __name__ == '__main__':
    vts = sys.argv[1] if len(sys.argv) > 1 else 'output/field.vts'
    X, Y, ni, nj, nci, ncj, F = read_vts(vts)
    os.makedirs('output', exist_ok=True)

    print(f"Grid: {nci}x{ncj} cells, {ni}x{nj} nodes")
    print(f"Mach range:    [{F['mach'].min():.3f}, {F['mach'].max():.3f}]")
    print(f"Density range: [{F['density'].min():.3f}, {F['density'].max():.3f}]")
    print(f"Cp range:      [{F['cp'].min():.3f}, {F['cp'].max():.3f}]")

    # ---- Figure 1: Mach + sonic line ----------------------------------------
    fig, axes = plt.subplots(1, 2, figsize=(16, 6), facecolor='#0d1117')

    # Full domain
    ax = axes[0]
    lvl_m = np.linspace(0.0, 1.6, 33)
    cf = style_ax(ax, 'Mach number — full domain', 'M', 'RdBu_r',
                  F['mach'], X, Y, levels=lvl_m)
    # Sonic line
    xc_all = 0.25*(X[:-1,:-1]+X[:-1,1:]+X[1:,:-1]+X[1:,1:])
    yc_all = 0.25*(Y[:-1,:-1]+Y[:-1,1:]+Y[1:,:-1]+Y[1:,1:])
    ax.contour(xc_all, yc_all, F['mach'], levels=[1.0],
               colors='yellow', linewidths=1.2, linestyles='--')
    ax.set_xlim(-2, 3); ax.set_ylim(-2, 2)

    # Closeup airfoil
    ax2 = axes[1]
    style_ax(ax2, 'Mach number — airfoil region', 'M', 'RdBu_r',
             F['mach'], X, Y, levels=lvl_m, clip_j=25)
    xc25 = 0.25*(X[:25,:-1]+X[:25,1:]+X[1:26,:-1]+X[1:26,1:])
    yc25 = 0.25*(Y[:25,:-1]+Y[:25,1:]+Y[1:26,:-1]+Y[1:26,1:])
    ax2.contour(xc25, yc25, F['mach'][:25,:], levels=[1.0],
                colors='yellow', linewidths=1.5, linestyles='--')
    ax2.set_xlim(-0.1, 1.15); ax2.set_ylim(-0.2, 0.35)

    fig.suptitle('RAE 2822  M=0.729  AoA=2.31°  (Euler MUSCL)', color='white', fontsize=13)
    fig.tight_layout()
    fig.savefig('output/field_mach.png', dpi=160, bbox_inches='tight',
                facecolor=fig.get_facecolor())
    plt.close(fig)
    print('Saved output/field_mach.png')

    # ---- Figure 2: Density --------------------------------------------------
    fig, axes = plt.subplots(1, 2, figsize=(16, 6), facecolor='#0d1117')

    lvl_d = np.linspace(0.7, 1.4, 29)
    style_ax(axes[0], 'Density — full domain', 'ρ / ρ∞', 'plasma',
             F['density'], X, Y, levels=lvl_d)
    axes[0].set_xlim(-2, 3); axes[0].set_ylim(-2, 2)

    style_ax(axes[1], 'Density — airfoil region', 'ρ / ρ∞', 'plasma',
             F['density'], X, Y, levels=lvl_d, clip_j=25)
    axes[1].set_xlim(-0.1, 1.15); axes[1].set_ylim(-0.2, 0.35)

    fig.suptitle('RAE 2822  M=0.729  AoA=2.31°  (Euler MUSCL)', color='white', fontsize=13)
    fig.tight_layout()
    fig.savefig('output/field_density.png', dpi=160, bbox_inches='tight',
                facecolor=fig.get_facecolor())
    plt.close(fig)
    print('Saved output/field_density.png')

    # ---- Figure 3: Velocity magnitude + vectors (coarse) --------------------
    fig, axes = plt.subplots(1, 2, figsize=(16, 6), facecolor='#0d1117')

    lvl_v = np.linspace(0.0, 1.1, 23)
    style_ax(axes[0], 'Velocity magnitude — full domain', '|u| / a∞', 'viridis',
             F['vel_mag'], X, Y, levels=lvl_v)
    axes[0].set_xlim(-2, 3); axes[0].set_ylim(-2, 2)

    style_ax(axes[1], 'Velocity magnitude — airfoil region', '|u| / a∞', 'viridis',
             F['vel_mag'], X, Y, levels=lvl_v, clip_j=25)
    axes[1].set_xlim(-0.1, 1.15); axes[1].set_ylim(-0.2, 0.35)

    # Add streamwise velocity vectors (coarsened)
    skip = 8
    xc25 = 0.25*(X[:25,:-1]+X[:25,1:]+X[1:26,:-1]+X[1:26,1:])
    yc25 = 0.25*(Y[:25,:-1]+Y[:25,1:]+Y[1:26,:-1]+Y[1:26,1:])
    axes[1].quiver(xc25[::skip, ::skip], yc25[::skip, ::skip],
                   F['vel_u'][:25:skip, ::skip], F['vel_v'][:25:skip, ::skip],
                   color='white', alpha=0.5, scale=20, width=0.002)

    fig.suptitle('RAE 2822  M=0.729  AoA=2.31°  (Euler MUSCL)', color='white', fontsize=13)
    fig.tight_layout()
    fig.savefig('output/field_velocity.png', dpi=160, bbox_inches='tight',
                facecolor=fig.get_facecolor())
    plt.close(fig)
    print('Saved output/field_velocity.png')

    # ---- Figure 4: Cp + pressure --------------------------------------------
    fig, axes = plt.subplots(1, 2, figsize=(16, 6), facecolor='#0d1117')

    lvl_cp = np.linspace(-1.6, 0.8, 25)
    style_ax(axes[0], 'Pressure coefficient Cp', 'Cp', 'seismic',
             F['cp'], X, Y, levels=lvl_cp, clip_j=25)
    axes[0].set_xlim(-0.1, 1.15); axes[0].set_ylim(-0.2, 0.35)

    lvl_p = np.linspace(0.45, 0.85, 25)
    style_ax(axes[1], 'Static pressure p / p∞·γ', 'p', 'coolwarm',
             F['pressure'], X, Y, levels=lvl_p, clip_j=25)
    axes[1].set_xlim(-0.1, 1.15); axes[1].set_ylim(-0.2, 0.35)

    fig.suptitle('RAE 2822  M=0.729  AoA=2.31°  (Euler MUSCL)', color='white', fontsize=13)
    fig.tight_layout()
    fig.savefig('output/field_cp_pressure.png', dpi=160, bbox_inches='tight',
                facecolor=fig.get_facecolor())
    plt.close(fig)
    print('Saved output/field_cp_pressure.png')
