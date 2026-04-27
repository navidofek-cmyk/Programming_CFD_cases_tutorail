#!/usr/bin/env python3
"""
Read output/field.vts (VTK StructuredGrid XML, ASCII) and plot 2D contour maps
using pcolormesh on the structured grid.

Usage:
    python3 scripts/plot_field.py [field=mach] [output/field.vts]

Available fields: density, pressure, mach, cp, temperature, entropy, mu_t_ratio
"""
import sys
import os
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import xml.etree.ElementTree as ET

FIELD_OPTS = {
    'mach':       dict(cmap='RdBu_r',   vmin=0.0,  vmax=1.4,  label='Mach number'),
    'cp':         dict(cmap='RdBu_r',   vmin=-1.5, vmax=1.0,  label=r'$C_p$'),
    'density':    dict(cmap='viridis',  vmin=None, vmax=None,  label=r'$\rho/\rho_\infty$'),
    'pressure':   dict(cmap='plasma',   vmin=None, vmax=None,  label=r'$p / p_\infty$'),
    'mu_t_ratio': dict(cmap='hot_r',    vmin=0,    vmax=200,   label=r'$\mu_t / \mu_{lam}$'),
    'entropy':    dict(cmap='coolwarm', vmin=None, vmax=None,  label=r'$\ln(p/\rho^\gamma)$'),
    'temperature':dict(cmap='inferno',  vmin=None, vmax=None,  label=r'$T / T_\infty$'),
}

def read_vts(path):
    tree = ET.parse(path)
    root = tree.getroot()
    sg = root.find('StructuredGrid')
    extent = [int(v) for v in sg.attrib['WholeExtent'].split()]
    ni = extent[1] - extent[0] + 1   # node count i
    nj = extent[3] - extent[2] + 1   # node count j

    piece = sg.find('Piece')

    pts_da = piece.find('Points/DataArray')
    pts = np.fromstring(pts_da.text, sep=' ').reshape(-1, 3)
    # nodes shaped (nj, ni) — outer j = farfield, j=0 = wall
    xn = pts[:, 0].reshape(nj, ni)
    yn = pts[:, 1].reshape(nj, ni)

    cell_data = {}
    cd = piece.find('CellData')
    if cd is not None:
        nci, ncj = ni - 1, nj - 1
        for da in cd.findall('DataArray'):
            name = da.attrib['Name']
            nc = int(da.attrib.get('NumberOfComponents', '1'))
            raw = np.fromstring(da.text, sep=' ' if ' ' in da.text.strip() else '\n')
            if nc == 3:
                raw = np.linalg.norm(raw.reshape(-1, 3), axis=1)
            cell_data[name] = raw.reshape(ncj, nci)

    return xn, yn, cell_data

def make_figure(xn, yn, vals2d, opts, title, xlim=(-0.2, 1.2), ylim=(-0.5, 0.5)):
    fig, ax = plt.subplots(figsize=(13, 6))

    # Clip to view region to avoid stretched outer cells
    xc = 0.25*(xn[:-1,:-1]+xn[1:,:-1]+xn[1:,1:]+xn[:-1,1:])
    yc = 0.25*(yn[:-1,:-1]+yn[1:,:-1]+yn[1:,1:]+yn[:-1,1:])
    mask = (xc >= xlim[0]) & (xc <= xlim[1]) & (yc >= ylim[0]) & (yc <= ylim[1])

    # Find j-range covering the view
    j_ok = np.where(mask.any(axis=1))[0]
    if len(j_ok) == 0:
        j0, j1 = 0, vals2d.shape[0]
    else:
        j0, j1 = j_ok[0], j_ok[-1]+1

    pcm = ax.pcolormesh(xn[j0:j1+1, :], yn[j0:j1+1, :], vals2d[j0:j1, :],
                        cmap=opts['cmap'], shading='flat',
                        vmin=opts.get('vmin'), vmax=opts.get('vmax'),
                        rasterized=True)
    plt.colorbar(pcm, ax=ax, label=opts['label'], fraction=0.03, pad=0.02)

    ax.set_xlim(xlim)
    ax.set_ylim(ylim)
    ax.set_aspect('equal')
    ax.set_xlabel('x/c', fontsize=12)
    ax.set_ylabel('y/c', fontsize=12)
    ax.set_title(title, fontsize=12)
    return fig, ax

def main():
    field_name = 'mach'
    vts_path = 'output/field.vts'
    for arg in sys.argv[1:]:
        if arg.endswith('.vts'):
            vts_path = arg
        else:
            field_name = arg

    if not os.path.exists(vts_path):
        print(f"File not found: {vts_path}")
        sys.exit(1)

    print(f"Reading {vts_path} ...")
    xn, yn, cell_data = read_vts(vts_path)
    print(f"Available fields: {list(cell_data.keys())}")

    if field_name not in cell_data:
        print(f"Field '{field_name}' not found. Available: {list(cell_data.keys())}")
        sys.exit(1)

    opts = FIELD_OPTS.get(field_name, dict(cmap='viridis', vmin=None, vmax=None, label=field_name))
    title = (f'RAE 2822  SA-RANS — {opts["label"]}'
             f'  M=0.729  AoA=2.31°  Re=6.5×10⁶')

    vals2d = cell_data[field_name]

    # Full domain view
    fig1, ax1 = make_figure(xn, yn, vals2d, opts, title, xlim=(-0.2, 1.3), ylim=(-0.5, 0.5))
    out1 = f"output/{field_name}.png"
    fig1.tight_layout()
    fig1.savefig(out1, dpi=150)
    print(f"Saved {out1}")

    # Near-airfoil zoom
    fig2, ax2 = make_figure(xn, yn, vals2d, opts, title + ' [zoom]',
                             xlim=(-0.05, 1.05), ylim=(-0.12, 0.15))
    out2 = f"output/{field_name}_zoom.png"
    fig2.tight_layout()
    fig2.savefig(out2, dpi=150)
    print(f"Saved {out2}")

    plt.show()

if __name__ == "__main__":
    main()
