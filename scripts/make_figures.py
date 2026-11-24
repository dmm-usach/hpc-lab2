#!/usr/bin/env python3
"""Genera las figuras del informe y la presentación.

Requiere que existan los .vti en data/figs (ver README). No necesita MPI.
"""
import os
import struct

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DATA = os.path.join(ROOT, "data", "figs")
OUT = os.path.join(ROOT, "docs", "figs")
os.makedirs(OUT, exist_ok=True)

DX_MM = 1.0          # tamaño de celda usado en las corridas de figura
DT_US = 1.64217      # paso temporal del neón con dx=1 mm
RHO_P0 = 3.0         # densidad lattice que corresponde a p0
P0_ATM = 1.0


def read_vti(path):
    raw = open(path, "rb").read()
    j = raw.index(b"_", raw.index(b"<AppendedData")) + 1
    n = struct.unpack("<I", raw[j:j + 4])[0]
    data = np.frombuffer(raw[j + 4:j + 4 + n], dtype=np.uint8)
    N = int(round(n ** (1.0 / 3.0)))
    if N ** 3 != n:
        raise ValueError(f"{path}: se esperaba un cubo, n={n}")
    # VTK guarda x rápido: arr[z, y, x]
    return data.reshape((N, N, N))


def block_mean(a, k=2):
    h, w = a.shape
    h -= h % k
    w -= w % k
    return a[:h, :w].reshape(h // k, k, w // k, k).mean(axis=(1, 3))


def save_slices(prefix, steps, path, ncols=2):
    n = len(steps)
    nrows = int(np.ceil(n / ncols))
    fig, axes = plt.subplots(nrows, ncols, figsize=(4.6 * ncols, 4.2 * nrows),
                             constrained_layout=True)
    axes = np.atleast_1d(axes).ravel()
    for ax, step in zip(axes, steps):
        a = read_vti(os.path.join(DATA, f"{prefix}_{step:05d}.vti"))
        N = a.shape[0]
        # promedio de una losa de 5 capas en z y suavizado 2x2
        slab = a[max(0, N // 2 - 2):N // 2 + 3, :, :].mean(axis=0)
        sl = block_mean(slab, 2)
        im = ax.imshow(sl, origin="lower", extent=[0, N * DX_MM, 0, N * DX_MM],
                       vmin=0, vmax=3, cmap="viridis", interpolation="nearest")
        ax.set_title(f"t = {step} pasos ({step * DT_US:.0f} µs)")
        ax.set_xlabel("x [mm]")
        ax.set_ylabel("y [mm]")
        ax.set_aspect("equal")
        if "chamber" in prefix:
            ax.axvline(0.25 * N * DX_MM, color="w", ls="--", lw=0.8, alpha=0.8)
    for ax in axes[n:]:
        ax.axis("off")
    cb = fig.colorbar(im, ax=axes.tolist(), shrink=0.85, pad=0.02)
    cb.set_label(r"$\rho$ [partículas/celda]")
    fig.savefig(path, dpi=180)
    plt.close(fig)


def save_chamber_profile():
    steps = [0, 50, 100, 150]
    fig, ax = plt.subplots(figsize=(7.2, 4.2), constrained_layout=True)
    for step in steps:
        a = read_vti(os.path.join(DATA, f"chamber2_{step:05d}.vti"))
        N = a.shape[0]
        prof = a.mean(axis=(0, 1))  # promedio en y,z -> perfil en x
        ax.plot(np.arange(N) * DX_MM, prof,
                label=f"t = {step} pasos ({step * DT_US:.0f} µs)")
    N = a.shape[0]
    ax.axvline(0.25 * N * DX_MM, color="k", ls="--", lw=1.0)
    ax.text(0.25 * N * DX_MM + 2, 2.55, "interfaz inicial", fontsize=8)
    ax.set_xlabel("x [mm]")
    ax.set_ylabel(r"$\rho$ [partículas/celda]")
    ax.set_xlim(0, N * DX_MM)
    ax.set_ylim(0, 3.1)
    ax.grid(alpha=0.3)
    ax.legend(fontsize=8)
    sec = ax.secondary_yaxis(
        "right",
        functions=(lambda r: r / RHO_P0 * P0_ATM,
                   lambda p: p / P0_ATM * RHO_P0))
    sec.set_ylabel("p [atm]")
    fig.savefig(os.path.join(OUT, "chamber_profile.png"), dpi=180)
    plt.close(fig)


def save_scaling():
    # N=256, 200 pasos, laptop 8C/16T (Docker), sin VTK.
    data = [
        ("1×1", 1, 1, 1, 21.268646, 0.010274),
        ("1×2", 1, 2, 2, 11.114637, 0.010950),
        ("1×4", 1, 4, 4, 6.187174, 0.012324),
        ("1×8", 1, 8, 8, 3.717493, 0.019234),
        ("1×16", 1, 16, 16, 4.040582, 0.016197),
        ("2×4", 2, 4, 8, 6.404070, 0.790496),
        ("4×2", 4, 2, 8, 6.760879, 3.049404),
        ("4×4", 4, 4, 16, 4.057614, 0.298536),
        ("8×2", 8, 2, 16, 4.182729, 0.418648),
    ]
    base = data[0][4]
    labels = [d[0] for d in data]
    times = [d[4] for d in data]
    threads = [d[3] for d in data]
    comm = [d[5] for d in data]
    cmap = {1: "tab:blue", 2: "tab:red", 4: "tab:orange", 8: "tab:green"}
    colors = [cmap[d[1]] for d in data]

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11, 4.2),
                                   constrained_layout=True)
    ax1.bar(labels, times, color=colors)
    ax1.set_ylabel("tiempo del loop [s]")
    ax1.set_title("N=256, 200 pasos (laptop 8C/16T)")
    ax1.tick_params(axis="x", rotation=45)
    ax1.grid(axis="y", alpha=0.3)
    for i, t in enumerate(times):
        ax1.text(i, t + 0.3, f"{t:.2f}", ha="center", fontsize=8)

    speed = [base / t for t in times]
    ax2.plot([1, 16], [1, 16], "k--", lw=1, label="ideal lineal")
    ax2.scatter(threads, speed, c=colors, s=65, zorder=3)
    for d, s in zip(data, speed):
        ax2.annotate(d[0], (d[3], s), textcoords="offset points",
                     xytext=(4, 4), fontsize=8)
    ax2.set_xlabel("hilos totales (P×T)")
    ax2.set_ylabel("speedup vs 1×1")
    ax2.set_xlim(0, 17)
    ax2.set_ylim(0, 17)
    ax2.grid(alpha=0.3)
    ax2.legend(fontsize=8)
    fig.savefig(os.path.join(OUT, "scaling.png"), dpi=180)
    plt.close(fig)

    # Fracción de comunicación (solo configuraciones con P>1)
    sel = [d for d in data if d[1] > 1]
    fig, ax = plt.subplots(figsize=(6.4, 3.4), constrained_layout=True)
    frac = [100.0 * d[5] / d[4] for d in sel]
    ax.bar([d[0] for d in sel], frac, color=[cmap[d[1]] for d in sel])
    ax.set_ylabel("comunicación / tiempo [%]")
    ax.set_title("Costo de bsp_sync + halos")
    ax.grid(axis="y", alpha=0.3)
    for i, f in enumerate(frac):
        ax.text(i, f + 1.0, f"{f:.1f}%", ha="center", fontsize=8)
    fig.savefig(os.path.join(OUT, "comm_fraction.png"), dpi=180)
    plt.close(fig)


def main():
    save_slices("chamber2", [0, 50, 100, 150],
                os.path.join(OUT, "chamber_slices.png"))
    save_slices("chamber2", [0, 50, 100],
                os.path.join(OUT, "chamber_slices_wide.png"), ncols=3)
    save_slices("chamber2", [0, 100],
                os.path.join(OUT, "chamber_slices_wide2.png"), ncols=2)
    save_slices("sphere2", [0, 25, 50, 75],
                os.path.join(OUT, "sphere_slices.png"))
    save_chamber_profile()
    save_scaling()
    print("Figuras en", OUT)


if __name__ == "__main__":
    main()
