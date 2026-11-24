#!/usr/bin/env python3
"""Grafica la evolución temporal a partir de los CSV de PROFILE_CSV.

Genera:
  - docs/figs/spacetime_chamber.png
  - docs/figs/evolution_chamber.png
  - docs/figs/evolution_chamber.gif
  - docs/figs/spacetime_sphere.png
"""
import os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation, PillowWriter

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DATA = os.path.join(ROOT, "data", "figs")
OUT = os.path.join(ROOT, "docs", "figs")
os.makedirs(OUT, exist_ok=True)


def read_profile(path):
    meta = {}
    with open(path) as f:
        line = f.readline().strip()
        if not line.startswith("#"):
            raise ValueError("falta la linea de metadatos")
        for tok in line[1:].split():
            k, v = tok.split("=", 1)
            meta[k] = v
        header = f.readline().strip().split(",")
        rows = [list(map(float, ln.strip().split(","))) for ln in f if ln.strip()]
    arr = np.array(rows)
    N = int(meta["N"])
    cols = {name: i for i, name in enumerate(header)}
    prof = arr[:, cols["rho_0"]:cols["rho_0"] + N]
    return meta, arr, prof


def spacetime(path_in, path_out, title):
    meta, arr, prof = read_profile(path_in)
    N = int(meta["N"])
    dx = float(meta["dx_mm"])
    t = arr[:, 1]
    x = np.arange(N) * dx

    fig, ax = plt.subplots(figsize=(7.6, 4.6), constrained_layout=True)
    im = ax.imshow(prof, origin="lower", aspect="auto",
                   extent=[0, x[-1], t[0], t[-1]], vmin=0, vmax=3,
                   cmap="viridis", interpolation="nearest")
    ax.axvline(0.25 * x[-1], color="w", ls="--", lw=1.0)
    ax.set_xlabel("x [mm]")
    ax.set_ylabel(r"$t$ [$\mu$s]")
    ax.set_title(title)
    cb = fig.colorbar(im, ax=ax)
    cb.set_label(r"$\rho$ [partículas/celda]")
    fig.savefig(path_out, dpi=180)
    plt.close(fig)


def metrics(path_in, path_out, title):
    meta, arr, prof = read_profile(path_in)
    N = int(meta["N"])
    dx = float(meta["dx_mm"])
    t = arr[:, 1]
    series = [
        (prof.max(axis=1), r"$\max_x \bar\rho(x)$ [partículas/celda]"),
        (arr[:, 7], r"centro de masa $x_{cm}$ [mm]"),
        (arr[:, 8], r"ancho RMS $x_{rms}$ [mm]"),
        (arr[:, 10], r"fracción de masa con $x \geq 0.25N$"),
    ]
    fig, axes = plt.subplots(2, 2, figsize=(10, 6), constrained_layout=True)
    for ax, (y, label) in zip(axes.ravel(), series):
        ax.plot(t, y, lw=1.6)
        ax.set_xlabel(r"$t$ [$\mu$s]")
        ax.set_ylabel(label)
        ax.grid(alpha=0.3)
    fig.suptitle(title)
    fig.savefig(path_out, dpi=180)
    plt.close(fig)


def animation(path_in, path_out):
    meta, arr, prof = read_profile(path_in)
    N = int(meta["N"])
    dx = float(meta["dx_mm"])
    t = arr[:, 1]
    x = np.arange(N) * dx

    fig, ax = plt.subplots(figsize=(6.4, 4.0), constrained_layout=True)
    ax.set_xlim(0, x[-1])
    ax.set_ylim(0, 3.2)
    ax.set_xlabel("x [mm]")
    ax.set_ylabel(r"$\rho$ [partículas/celda]")
    ax.grid(alpha=0.3)
    ax.axvline(0.25 * x[-1], color="k", ls="--", lw=1.0)
    line, = ax.plot([], [], lw=1.8)
    text = ax.text(0.02, 0.92, "", transform=ax.transAxes)

    step = max(1, len(t) // 60)
    frames = list(range(0, len(t), step))

    def update(i):
        line.set_data(x, prof[i])
        text.set_text(f"t = {t[i]:.0f} µs")
        return line, text

    ani = FuncAnimation(fig, update, frames=frames, blit=False)
    ani.save(path_out, writer=PillowWriter(fps=15))
    plt.close(fig)


def main():
    spacetime(os.path.join(DATA, "chamber_profile.csv"),
              os.path.join(OUT, "spacetime_chamber.png"),
              "Cámara 10:1, RAMP=8 — evolución espacio-temporal")
    metrics(os.path.join(DATA, "chamber_profile.csv"),
            os.path.join(OUT, "evolution_chamber.png"),
            "Cámara 10:1 — métricas globales")
    animation(os.path.join(DATA, "chamber_profile.csv"),
              os.path.join(OUT, "evolution_chamber.gif"))
    spacetime(os.path.join(DATA, "sphere_profile.csv"),
              os.path.join(OUT, "spacetime_sphere.png"),
              "Esfera densa — evolución espacio-temporal")
    print("Figuras en", OUT)


if __name__ == "__main__":
    main()
