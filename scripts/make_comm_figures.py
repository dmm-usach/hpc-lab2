#!/usr/bin/env python3
"""Esquemas de descomposición de dominio y comunicación.

Genera:
  docs/figs/decomposition.png
  docs/figs/halo_exchange.png
  docs/figs/communication_volume.png
  docs/figs/hybrid_threads.png
"""
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle, FancyArrowPatch, Patch, Polygon

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT = os.path.join(ROOT, "docs", "figs")
os.makedirs(OUT, exist_ok=True)

C_PROC = "#cfe2f3"
C_PROC_EDGE = "#3d6b8e"
C_GHOST = "#f6b26b"
C_GHOST_EDGE = "#b45f06"
C_COMM = "#cc4125"
C_LOCAL = "#93c47d"
C_LOCAL_EDGE = "#38761d"
C_COMP = "#9fc5e8"
C_BARRIER = "#cc0000"
C_TEXT = "#222222"

plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 10})


def box(ax, x, y, w, h, fc, ec, lw=1.2, ls="-", alpha=1.0):
    ax.add_patch(Rectangle((x, y), w, h, facecolor=fc, edgecolor=ec,
                           linewidth=lw, linestyle=ls, alpha=alpha))


def arrow(ax, p0, p1, color=C_COMM, lw=1.8, style="-|>", ls="-"):
    ax.add_patch(FancyArrowPatch(p0, p1, arrowstyle=style, mutation_scale=15,
                                 linewidth=lw, color=color, linestyle=ls,
                                 shrinkA=0, shrinkB=0))


def clean(ax, xlim, ylim):
    ax.set_xlim(*xlim)
    ax.set_ylim(*ylim)
    ax.set_aspect("equal")
    ax.axis("off")


# ---------------------------------------------------------------------------
def decomposition():
    fig, ax = plt.subplots(figsize=(12, 4.3))
    P, w, h, gap, x0, y0 = 4, 2.45, 1.55, 0.35, 0.4, 1.35
    xr = x0 + P * (w + gap) - gap

    for r in range(P):
        x = x0 + r * (w + gap)
        box(ax, x, y0, w, h, C_PROC, C_PROC_EDGE)
        box(ax, x, y0, 0.18, h, C_GHOST, C_GHOST_EDGE, lw=0.8)
        box(ax, x + w - 0.18, y0, 0.18, h, C_GHOST, C_GHOST_EDGE, lw=0.8)
        ax.text(x + w / 2, y0 + h / 2 + 0.10, f"rank {r}",
                ha="center", va="center", fontsize=12)
        ax.text(x + w / 2, y0 + h / 2 - 0.24, "nx capas",
                ha="center", va="center", fontsize=8, color="#555")
        ax.text(x + w / 2, y0 - 0.26,
                f"[{r * 64}, {(r + 1) * 64})",
                ha="center", va="center", fontsize=8, color="#555")
        if r < P - 1:
            arrow(ax, (x + w + 0.02, y0 + h + 0.12),
                  (x + w + gap - 0.02, y0 + h + 0.12), lw=1.2)
            arrow(ax, (x + w + gap - 0.02, y0 + h + 0.40),
                  (x + w + 0.02, y0 + h + 0.40), lw=1.2)

    ax.annotate("halo\n(copia)", xy=(x0 + 0.09, y0 + h / 2),
                xytext=(x0 - 0.55, y0 + h + 0.55), fontsize=8.5,
                color=C_GHOST_EDGE, ha="center",
                arrowprops=dict(arrowstyle="-", color=C_GHOST_EDGE, lw=0.8))
    ax.text((x0 + xr) / 2, y0 + h + 0.95,
            "cada capa = NY × NZ celdas; solo se comunican las caras naranjas",
            ha="center", fontsize=9, color="#555")


    # borde periódico
    ax.text((x0 + xr) / 2, y0 - 1.15,
            "borde periódico: rank 3 se comunica con rank 0",
            ha="center", fontsize=9, color="#6a329f")

    # eje x
    arrow(ax, (x0, y0 - 0.60), (xr + 0.5, y0 - 0.60), color="#444", lw=1.1)
    ax.text(xr + 0.62, y0 - 0.60, "x", fontsize=11, va="center")
    ax.text((x0 + xr) / 2, y0 - 0.88,
            "N celdas globales (ejemplo N=256, P=4 → nx=64)",
            ha="center", fontsize=8.5, color="#555")

    clean(ax, (x0 - 0.7, xr + 0.95), (0.3, y0 + h + 1.35))
    ax.set_title("Descomposición del dominio entre procesos (P = 4)", fontsize=14,
                 weight="bold", pad=6)
    fig.savefig(os.path.join(OUT, "decomposition.png"), dpi=200,
                bbox_inches="tight")
    plt.close(fig)


# ---------------------------------------------------------------------------
def halo_exchange():
    fig, ax = plt.subplots(figsize=(12, 5.4))
    y0, h = 2.2, 2.0
    # rank r
    box(ax, 0.5, y0, 4.2, h, C_PROC, C_PROC_EDGE)
    box(ax, 0.5, y0, 0.55, h, C_GHOST, C_GHOST_EDGE, lw=0.9)
    box(ax, 4.15, y0, 0.55, h, C_GHOST, C_GHOST_EDGE, lw=0.9)
    ax.text(2.4, y0 + h / 2 + 0.12, "rank r", ha="center", fontsize=12)
    ax.text(2.4, y0 + h / 2 - 0.25, "interior (cómputo local)",
            ha="center", fontsize=8.5, color="#555")
    ax.text(0.78, y0 + h + 0.16, "ghost_hi", fontsize=8.5,
            color=C_GHOST_EDGE, ha="center")
    ax.text(4.42, y0 + h + 0.16, "face_hi", fontsize=8.5,
            color=C_GHOST_EDGE, ha="center")

    # rank r+1
    box(ax, 7.3, y0, 4.2, h, C_PROC, C_PROC_EDGE)
    box(ax, 7.3, y0, 0.55, h, C_GHOST, C_GHOST_EDGE, lw=0.9)
    box(ax, 10.95, y0, 0.55, h, C_GHOST, C_GHOST_EDGE, lw=0.9)
    ax.text(9.4, y0 + h / 2 + 0.12, "rank r+1", ha="center", fontsize=12)
    ax.text(9.4, y0 + h / 2 - 0.25, "interior (cómputo local)",
            ha="center", fontsize=8.5, color="#555")
    ax.text(7.58, y0 + h + 0.16, "ghost_lo", fontsize=8.5,
            color=C_GHOST_EDGE, ha="center")
    ax.text(11.22, y0 + h + 0.16, "face_lo", fontsize=8.5,
            color=C_GHOST_EDGE, ha="center")

    # flechas cortas
    arrow(ax, (4.75, y0 + h - 0.45), (7.70, y0 + h - 0.45), lw=2.0)
    ax.text(6.22, y0 + h - 0.30, "face_hi → ghost_lo",
            ha="center", fontsize=8.5, color=C_COMM,
            bbox=dict(facecolor="white", edgecolor="none", pad=1.0))
    arrow(ax, (7.70, y0 + 0.45), (4.75, y0 + 0.45), lw=2.0)
    ax.text(6.22, y0 + 0.22, "face_lo → ghost_hi",
            ha="center", fontsize=8.5, color=C_COMM,
            bbox=dict(facecolor="white", edgecolor="none", pad=1.0))

    # código
    code = ("bsp_put(right, face_hi, ghost_lo, 0, bytes);\n"
            "bsp_put(left,  face_lo, ghost_hi, 0, bytes);\n"
            "bsp_sync();")
    ax.text(6.0, 0.95, code, ha="center", va="center",
            family="monospace", fontsize=9,
            bbox=dict(facecolor="#f7f7f7", edgecolor="#999999",
                      boxstyle="round,pad=0.35"))

    ax.text(6.0, y0 + h + 0.85,
            "Payload por mensaje: una capa de (N+2) × (N+2) bytes",
            ha="center", fontsize=9.5, color="#555")
    ax.text(6.0, y0 + h + 1.20,
            "No se envía el volumen: solo una cara por vecino",
            ha="center", fontsize=11.5, weight="bold")

    clean(ax, (0.0, 12.0), (0.2, y0 + h + 1.40))
    ax.set_title("Intercambio de halos entre dos procesos", fontsize=14,
                 weight="bold", pad=4)
    fig.savefig(os.path.join(OUT, "halo_exchange.png"), dpi=200,
                bbox_inches="tight")
    plt.close(fig)


# ---------------------------------------------------------------------------
def communication_volume():
    fig, ax = plt.subplots(figsize=(11, 5.4))
    x0, y0, w, h = 2.4, 1.9, 7.0, 2.6

    # interior
    box(ax, x0, y0, w, h, C_COMP, C_PROC_EDGE, lw=1.2)
    # caras x (comunicadas)
    box(ax, x0 - 0.28, y0, 0.28, h, C_GHOST, C_GHOST_EDGE, lw=1.0)
    box(ax, x0 + w, y0, 0.28, h, C_GHOST, C_GHOST_EDGE, lw=1.0)
    # bordes y (copias locales)
    box(ax, x0, y0 + h, w, 0.22, C_LOCAL, C_LOCAL_EDGE, lw=1.0)
    box(ax, x0, y0 - 0.22, w, 0.22, C_LOCAL, C_LOCAL_EDGE, lw=1.0)

    ax.text(x0 + w / 2, y0 + h / 2 + 0.10,
            "interior: nx-2 capas", ha="center", fontsize=10, color="#1c3f5e")
    ax.text(x0 + w / 2, y0 + h / 2 - 0.30,
            "cómputo local (OpenMP)", ha="center", fontsize=8.5,
            color="#1c3f5e")

    # etiquetas caras x
    ax.annotate("cara x = 0\n(se envía)", xy=(x0 - 0.14, y0 + h * 0.55),
                xytext=(x0 - 1.55, y0 + h * 0.85), fontsize=9,
                color=C_GHOST_EDGE, ha="center",
                arrowprops=dict(arrowstyle="-", color=C_GHOST_EDGE, lw=0.8))
    ax.annotate("cara x = nx-1\n(se envía)", xy=(x0 + w + 0.14, y0 + h * 0.55),
                xytext=(x0 + w + 1.45, y0 + h * 0.85), fontsize=9,
                color=C_GHOST_EDGE, ha="center",
                arrowprops=dict(arrowstyle="-", color=C_GHOST_EDGE, lw=0.8))

    ax.text(x0 + w / 2, y0 + h + 0.36, "borde y: copia local",
            ha="center", fontsize=8.5, color=C_LOCAL_EDGE)
    ax.text(x0 + w / 2, y0 - 0.45, "borde y: copia local",
            ha="center", fontsize=8.5, color=C_LOCAL_EDGE)

    # dimensiones
    ax.annotate("", xy=(x0, y0 - 0.85), xytext=(x0 + w, y0 - 0.85),
                arrowprops=dict(arrowstyle="<->", color="#444", lw=1.0))
    ax.text(x0 + w / 2, y0 - 1.10, "nx capas (dirección de la partición)",
            ha="center", fontsize=8.5, color="#444")
    ax.text(x0 + w / 2, y0 - 1.55,
            "z: NZ capas (no dibujado); los bordes y/z son copias locales",
            ha="center", fontsize=8.5, color="#555")

    # leyenda
    handles = [
        Patch(facecolor=C_GHOST, edgecolor=C_GHOST_EDGE,
              label="2 caras x: bsp_put"),
        Patch(facecolor=C_LOCAL, edgecolor=C_LOCAL_EDGE,
              label="bordes y/z: copia local"),
        Patch(facecolor=C_COMP, edgecolor=C_PROC_EDGE,
              label="interior: cómputo local"),
    ]
    fig.legend(handles=handles, loc="lower center", ncol=3, frameon=False,
               fontsize=8.5)

    ax.text(x0 + w / 2, y0 + h + 1.15,
            "Qué se comunica en cada paso: 2 caras, no el volumen",
            ha="center", fontsize=13, weight="bold")
    ax.text(x0 + w / 2, y0 + h + 0.82,
            "cada cara es una capa de NY × NZ celdas",
            ha="center", fontsize=9.5, color="#555")

    clean(ax, (0.5, x0 + w + 1.9), (0.2, y0 + h + 1.35))
    fig.savefig(os.path.join(OUT, "communication_volume.png"), dpi=200,
                bbox_inches="tight")
    plt.close(fig)


# ---------------------------------------------------------------------------
def decomposition_uneven():
    fig, ax = plt.subplots(figsize=(12, 5.2))
    N, P = 10, 3
    base, rem = N // P, N % P
    nx = [base + (1 if r < rem else 0) for r in range(P)]
    x0 = [r * base + min(r, rem) for r in range(P)]

    scale, h, gap, y0 = 0.62, 1.5, 0.45, 1.7
    xs = []
    x = 0.5
    for r in range(P):
        w = nx[r] * scale
        xs.append((x, w))
        box(ax, x, y0, w, h, C_PROC, C_PROC_EDGE)
        box(ax, x, y0, 0.16, h, C_GHOST, C_GHOST_EDGE, lw=0.8)
        box(ax, x + w - 0.16, y0, 0.16, h, C_GHOST, C_GHOST_EDGE, lw=0.8)
        ax.text(x + w / 2, y0 + h / 2 + 0.12, f"rank {r}",
                ha="center", fontsize=11)
        ax.text(x + w / 2, y0 + h / 2 - 0.22, f"nx = {nx[r]}",
                ha="center", fontsize=9, color="#555")
        ax.text(x + w / 2, y0 - 0.28, f"x0 = {x0[r]}",
                ha="center", fontsize=8.5, color="#555")
        ax.text(x + w / 2, y0 - 0.55, f"[{x0[r]}, {x0[r] + nx[r]})",
                ha="center", fontsize=8.5, color="#555")
        if r < P - 1:
            xr = x + w
            arrow(ax, (xr + 0.02, y0 + h + 0.12),
                  (xr + gap - 0.02, y0 + h + 0.12), lw=1.2)
            arrow(ax, (xr + gap - 0.02, y0 + h + 0.40),
                  (xr + 0.02, y0 + h + 0.40), lw=1.2)
        x += w + gap

    ax.text(0.5, y0 + h + 0.95,
            f"N = {N}, P = {P}  →  base = {base}, resto = {rem}",
            fontsize=11.5, weight="bold")
    ax.text(0.5, y0 + h + 0.65,
            "el resto se reparte en los primeros rem ranks:  nx(r) = base + [r < rem]",
            fontsize=9, color="#555")
    ax.text(0.5, y0 + h + 1.25,
            "Reparto cuando N no es divisible por P",
            fontsize=13.5, weight="bold")
    ax.text((xs[0][0] + xs[-1][0] + xs[-1][1]) / 2, 0.55,
            "borde periódico: rank 2 se comunica con rank 0",
            ha="center", fontsize=9, color="#6a329f")

    clean(ax, (0.1, xs[-1][0] + xs[-1][1] + 0.3), (0.3, y0 + h + 1.45))
    fig.savefig(os.path.join(OUT, "decomposition_uneven.png"), dpi=200,
                bbox_inches="tight")
    plt.close(fig)


# ---------------------------------------------------------------------------
def iso(x, y, z):
    return (x + 0.55 * y, z + 0.32 * y)


def poly3(ax, pts, fc, ec, lw=1.0, alpha=1.0, ls="-", z=2):
    ax.add_patch(Polygon([iso(*p) for p in pts], closed=True, facecolor=fc,
                         edgecolor=ec, linewidth=lw, alpha=alpha,
                         linestyle=ls, zorder=z))


def communication_volume_3d():
    fig, ax = plt.subplots(figsize=(10.5, 6.2))
    nx, ny, nz = 4.0, 2.2, 2.2

    # volumen interior
    poly3(ax, [(0, 0, 0), (nx, 0, 0), (nx, 0, nz), (0, 0, nz)],
          C_COMP, C_PROC_EDGE, z=3)
    poly3(ax, [(0, 0, nz), (nx, 0, nz), (nx, ny, nz), (0, ny, nz)],
          "#7faac8", C_PROC_EDGE, z=2)
    poly3(ax, [(nx, 0, 0), (nx, ny, 0), (nx, ny, nz), (nx, 0, nz)],
          "#6f9bbd", C_PROC_EDGE, z=2)

    # caras x explotadas
    xl = -0.9
    poly3(ax, [(xl, 0, 0), (xl, ny, 0), (xl, ny, nz), (xl, 0, nz)],
          C_GHOST, C_GHOST_EDGE, z=4)
    xr = nx + 0.9
    poly3(ax, [(xr, 0, 0), (xr, ny, 0), (xr, ny, nz), (xr, 0, nz)],
          C_GHOST, C_GHOST_EDGE, z=4)

    # conectores
    for (x0_, y0_, z0_) in [(0, 0, 0), (0, ny, 0), (0, ny, nz), (0, 0, nz)]:
        ax.plot([iso(xl, y0_, z0_)[0], iso(x0_, y0_, z0_)[0]],
                [iso(xl, y0_, z0_)[1], iso(x0_, y0_, z0_)[1]],
                ls="--", color=C_GHOST_EDGE, lw=0.8, zorder=1)
    for (x0_, y0_, z0_) in [(nx, 0, 0), (nx, ny, 0), (nx, ny, nz),
                            (nx, 0, nz)]:
        ax.plot([iso(xr, y0_, z0_)[0], iso(x0_, y0_, z0_)[0]],
                [iso(xr, y0_, z0_)[1], iso(x0_, y0_, z0_)[1]],
                ls="--", color=C_GHOST_EDGE, lw=0.8, zorder=1)

    # etiquetas
    ax.text(*iso(nx / 2, 0, nz / 2), "interior\ncómputo local",
            ha="center", va="center", fontsize=9, color="#12395a")


    # dimensiones
    ax.annotate("", xy=iso(0, 0, -0.2), xytext=iso(nx, 0, -0.2),
                arrowprops=dict(arrowstyle="<->", color="#444", lw=1.0))
    ax.text(*iso(nx / 2, 0, -0.55), "nx", ha="center", fontsize=9,
            color="#444", zorder=10)


    handles = [
        Patch(facecolor=C_GHOST, edgecolor=C_GHOST_EDGE,
              label="caras x (se envían): bsp_put"),
        Patch(facecolor=C_COMP, edgecolor=C_PROC_EDGE,
              label="volumen interior: cómputo local"),
    ]
    fig.legend(handles=handles, loc="lower center", ncol=2, frameon=False,
               fontsize=8.5)

    clean(ax, (-2.4, 7.6), (-0.8, 3.6))
    fig.savefig(os.path.join(OUT, "communication_volume_3d.png"), dpi=200,
                bbox_inches="tight")
    plt.close(fig)


# ---------------------------------------------------------------------------
def hybrid_threads():
    fig, ax = plt.subplots(figsize=(10, 5.2))
    P, T = 4, 4
    x0, y0 = 0.6, 0.8
    bw, bh = 2.0, 3.4
    gap = 0.28
    shades = ["#cfe2f3", "#b6d7e8", "#9fc5e8", "#8ab6d6"]
    for p in range(P):
        x = x0 + p * (bw + gap)
        box(ax, x, y0, 0.16, bh, C_GHOST, C_GHOST_EDGE, lw=0.7)
        box(ax, x + bw - 0.16, y0, 0.16, bh, C_GHOST, C_GHOST_EDGE, lw=0.7)
        th = bh / T
        for t in range(T):
            box(ax, x + 0.16, y0 + t * th, bw - 0.32, th, shades[t],
                "#7a9cbb", lw=0.5)
        ax.text(x + bw / 2, y0 + bh + 0.18, f"proceso {p}",
                ha="center", fontsize=9.5)
        if p == 0:
            for t in range(T):
                ax.text(x - 0.12, y0 + (t + 0.5) * th, f"h{t}",
                        ha="right", va="center", fontsize=7.5, color="#555")

    ax.text(x0, y0 + bh + 0.75, "P = 4 procesos × T = 4 hilos",
            fontsize=13, weight="bold")
    ax.text(x0, y0 + bh + 0.42,
            "OpenMP reparte el interior; el hilo maestro hace bsp_put y bsp_sync",
            fontsize=9, color="#555")

    arrow(ax, (x0, y0 - 0.35), (x0 + P * (bw + gap) - gap, y0 - 0.35),
          color="#444", lw=1.0)
    ax.text(x0 + P * (bw + gap) - gap + 0.12, y0 - 0.35, "x",
            fontsize=10, va="center")

    clean(ax, (0.0, x0 + P * (bw + gap) + 0.2), (0.2, y0 + bh + 0.95))
    ax.set_title("Paralelismo híbrido: procesos e hilos",
                 fontsize=14, weight="bold", pad=4)
    fig.savefig(os.path.join(OUT, "hybrid_threads.png"), dpi=200,
                bbox_inches="tight")
    plt.close(fig)


def main():
    decomposition()
    decomposition_uneven()
    halo_exchange()
    communication_volume()
    communication_volume_3d()
    hybrid_threads()
    print("Figuras en", OUT)


if __name__ == "__main__":
    main()
