# hpc-lab2 — Expansión de un gas 3D con autómata celular híbrido (BSPlib + OpenMP)

**Autores:** Gonzalo Ahumada y Daniel Muñoz.

Simulación de un gas que se expande desde una esfera densa o desde una cámara
de dos presiones en una malla cúbica, usando un **autómata celular de gas
reticular (LGA)** en 3D y paralelismo híbrido:

- **Procesos:** BSPlib (BSPonMPI) con descomposición de dominio 1D y
  superpasos BSP (`bsp_push_reg` + `bsp_put` + `bsp_sync`).
- **Hilos:** OpenMP sobre el interior de cada subdominio.
- **Salida:** archivos `.vti` para ParaView con la densidad por celda.
- **Esquemas:** descomposición, halos y qué se comunica en
  `docs/figs/decomposition.png`, `halo_exchange.png`,
  `communication_volume.png` y `hybrid_threads.png`.

Documentos:

- [`docs/informe.pdf`](docs/informe.pdf) — informe con modelo, validación y
  resultados.
- [`docs/presentacion_clase.pdf`](docs/presentacion_clase.pdf) — **presentación
  para la clase** (18 páginas): solo lo que pidió el profesor, es decir el
  sistema, las reglas del autómata y la paralelización híbrida. **Es la
  versión que se presenta.**
- [`docs/presentacion.pdf`](docs/presentacion.pdf) — presentación completa de la
  propuesta (29 páginas: además incluye validación, experimentos, riesgos y
  alternativas). Queda como referencia.

> **Importante: las dos presentaciones son archivos independientes.** Las
> slides de `presentacion_clase.tex` se copiaron desde `presentacion.tex`; no
> están enlazadas. Un cambio en uno **no** aparece en el otro. Acuerdo: los
> cambios para la clase se hacen **solo en `presentacion_clase.tex`**.
- [`docs/JUSTIFICACION_GAS_IDEAL.md`](docs/JUSTIFICACION_GAS_IDEAL.md) —
  justificación del gas ideal y mapeo a neón/helio.
- [`AI_HANDOFF.md`](AI_HANDOFF.md) — contexto completo para revisión por otra
  IA.

## Modelo

- Malla `N × N × N`. Cada celda guarda una máscara de 6 bits: una partícula
  por dirección (±x, ±y, ±z), con exclusión (máx. 1 por dirección).
- **Colisión:** tabla precomputada que mapea cada estado a otro con la misma
  masa y el mismo momento. Hay 8 tablas y se elige una por celda con un hash,
  para romper simetrías sin romper los invariantes.
- **Streaming:** cada partícula se mueve a la celda vecina en su dirección.
- **Condición inicial:** dos geometrías:
  - `sphere` (default): esfera central densa; vacío afuera.
  - `chamber`: cámara de alta presión (`x < N/4`) contra una de baja presión,
    con razón de presiones configurable y gradiente suavizable (`RAMP`).
- **Bordes:** periódicos en x (anillo entre procesos) y en y/z (copia local).
- **Invariantes:** masa y momento total se conservan exactamente.

## Gas y unidades físicas

La dinámica del LGA no depende de la especie: es un gas reticular idealizado.
Para reportar magnitudes físicas se elige un gas concreto y se mapean las
unidades lattice a SI. El caso base es **neón** (monoatómico, `γ = 5/3`,
`c_s ≈ 352 m/s` isotérmica a 300 K); también se soporta **helio**.

Variables de entorno (el programa imprime el mapeo al inicio):

| Variable | Default | Significado |
|---|---|---|
| `GAS` | `neon` | `neon` o `helium` |
| `DX_MM` | `1` | tamaño de celda en mm; `0` = 1 partícula lattice = 1 molécula |
| `T0_K` | `300` | temperatura de referencia |
| `P0_ATM` | `1` | presión de referencia |
| `GEOM` | `sphere` | `sphere` o `chamber` |
| `RHO_H` | `3.0` | densidad alta (partículas/celda) |
| `RHO_L` | `0` (`0.3` en chamber) | densidad baja/fondo |
| `P_RATIO` | — | si se define, `RHO_L = RHO_H / P_RATIO` |
| `IFACE_FRAC` | `0.25` | posición de la interfaz (fracción de `N`) |
| `RAMP` | `0` | ancho del gradiente en celdas (0 = escalón) |
| `PROFILE_CSV` | — | si se define, escribe el perfil $\rho(x,t)$ y métricas por paso |
| `PROFILE_EVERY` | `1` | cada cuántos pasos muestrear el perfil |

Ejemplo con neón, celdas de 1 mm:

```bash
GAS=neon DX_MM=1 OMP_NUM_THREADS=4 mpirun -np 4 --bind-to none \
    ./build/gas 256 1000 4 -1
```

Salida típica:

```
gas=neon m=3.35092e-26 kg gamma=1.66667
T0=300 K p0=101325 Pa dx=0.001 m
cs_iso=351.577 m/s cs_ad=453.884 m/s dt=1.64217e-06 s v_lat=608.949 m/s
rho_mass=0.81974 kg/m3 molecules_per_particle=8.15438e+15
domain=0.256 m mass_total=0.000112677 kg
```

La justificación completa (gas ideal, gas noble, corrección del virial,
mapeos macroscópico/atómico y conversión del VTK) está en
[`docs/JUSTIFICACION_GAS_IDEAL.md`](docs/JUSTIFICACION_GAS_IDEAL.md).

## Paralelización

```
rank 0        rank 1        rank 2        rank 3
[ nx capas ] [ nx capas ] [ nx capas ] [ nx capas ]
     ↕ halo        ↕ halo        ↕ halo
        bsp_put + bsp_sync  (superpaso BSP)

Dentro de cada rank:
  #pragma omp parallel for collapse(3)   // collide y stream
```

- `N` no necesita ser divisible por `P`: se reparte `N/P` y el resto.
- El halo se intercambia con dos `bsp_put` por paso (una cara a cada vecino) y
  se materializa en `bsp_sync`.
- Los halos en y/z son periódicos y se copian localmente, sin comunicación.

## Requisitos

- C++17, OpenMP, MPI (OpenMPI o MPICH) y **BSPonMPI**.
- BSPonMPI no está en los repositorios: hay que compilarlo. El
  [`docker/Dockerfile`](docker/Dockerfile) lo hace solo.
- Con OpenMPI 4.x es obligatorio compilar con `-DOMPI_SKIP_MPICXX` (el
  `CMakeLists.txt` ya lo define).

### Opción A: Docker (recomendado)

```bash
docker build -t hpc-lab2:latest -f docker/Dockerfile .

# compilar y correr 4 procesos × 4 hilos
./scripts/docker_shell.sh bash -lc './scripts/build.sh && ./scripts/run.sh 4 4 256 500'
```

### Opción B: máquina con MPI + BSPonMPI instalados

```bash
# Ajustar BSPONMPI_LIBRARY / BSPONMPI_INCLUDE_DIR si no están en /usr/local
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j

OMP_NUM_THREADS=4 mpirun -np 4 --bind-to none \
    ./build/gas 256 1000 4 0 gas
```

`--bind-to none` es importante: si OpenMPI ata cada rank a un núcleo, los
hilos de OpenMP se pelean por él.

## Uso

```
./gas [N] [pasos] [hilos] [cada_cuanto_vtk] [prefijo]

  N            lado de la malla (default 256)
  pasos        pasos de simulación (default 1000)
  hilos        hilos OpenMP por proceso (default 1)
  cada_cuanto  0 = VTK solo al inicio y al final; >0 = cada n pasos;
               <0 = sin VTK (default 0)
  prefijo      prefijo de los .vti (default gas)
```

Ejemplo con frames cada 100 pasos:

```bash
OMP_NUM_THREADS=4 mpirun -np 4 --bind-to none ./build/gas 256 1000 4 100 gas
```

Cámara de dos presiones (10:1) con gradiente suave de 8 celdas:

```bash
GEOM=chamber P_RATIO=10 RAMP=8 OMP_NUM_THREADS=4 \
  mpirun -np 4 --bind-to none ./build/gas 256 1000 4 100 chamber
```

## Salida

- `gas_00000.vti`, `gas_00100.vti`, … : densidad (0..6) en formato VTK
  ImageData, abribles en ParaView.
- Resumen por stdout (solo rank 0):

```
N=256 P=4 T=4 steps=1000 geom=sphere
mass0=412365 momentum0=(-273,82,-506)
mass1=412365 (delta=0)
momentum1=(-273,82,-506)
time_loop=... s compute_max=... s comm_max=... s
```

## Evolución temporal

El diagnóstico `PROFILE_CSV` vuelca el perfil $\rho(x,t)$ y métricas globales
por paso. Ejemplo para la cámara:

```bash
GEOM=chamber P_RATIO=10 RAMP=8 \
PROFILE_CSV=data/figs/chamber_profile.csv PROFILE_EVERY=1 \
OMP_NUM_THREADS=4 mpirun -np 4 --bind-to none ./build/gas 128 150 4 -1

python3 scripts/plot_evolution.py
```

Genera `docs/figs/spacetime_chamber.png` (diagrama espacio-temporal),
`docs/figs/evolution_chamber.png` (métricas) y el GIF
`docs/figs/evolution_chamber.gif`. El muestreo agrega un `bsp_sync` por paso:
usarlo para diagnóstico, no para medir rendimiento.

## Esquemas de comunicación y paralelización

```bash
python3 scripts/make_comm_figures.py
```

Genera:

- `docs/figs/decomposition.png`: reparto 1D entre procesos y halos.
- `docs/figs/decomposition_uneven.png`: caso $N$ no divisible por $P$.
- `docs/figs/communication_volume.png`: qué se comunica por paso (2D).
- `docs/figs/communication_volume_3d.png`: vista 3D del subdominio local.
- `docs/figs/halo_exchange.png`: detalle de `bsp_put`/`bsp_sync`.
- `docs/figs/hybrid_threads.png`: procesos × hilos OpenMP.

## Validación

- Masa y momento conservados exactamente en `P = 1, 2, 3, 4`, con `N`
  divisible y no divisible por `P`, en ambas geometrías.
- La suma de la densidad del `.vti` coincide con la masa global.
- `P=1` vs `P=4` produce archivos `.vti` **idénticos byte a byte** (esfera y
  cámara), lo que confirma que los halos y la descomposición no alteran la
  física.
- Cámara: `p_H = 1 atm`, `p_L = 0.1 atm` (razón 10:1); la densidad inicial
  muestra los dos niveles y la expansión rompe la interfaz.
- Test de halos: `tests/halo.cpp` (0 errores).
- Test híbrido: `tests/hello.cpp` (4 procesos × 4 hilos).

## Rendimiento medido (referencia)

Laptop AMD Ryzen 7 PRO 4750U (8 núcleos / 16 hilos), `N=256`, 200 pasos,
sin VTK, dentro de Docker:

| P | T | hilos totales | tiempo (s) | speedup vs 1×1 |
|---|---|---:|---:|---:|
| 1 | 1 | 1 | 21.27 | 1.00 |
| 1 | 2 | 2 | 11.11 | 1.91 |
| 1 | 4 | 4 | 6.19 | 3.44 |
| 1 | 8 | 8 | 3.72 | 5.72 |
| 1 | 16 | 16 | 4.04 | 5.26 |
| 2 | 4 | 8 | 6.40 | 3.32 |
| 4 | 2 | 8 | 6.76 | 3.15 |
| 4 | 4 | 16 | 4.06 | 5.24 |
| 8 | 2 | 16 | 4.18 | 5.09 |

Lectura: el problema es **memory-bound** y en un solo nodo de 8 núcleos el
mejor resultado es OpenMP puro (1×8). Los procesos no mejoran el tiempo en un
nodo porque el costo de `bsp_sync`/alltoall de BSPonMPI compensa la división
del trabajo. La capa de procesos se justifica al cruzar nodos o con mallas más
grandes. Para el objetivo 4×8 hace falta un nodo de 32 hilos.

## Estructura

```
.
├── CMakeLists.txt
├── src/gas.cpp          # autómata + driver BSP/OpenMP + VTK
├── tests/hello.cpp      # smoke test BSPlib + OpenMP
├── tests/halo.cpp       # test de intercambio de halos
├── docker/Dockerfile    # Ubuntu + OpenMPI + BSPonMPI
├── scripts/             # build, run, scaling, figuras, shell docker
│   ├── make_comm_figures.py
│   ├── make_figures.py
│   └── plot_evolution.py
├── data/                # salida .vti y CSV (ignorada por git)
├── results/             # CSV de scaling (ignorado por git)
├── docs/
│   ├── informe.tex / informe.pdf
│   ├── presentacion_clase.tex / .pdf   # versión para la clase (se presenta)
│   ├── presentacion.tex / presentacion.pdf   # propuesta completa (referencia)
│   ├── JUSTIFICACION_GAS_IDEAL.md
│   └── figs/            # figuras del informe y la presentación
├── ejemplos/            # material original del curso
└── AI_HANDOFF.md        # contexto completo para revisión por otra IA
```

## Limitaciones y trabajo futuro

- Modelo LGA 6-velocidades: anisotropía cúbica y poca termalización; se puede
  mejorar con partícula en reposo, más direcciones o reglas de colisión
  específicas.
- La cámara ya está implementada (escalón o `RAMP`); falta la opción de paredes
  reflectantes para un recipiente cerrado (hoy todo es periódico).
- Es un modelo isotérmico: no captura el enfriamiento adiabático.
- El VTK solo exporta densidad; se puede agregar momento/velocidad.
- Optimizaciones: `-march=native`, eliminar branches en el streaming, tablas
  de colisión sin hash, comunicación solapada (fuera del modelo BSP puro).
- Estudio de scaling en un nodo de 32 hilos para el 4×8 real.
