# AI_HANDOFF — hpc-lab2

Documento de traspaso para que otra IA (o una persona) pueda entender, revisar,
corregir o extender este proyecto sin reconstruir el contexto.

- **Repo:** `/home/demm/Projects/usach/hpc-lab2`
- **Autores:** Gonzalo Ahumada y Daniel Muñoz.
- **Curso:** Laboratorio 2 — High Performance Computing, Universidad de
  Santiago de Chile.
- **Última actualización:** 2026-11-26.
- **Estado:** implementación, validación, informe y presentación terminados.
  Quedan pendientes experimentales (nodo de 32 hilos, paredes reflectantes,
  modelo térmico) que no bloquean la entrega.

> **Advertencia importante para quien revise:**
> - El **informe** (`docs/informe.pdf`) **sí incluye resultados**.
> - La **presentación que se expone en clase** es `docs/presentacion_clase.pdf`
>   (fuente `docs/presentacion_clase.tex`, 15 páginas): solo sistema, reglas
>   del CA y paralelización híbrida, que es lo que pidió el profesor. Sus
>   slides son **copias** de `presentacion.tex`, no están enlazadas: editar
>   **solo** `presentacion_clase.tex` para la clase.
> - La **presentación completa** (`docs/presentacion.pdf`) es una **propuesta**: no
>   incluye resultados, solo problema, modelo, decisiones, plan de validación,
>   plan de experimentos, riesgos y alternativas. No agregar métricas de
>   corridas a la presentación sin consultar.

---

## 0. Cómo usar este documento

1. Leer §1 y §2 para el contexto y el estado.
2. Seguir §3 para compilar y correr.
3. Leer §4 y §5 antes de tocar la comunicación (BSPonMPI tiene varias trampas).
4. Leer §6 y §7 para el modelo y la paralelización.
5. Usar §8 como mapa del código.
6. Revisar §9 y §10 para la evidencia y los números.
7. Usar §14 para reproducir todo desde cero.
8. Si se cambia el modelo o la comunicación, **volver a validar** (§9).

---

## 1. Contexto del curso y objetivo

Proyecto de HPC que simula la **expansión de un gas en 3D** con un **autómata
celular de gas reticular (LGA)** y **paralelismo híbrido**:

- **Procesos:** BSPlib sobre MPI (BSPonMPI), descomposición de dominio 1D y
  superpasos BSP.
- **Hilos:** OpenMP dentro de cada subdominio.
- **Salida:** VTK `.vti` para ParaView y CSV de diagnóstico temporal.

Requisitos originales del profesor:

- C++.
- BSP para los procesos y OpenMP/OpenMPI para los hilos (la frase original era
  ambigua; la implementación real usa **BSPlib/BSPonMPI + OpenMP**).
- Ejecución exigente pero posible en un PC; referencia **4 procesos × 8 hilos**.
- Entregables: PDF con la descripción del problema, presentación Beamer con
  reglas del CA y esquema de paralelización, y archivos fuente.

La documentación original del curso está en `ejemplos/`:

- `Directivas BSPonMPI.pdf`: en realidad es el manual del **Oxford BSP
  Toolset**, no de BSPonMPI. Contiene comandos y colectivas que BSPonMPI no
  implementa.
- `Instalar_BsponMpi.pdf`: instalación de BSPonMPI 0.2 con autotools,
  `mpiCC`, `-lbsponmpi` y `LD_LIBRARY_PATH`.
- `Ejemplo2.cc.txt`, `Ejemplo2_v2.cc.txt`, `Ejemplo3.cc.txt`,
  `Ejemplo3_v3.txt`: ejemplos de BSPlib. Algunos tienen errores didácticos
  (serialización de punteros, buffers fijos, falta de llaves).
- Los textos extraídos de los PDFs (`directivas.txt`, `instalar.txt`) no se
  versionan; se pueden regenerar con `pdftotext` si hacen falta.

---

## 2. Estado actual

### 2.1 Funciona

- LGA 3D de **6 direcciones**, sin partícula en reposo, con exclusión.
- Colisión por **8 tablas** que conservan masa y momento exactamente.
- Streaming gather; formulación segura para hilos.
- Dos geometrías iniciales:
  - **`sphere`**: esfera densa (default).
  - **`chamber`**: cámara de dos presiones con gradiente configurable.
- Descomposición **1D en x** entre procesos, con halos:
  - `bsp_push_reg` + `bsp_put` + `bsp_sync`.
  - Soporta `N` no divisible por `P`.
- OpenMP en colisión, streaming, métricas y densidad.
- Salida VTK `.vti` (densidad por celda) escrita por rank 0.
- Diagnóstico temporal `PROFILE_CSV`: perfil $\rho(x,t)$ y métricas por paso.
- Mapeo físico a **neón** o **helio** (`GAS`, `DX_MM`, `T0_K`, `P0_ATM`).
- Informe PDF y presentación Beamer compilados.

### 2.2 Validado

- Masa y momento conservados exactamente en `P = 1, 2, 3, 4`, con `N`
  divisible y no divisible por `P`, en esfera y cámara.
- `P=1` vs `P=4` produce archivos `.vti` **idénticos byte a byte**:
  - esfera, `N=64`, 50 pasos;
  - cámara 10:1 con `RAMP=8`, `N=64`, 20 pasos.
- La suma de densidad del `.vti` coincide con la masa global.
- Test de halos (`tests/halo.cpp`): 0 errores en 4 ranks.
- Test híbrido (`tests/hello.cpp`): 4 procesos × 4 hilos.

### 2.3 Pendiente / no bloqueante

- Corrida **4×8 real** en un nodo de 32 hilos (la laptop es 8C/16T).
- Paredes reflectantes (hoy todo es periódico).
- Modelo térmico para expansión adiabática (hoy isotérmico).
- Exportar presión/velocidad al VTK (hoy la presión se deriva de la densidad).
- Optimizaciones de rendimiento y tests automáticos CTest.

### 2.4 Estado del repositorio

- No se hizo commit desde esta sesión; el usuario maneja git.
- `git status` esperado: `.gitignore` y `README.md` modificados; el resto
  (`AI_HANDOFF.md`, `CMakeLists.txt`, `docker/`, `docs/`, `scripts/`, `src/`,
  `tests/`, `ejemplos/`) sin trackear.
- `ejemplos/` es material del curso; no borrarlo.
- `data/`, `results/`, `build/` están ignorados por git.
- `./gas` en la raíz era un binario residual de una compilación manual; fue
  eliminado y está ignorado.
- Los textos extraídos de los PDFs y las figuras auxiliares no usadas están
  en `.gitignore`; no se versionan.
- Los auxiliares de LaTeX (`docs/*.aux`, `.log`, `.out`, `.toc`, `.nav`,
  `.snm`, `.vrb`, `.fls`, `.fdb_latexmk`) son generados; se pueden limpiar con
  `latexmk -c` en `docs/`.

---

## 3. Quick start

### 3.1 Docker (recomendado)

```bash
docker build -t hpc-lab2:latest -f docker/Dockerfile .

./scripts/docker_shell.sh bash -lc './scripts/build.sh'

# test híbrido y de halos
./scripts/docker_shell.sh bash -lc \
  'mpirun -np 4 --bind-to none ./build/hello 4'
./scripts/docker_shell.sh bash -lc \
  'mpirun -np 4 --bind-to none ./build/halo_test 4'

# corrida corta de la cámara
./scripts/docker_shell.sh bash -lc \
  'GEOM=chamber P_RATIO=10 RAMP=8 OMP_NUM_THREADS=4 \
   mpirun -np 4 --bind-to none ./build/gas 128 100 4 -1'
```

### 3.2 Máquina con MPI + BSPonMPI nativos

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j

OMP_NUM_THREADS=4 mpirun -np 4 --bind-to none \
  ./build/gas 256 1000 4 0 gas
```

### 3.3 Documentos

```bash
cd docs
latexmk -pdf informe.tex
latexmk -pdf presentacion.tex
latexmk -pdf presentacion_clase.tex
```

---

## 4. Entorno y build

### 4.1 Máquina de desarrollo

- AMD Ryzen 7 PRO 4750U: **8 núcleos / 16 hilos**, 7.1 GB RAM.
- Linux, gcc 16, Docker disponible, **sin sudo**, **sin MPI ni BSPonMPI
  nativos**.
- Todo se validó dentro de Docker.

### 4.2 Imagen Docker

`docker/Dockerfile` construye `hpc-lab2:latest` con:

- Ubuntu 24.04.
- `build-essential`, `cmake`.
- OpenMPI **4.1.6** (`openmpi-bin`, `libopenmpi-dev`).
- **BSPonMPI 0.3** descargado de SourceForge y compilado a mano.

BSPonMPI 0.3 usa un `SConstruct` con sintaxis **Python 2**, que SCons moderno
no parsea. Por eso el Dockerfile compila a mano:

```bash
mpicc -O2 -Iinclude -c src/*.c -o obj/...
ar rcs lib/libbsponmpi.a obj/*.o
cp include/*.h /usr/local/include/
cp lib/libbsponmpi.a /usr/local/lib/
```

La librería queda estática; no se necesita `LD_LIBRARY_PATH`. La guía del
curso describe BSPonMPI 0.2 (compartida), así que en ese entorno puede hacer
falta `LD_LIBRARY_PATH` y enlazar `-lbsponmpi`.

### 4.3 CMake

`CMakeLists.txt`:

- Busca `bsp.h` y `libbsponmpi`.
- Enlaza `MPI::MPI_CXX` y `OpenMP::OpenMP_CXX`.
- Define `OMPI_SKIP_MPICXX`.

**Clave:** `-DOMPI_SKIP_MPICXX` es obligatorio con OpenMPI 4.x. Sin eso,
OpenMPI incluye sus bindings C++ obsoletos dentro del
`extern "C" { #include <bsp.h> }` y falla con GCC 14.

### 4.4 Ejecución

```bash
OMP_NUM_THREADS=4 mpirun -np 4 --bind-to none ./build/gas 256 1000 4 0 gas
```

- `--bind-to none` evita que OpenMPI ate cada rank a un núcleo y ahogue a los
  hilos de OpenMP.
- En Docker como root hace falta `--allow-run-as-root`; los scripts usan
  `--user $(id -u):$(id -g)` para no dejar archivos root-owned.
- `bsp_begin(bsp_nprocs())` toma el tamaño de `MPI_COMM_WORLD`, fijado por
  `mpirun -np P`.

---

## 5. BSPonMPI: API real y trampas

BSPonMPI implementa BSPlib sobre MPI. No es el Oxford BSP Toolset, aunque la
API sea compatible en su mayor parte.

### 5.1 Arranque

```cpp
extern "C" {
#include <bsp.h>
}

void bsp_main() {
    bsp_begin(bsp_nprocs());
    // ...
    bsp_end();
}

int main(int argc, char **argv) {
    bsp_init(bsp_main, argc, argv);
    bsp_main();
}
```

- `bsp_init` es **obligatorio** (única desviación declarada del estándar).
- `bsp_init` hace `MPI_Init`, setea `nprocs = MPI_Comm_size`.
- En ranks distintos de 0, `bsp_init` llama a `bsp_main` y termina; en rank 0
  vuelve a `main`, que debe llamar a `bsp_main`.
- `bsp_begin(bsp_nprocs())` usa todos los ranks de `mpirun`.
- `bsp_begin(k)` con `k < P` deja fuera a los ranks sobrantes; para el proyecto
  se usan todos.

### 5.2 Funciones disponibles

El `bsp.h` de BSPonMPI 0.3 exporta:

```
bsp_abort bsp_begin bsp_end bsp_init bsp_nprocs bsp_pid bsp_time bsp_sync
bsp_push_reg bsp_pop_reg bsp_put bsp_get bsp_hpput bsp_hpget
bsp_send bsp_qsize bsp_get_tag bsp_move bsp_hpmove bsp_set_tagsize bsp_config
```

**No** existen `bsp_bcast`, `bsp_fold`, `bsp_l`, `bsp_s`, `bsp_g` ni los
comandos `bspcc`/`bsprun`/`bsparch` del manual de Oxford. Las reducciones se
hacen a mano con `bsp_send` o `bsp_put`.

### 5.3 Firmas importantes

Usa `size_t` en varios lugares, no `int`:

```c
void bsp_push_reg(const void *, size_t);
void bsp_put(int, const void *, void *, long int, size_t);
void bsp_send(int, const void *, const void *, size_t);
void bsp_qsize(int *, size_t *);
void bsp_get_tag(int *, void *);
void bsp_move(void *, size_t);
void bsp_set_tagsize(size_t *);
```

### 5.4 Quirk crítico: `bsp_get_tag` devuelve payload + tag

En BSPonMPI, `bsp_send` guarda `element.size = payload_nbytes + tag_size`, y
`bsp_get_tag` devuelve ese `element.size`, **no** el tamaño del payload. El
código lo compensa:

```cpp
constexpr int TAG_BYTES = sizeof(int);
int size = -1, tag = -1;
bsp_get_tag(&size, &tag);
if (size <= 0) break;
const int payload = size - TAG_BYTES;   // necesario
bsp_move(dst, payload);                 // pasar payload, no `size`
```

Si se pasa `size` a `bsp_move`, se copian 4 bytes de más y se corrompe el
destino. Esto está resuelto en `gather_stats`, `sample_profile` y
`gather_and_write`.

### 5.5 Semántica de superpaso

- `bsp_put`, `bsp_get`, `bsp_send` son **diferidos**: se materializan en el
  próximo `bsp_sync`.
- `bsp_push_reg`/`bsp_pop_reg` toman efecto en el próximo `bsp_sync`.
- El receptor no puede leer el halo hasta después del `bsp_sync`.
- `bsp_send` copia tag+payload; el orden de mensajes no está garantizado:
  despachar por tag.
- `bsp_time()` devuelve segundos desde `bsp_begin`.
- Solo rank 0 debería hacer I/O de archivos de forma confiable.

### 5.6 Otras notas

- BSPonMPI no es thread-safe: todas las llamadas desde el hilo maestro, fuera
  de regiones OpenMP.
- `bsp_put`/`bsp_send` no solapan comunicación con cómputo (modelo BSP puro).
- El `bsp_sync` es una barrera global; su costo crece con `P`.

---

## 6. Modelo físico y autómata

### 6.1 Gas ideal

Supuestos:

- Partículas puntuales, sin fuerzas de largo alcance.
- Colisiones binarias elásticas y caos molecular.
- Gas noble monoatómico: capa cerrada, $\gamma=5/3$, sin modos internos.
- A 1 atm y 300 K:
  $n=2.446\times10^{25}\,\text{m}^{-3}$ y
  $B_2(\text{Ne})\approx 11.4\,\text{cm}^3/\text{mol}$,
  por lo que $Z=1+B_2n\approx1.00046$; corrección del **0.05%**.

**Limitación declarada:** el autómata es **isotérmico** ($p=\rho/3$,
$c_s^2=1/3$). No captura el enfriamiento adiabático. Se justifica porque el
foco es HPC y morfología, no termodinámica fina.

### 6.2 Mapeo a neón/helio

Se fija $T_0=300$ K, $p_0=1$ atm y $\rho_0=3$ partículas/celda como estado de
referencia:

```text
n0        = p0 / (kB T0)
rho_mass0 = n0 * m
cs_iso    = sqrt(kB T0 / m)
dt        = (1/sqrt(3)) * dx / cs_iso
p [Pa]    = (rho/3) * p0
```

| Magnitud | Neón | Helio |
|---|---:|---:|
| `m` | $3.3509\times10^{-26}$ kg | $6.6465\times10^{-27}$ kg |
| `c_s` isotérmica | 351.6 m/s | 789.4 m/s |
| `c_s` adiabática | 453.9 m/s | 1019.1 m/s |
| densidad a 1 atm | 0.8197 kg/m³ | 0.1626 kg/m³ |
| `Δt` con `Δx=1 mm` | 1.642 µs | 0.731 µs |

Con `DX_MM=1` (1 mm) y `N=256`, cada partícula lattice representa
$N_c\approx8.15\times10^{15}$ moléculas; el dominio mide 25.6 cm. La masa total
de la esfera es 0.1127 g de neón; la cámara 10:1 tiene 4.47 g (incluye la
cámara de baja presión).

Con `DX_MM=0` se fuerza **1 partícula = 1 molécula**:
$\Delta x\approx4.97$ nm, dominio $\approx1.27\,\mu$m y
$\Delta t\approx8.16$ ps (Ne) / 3.63 ps (He).

Detalle completo: `docs/JUSTIFICACION_GAS_IDEAL.md`.

### 6.3 Estado del LGA

- Malla cúbica $N\times N\times N$.
- Cada celda: `uint8_t` con 6 bits, una partícula por dirección $\pm x,\pm y,
  \pm z$; exclusión.
- 64 estados posibles.

### 6.4 Colisión

Se precomputan 8 tablas `g_table[8][64]`. Para cada estado:

1. Calcular masa $M(s)$ y momento $\vec P(s)$.
2. Agrupar por $(M,P_x,P_y,P_z)$.
3. Mapear cada estado a otro elegido al azar dentro de su grupo.

Como el mapeo no sale del grupo, se conservan masa y momento exactamente. La
tabla se genera con `std::mt19937` y semilla fija `0x5eed1234`, iterando un
`std::map` (orden determinista). **Todos los ranks deben generar la misma
tabla**; no cambiar a `std::unordered_map` sin sincronizar.

La tabla se elige por celda con `coord_hash(x_global, y, z, step) & 7`.

### 6.5 Streaming

Formulación *gather*:

```text
s^{t+1}_{x,y,z} =
    (s^t_{x-1,y,z} & B_{+x}) |
    (s^t_{x+1,y,z} & B_{-x}) |
    (s^t_{x,y-1,z} & B_{+y}) |
    (s^t_{x,y+1,z} & B_{-y}) |
    (s^t_{x,y,z-1} & B_{+z}) |
    (s^t_{x,y,z+1} & B_{-z})
```

Cada hilo escribe celdas distintas y solo lee de `tmp`: no hay carreras.

### 6.6 Condiciones iniciales

**`sphere` (default):**

- Esfera central, radio $R=\max(6,N/8)$.
- $\rho_H=3$ adentro, $\rho_L=0$ afuera.
- Cada dirección se ocupa con probabilidad $\rho/6$.

**`chamber`:**

- $\rho_H=3$ ($p_H=1$ atm) para $x<0.25N$.
- $\rho_L=0.3$ ($p_L=0.1$ atm) para $x\ge0.25N$.
- Razón $10{:}1$; `RAMP` suaviza la interfaz con un gradiente lineal.
- Con `Δx=1 mm` y `RAMP=8`, $|\nabla p|\approx0.1125$ atm/mm.

La densidad objetivo se muestrea dirección por dirección usando un hash de
**coordenadas globales**, así que el estado inicial no depende de `P`.

### 6.7 Bordes

- **x:** periódico, implementado como anillo entre ranks.
- **y,z:** periódico, copiando localmente los planos opuestos en los halos.

No hay paredes reflectantes. Para un recipiente cerrado habría que reflejar
las partículas en las caras (invertir el bit normal).

### 6.8 Invariantes

- Masa total $=\sum \operatorname{popcount}(s)$.
- Momento total $=\sum \vec P(s)$.
- Ambos se conservan exactamente: estado entero, reglas locales y
  conservativas.

---

## 7. Paralelización

### 7.1 Descomposición 1D

- $P=$ `bsp_nprocs()`.
- `base = N/P`, `rem = N%P`.
- Rank `r` recibe:

```text
nx(r) = base + (r < rem ? 1 : 0)
x0(r) = r*base + min(r, rem)
```

- No se exige que `N` sea divisible por `P`.
- Ejemplo `N=10`, `P=3`: `nx = [4,3,3]`, `x0 = [0,4,7]`.
  - rank 0: `[0,4)`, vecinos 2 (izq.) y 1 (der.).
  - rank 1: `[4,7)`, vecinos 0 y 2.
  - rank 2: `[7,10)`, vecinos 1 y 0.

### 7.2 Malla local y halos

```text
x ∈ [-1, nx]     halo real (viene del vecino)
y ∈ [-1, N]      halo periódico local
z ∈ [-1, N]      halo periódico local

idx(x,y,z) = ((x+1)*(N+2) + (y+1))*(N+2) + (z+1)
```

- `LY = LZ = N+2`, `layer = LY*LZ`, `cells = (nx+2)*layer`.
- Buffers: `state`, `tmp` (colisión), `ghost_lo`, `ghost_hi`.

### 7.3 Intercambio de halos

Cada rank registra `ghost_lo` y `ghost_hi` con `bsp_push_reg`.

```cpp
bsp_put(left,  &tmp[idx(0, -1, -1)],      ghost_hi.data(), 0, bytes);
bsp_put(right, &tmp[idx(nx - 1, -1, -1)], ghost_lo.data(), 0, bytes);
bsp_sync();
memcpy(&tmp[idx(-1, -1, -1)], ghost_lo.data(), bytes);
memcpy(&tmp[idx(nx, -1, -1)], ghost_hi.data(), bytes);
```

- `ghost_lo` (x=-1) recibe la frontera alta del vecino izquierdo.
- `ghost_hi` (x=nx) recibe la frontera baja del vecino derecho.
- Para `P=1` se copia la frontera opuesta localmente (anillo de un proceso).

### 7.4 Orden del superpaso

1. `collide` (OpenMP) sobre celdas reales.
2. `fill_yz_ghosts` (copia local en y/z).
3. `bsp_put` de las dos caras + `bsp_sync`.
4. Copia de `ghost_lo`/`ghost_hi` a los halos de `tmp`.
5. `stream` (OpenMP) hacia `state`.

No se llama a BSPlib dentro de `omp parallel`: todas las llamadas desde el
hilo maestro.

### 7.5 OpenMP

- `omp_set_num_threads(T)` al inicio.
- `#pragma omp parallel for collapse(3) schedule(static)` en `collide`,
  `stream`, `init_state`, `local_stats`, densidad.
- `fill_yz_ghosts` paraleliza sobre x.
- `local_stats` usa `reduction(+:mass,px,py,pz)`.
- `sample_profile` paraleliza sobre x (no sobre y,z) para evitar carreras en
  el perfil local.

### 7.6 Salida y reducciones

- Densidad local con `popcount` → `bsp_send` a rank 0 → `.vti`.
- Masa/momento/tiempos → `Stats` por `bsp_send` a rank 0.
- Diagnóstico temporal → CSV en rank 0.
- Solo rank 0 escribe archivos.

### 7.7 Volumen de comunicación

Por paso, cada rank envía **2 caras** de $(N+2)^2$ bytes:

```text
V_rank  = 2 (N+2)^2 bytes
V_total = 2 P (N+2)^2 bytes
```

Con `N=256`, `P=4`:

- 66.6 KB por cara.
- 133 KB por rank/paso.
- 0.53 MB por paso total.
- 0.53 GB en 1000 pasos.

El volumen local por rank es $n_x N^2\approx4.2$ MB, así que la comunicación
es $\approx3\%$ del volumen por paso. El costo dominante es la barrera y el
alltoall de BSPonMPI, no el ancho de banda.

---

## 8. Mapa del código

### 8.1 Árbol

```text
.
├── CMakeLists.txt
├── README.md
├── AI_HANDOFF.md
├── src/gas.cpp
├── tests/hello.cpp
├── tests/halo.cpp
├── docker/Dockerfile
├── scripts/
│   ├── build.sh
│   ├── run.sh
│   ├── scaling.sh
│   ├── docker_shell.sh
│   ├── make_comm_figures.py
│   ├── make_figures.py
│   └── plot_evolution.py
├── docs/
│   ├── informe.tex / informe.pdf
│   ├── presentacion_clase.tex / .pdf   # versión para la clase
│   ├── presentacion.tex / presentacion.pdf
│   ├── JUSTIFICACION_GAS_IDEAL.md
│   └── figs/
├── data/          # VTK/CSV de corridas (ignorado)
├── results/       # CSV de scaling (ignorado)
├── ejemplos/      # material del curso
└── build/         # compilación (ignorado)
```

### 8.2 `src/gas.cpp`

| Sección | Qué hace |
|---|---|
| `Config` | `N`, `steps`, `threads`, `out_every`, `prefix`, gas y geometría. |
| `GasInfo`, `select_gas` | Neón/helio para el mapeo. |
| `load_env` | Lee todas las variables de entorno. |
| `build_collision_table` | Genera las 8 tablas conservativas. |
| `target_density` | Densidad objetivo según esfera/cámara y `RAMP`. |
| `sample_state` | Muestrea las 6 direcciones con `p=ρ/6` y hash global. |
| `initial_state` | Combina `target_density` + `sample_state`. |
| `Grid` | Malla local, `idx()` y `row()`. |
| `init_state` | Llena el estado inicial con OpenMP. |
| `collide` | Colisión local `state -> tmp`. |
| `fill_yz_ghosts` | Halos periódicos en y/z. |
| `exchange_x` | `bsp_put` + `bsp_sync` + copia de halos. |
| `stream` | Gather desde `tmp` hacia `state`. |
| `Stats`, `local_stats`, `gather_stats` | Masa/momento/tiempos. |
| `physical_dx_mm`, `physical_dt_us` | Escalas físicas para el CSV. |
| `sample_profile` | Perfil global $\rho(x,t)$ y métricas. |
| `write_vti` | Transpone a orden VTK y escribe `.vti` raw. |
| `gather_and_write` | Reúne densidades en rank 0 y escribe. |
| `bsp_main` | Driver: begin, registros, loop, salida, end. |
| `main`, `parse_args`, `usage` | CLI y `bsp_init`. |

### 8.3 CLI

```text
./gas [N] [pasos] [hilos] [cada_cuanto_vtk] [prefijo]

N            lado de la malla (default 256)
pasos        pasos de simulación (default 1000)
hilos        hilos OpenMP por proceso (default 1)
cada_cuanto  0 = VTK inicio/final; >0 = cada n pasos; <0 = sin VTK
prefijo      prefijo de los .vti (default gas)
```

Ejemplo:

```bash
GEOM=chamber P_RATIO=10 RAMP=8 \
GAS=neon DX_MM=1 T0_K=300 P0_ATM=1 \
OMP_NUM_THREADS=4 mpirun -np 4 --bind-to none \
  ./build/gas 256 1000 4 100 chamber
```

### 8.4 Variables de entorno

| Variable | Default | Significado |
|---|---|---|
| `GAS` | `neon` | `neon` o `helium` |
| `DX_MM` | `1` | tamaño de celda en mm; `0` = 1 partícula = 1 molécula |
| `T0_K` | `300` | temperatura de referencia |
| `P0_ATM` | `1` | presión de referencia |
| `GEOM` | `sphere` | `sphere` o `chamber` |
| `RHO_H` | `3.0` | densidad alta (partículas/celda) |
| `RHO_L` | `0` (`0.3` en chamber) | densidad baja/fondo |
| `P_RATIO` | — | si se define, `RHO_L = RHO_H / P_RATIO` |
| `IFACE_FRAC` | `0.25` | posición de la interfaz (fracción de `N`) |
| `RAMP` | `0` | ancho del gradiente en celdas (0 = escalón) |
| `PROFILE_CSV` | — | escribe perfil $\rho(x,t)$ y métricas por paso |
| `PROFILE_EVERY` | `1` | cada cuántos pasos muestrear |

### 8.5 Formatos de salida

**VTK `.vti`:**

- XML ImageData, `UInt8` “density”.
- Datos binarios *appended* raw, little-endian, prefijo de longitud `UInt32`.
- Orden VTK: `x + N*(y + N*z)`; el código transpone desde el orden de
  simulación `(x*N + y)*N + z`.
- Suma de densidad = masa global.

**CSV de diagnóstico (`PROFILE_CSV`):**

```text
# gas=... dx_mm=... dt_us=... N=... rho_h=... rho_l=... iface_frac=...
step,time_us,mass,mx,my,mz,max_rho,com_x_mm,rms_x_mm,front_x_mm,right_frac,rho_0,...,rho_{N-1}
```

`right_frac` es la fracción de masa con $x\ge$ `iface_frac*N` (útil para la
cámara; para la esfera es solo una referencia).

---

## 9. Validación

### 9.1 Criterios

1. Masa y momento exactos antes/después.
2. Halos: test sintético con 0 errores.
3. Serial vs. paralelo: `P=1` y `P=4` con `.vti` idénticos byte a byte.
4. VTK: suma de densidad = masa.
5. Híbrido: 4 procesos × 4 hilos.

### 9.2 Evidencia observada

- `tests/hello.cpp`: `pid 0..3: 4 hilos OpenMP`.
- `tests/halo.cpp`: `errores=0` en los 4 ranks.
- Esfera `N=64`, 50 pasos: `mass0=mass1=6589`; momento igual.
- Esfera `N=256`: `mass0=mass1=412365`; momento igual.
- Cámara `N=256`, 10:1: `mass0=mass1=16357802`; momento igual;
  `p_H=101325 Pa`, `p_L=10132.5 Pa`.
- `P=3`, `N=64` (desbalance): masa/momento conservados.
- `P=2`, `N=63` (impar): masa/momento conservados.
- `P=1` vs `P=4`, esfera `N=64`, 50 pasos: `cmp` de `.vti` → idénticos.
- `P=1` vs `P=4`, cámara `N=64`, 20 pasos, `RAMP=8`: idénticos.
- VTK: prefijo de longitud $N^3$, suma de densidad = masa, XML parsea.

### 9.3 Comandos de validación

```bash
# tests
./scripts/docker_shell.sh bash -lc \
  'mpirun -np 4 --bind-to none ./build/hello 4'
./scripts/docker_shell.sh bash -lc \
  'mpirun -np 4 --bind-to none ./build/halo_test 4'

# comparación P=1 vs P=4
mkdir -p data/par1 data/par4
docker run --rm --user "$(id -u):$(id -g)" -e HOME=/tmp \
  -v "$PWD":/work -w /work/data/par1 hpc-lab2:latest \
  /work/build/gas 64 50 1 0 g
docker run --rm --user "$(id -u):$(id -g)" -e HOME=/tmp \
  -v "$PWD":/work -w /work/data/par4 hpc-lab2:latest \
  bash -lc 'OMP_NUM_THREADS=4 mpirun -np 4 --bind-to none \
            /work/build/gas 64 50 4 0 g'
cmp data/par1/g_00050.vti data/par4/g_00050.vti
```

---

## 10. Rendimiento

Laptop 8C/16T, Docker, `N=256`, 200 pasos, sin VTK:

| P | T | hilos | tiempo (s) | speedup vs 1×1 |
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

Interpretación:

- El problema es **memory-bound**.
- En un solo nodo, el mejor tiempo es **OpenMP puro** (1×8); 1×16 no mejora
  porque la máquina tiene 8 núcleos físicos y SMT no aporta.
- La capa de procesos no mejora el tiempo en un solo nodo: el
  `bsp_sync`/alltoall de BSPonMPI compensa la división.
- La fracción de comunicación es ~7% en 4×4 y ~45% en 4×2 (poco trabajo por
  rank).
- El objetivo 4×8 requiere un nodo de 32 hilos.

**Ojo:** las corridas con `PROFILE_CSV` agregan un `bsp_sync` por paso; no
usarlas para medir rendimiento. La tabla se midió sin perfil.

---

## 11. Documentos y figuras

### 11.1 Entregables

- `docs/informe.pdf` (15 páginas, fuente `docs/informe.tex`):
  modelo, reglas, paralelización, validación, resultados, evolución temporal
  y apéndices. **Incluye resultados.**
- `docs/presentacion_clase.pdf` (15 páginas, fuente
  `docs/presentacion_clase.tex`): **versión que se presenta en clase**. Solo
  el sistema, las reglas del CA y la paralelización híbrida (4 slides:
  dos niveles, superpaso BSP, modelo de costo `w + h·g + L` aplicado y
  dimensionamiento N=1024 para 4×8). Slides copiadas desde
  `presentacion.tex`; los cambios para la clase se hacen solo aquí.
- `docs/presentacion.pdf` (29 páginas con divisores de sección, fuente
  `docs/presentacion.tex`):
  **propuesta** con problema, modelo, decisiones, plan de validación, plan de
  experimentos, riesgos y alternativas. **No incluye resultados.**
- `docs/JUSTIFICACION_GAS_IDEAL.md`: justificación física y mapeo de unidades.
- `src/gas.cpp`, `tests/`, `CMakeLists.txt`, `docker/`, `scripts/`.
- `README.md` y este documento.

Compilar:

```bash
cd docs
latexmk -pdf informe.tex
latexmk -pdf presentacion.tex
latexmk -pdf presentacion_clase.tex
```

### 11.2 Figuras de diseño (sin datos de corridas)

| Archivo | Qué muestra | Dónde se usa |
|---|---|---|
| `decomposition.png` | Reparto 1D y halos, `P=4` | Informe y presentación |
| `decomposition_uneven.png` | Caso `N=10`, `P=3`, resto | Informe y presentación |
| `communication_volume.png` | Qué se comunica (2D) | Informe y presentación |
| `communication_volume_3d.png` | Vista 3D del subdominio | Informe y presentación |
| `halo_exchange.png` | Detalle `bsp_put`/`bsp_sync` | Informe y presentación |
| `hybrid_threads.png` | Procesos × hilos | Informe y presentación |

### 11.3 Figuras de resultados

| Archivo | Qué muestra | Dónde |
|---|---|---|
| `chamber_slices.png` | Cortes de la cámara en 4 tiempos | Informe |
| `spacetime_chamber.png` | Diagrama espacio-temporal | Informe |
| `evolution_chamber.png` | Métricas globales | Informe |
| `evolution_chamber.gif` | Animación del perfil | Informe (texto) |
| `scaling.png` | Tiempo y speedup | Informe |
| `comm_fraction.png` | Fracción de comunicación | Informe |

### 11.4 Figuras auxiliares no versionadas

Estas figuras se pueden regenerar con los scripts, pero no se versionan para
no cargar el repositorio: `chamber_profile.png`,
`chamber_slices_wide.png`, `chamber_slices_wide2.png`, `sphere_slices.png`,
`spacetime_sphere.png`. Están incluidas en `.gitignore`.

### 11.5 Regenerar figuras

Los comandos que usan `/work/build/gas` se ejecutan **dentro del contenedor**
(`/work` es el repo). En nativo, cambiar `/work/build/gas` por la ruta
equivalente a `build/gas`.

```bash
# esquemas de diseño (no necesitan datos)
python3 scripts/make_comm_figures.py

# figuras de simulación (requieren data/figs)
# 1) cámara y esfera con frames
cd data/figs
GEOM=chamber P_RATIO=10 RAMP=8 OMP_NUM_THREADS=4 \
  mpirun -np 4 --bind-to none /work/build/gas 128 150 4 50 chamber2
GEOM=sphere OMP_NUM_THREADS=4 \
  mpirun -np 4 --bind-to none /work/build/gas 128 75 4 25 sphere2
# 2) perfiles temporales
GEOM=chamber P_RATIO=10 RAMP=8 \
  PROFILE_CSV=chamber_profile.csv PROFILE_EVERY=1 \
  OMP_NUM_THREADS=4 mpirun -np 4 --bind-to none \
  /work/build/gas 128 150 4 -1
GEOM=sphere PROFILE_CSV=sphere_profile.csv PROFILE_EVERY=1 \
  OMP_NUM_THREADS=4 mpirun -np 4 --bind-to none \
  /work/build/gas 128 100 4 -1
# 3) figuras
cd ../..
python3 scripts/make_figures.py
python3 scripts/plot_evolution.py
```

---

## 12. Riesgos, limitaciones y pendientes

### 12.1 Riesgos técnicos

1. **BSPonMPI 0.3** es antiguo y su build original usa SCons/Python 2.
   Mitigación: compilación manual y Docker reproducible.
2. **`bsp_get_tag` incluye el tag**: restar `sizeof(int)`.
3. **BSPlib no thread-safe**: llamadas solo desde el hilo maestro.
4. **Costo de comunicación** con `P` grande: medir; elegir `P` según el tamaño
   local.
5. **No hay nodo de 32 hilos** en la laptop: correr 4×4 local; 4×8 en cluster.
6. **Modelo isotérmico**: no captura enfriamiento adiabático.
7. **Anisotropía del LGA**: comparar tendencias cualitativas; evaluar partícula
   en reposo.
8. **Bordes periódicos**: la cámara es un toro; paredes reflectantes quedan
   como opción.

### 12.2 Limitaciones físicas

- LGA 6-velocidades: anisotropía cúbica, poca termalización, ondas
  artificiales.
- Gas ideal: sin virial; a 1 atm la corrección es 0.05%, aceptable.
- Sin ecuación de energía: temperatura constante.
- Coeficientes de transporte no son los del neón real.

### 12.3 Pendientes priorizados

1. Correr 4×8 en un nodo de 32 hilos y completar strong/weak scaling.
2. Paredes reflectantes opcionales.
3. LGA térmico (reposo/energía) para expansión adiabática.
4. Exportar presión/velocidad al VTK.
5. Optimizar streaming (`-march=native`, sin ramas) y evaluar 2D/3D
   decomposition.
6. Tests CTest con aserciones de invariantes y comparación P=1 vs P=4.
7. Revisión final de PDF/Beamer (contenido, estilo, autoría).

---

## 13. Decisiones de diseño (y por qué)

| Decisión | Motivo |
|---|---|
| LGA 6 bits, sin reposo | CA estricto, entero, conservativo, 1 byte/celda; sin condensación por reposo. |
| 8 tablas + hash global | Rompe simetrías sin perder invariantes; independiente de `P`. |
| Descomposición 1D en x | Solo 2 vecinos, simple; suficiente para 4–8 ranks. |
| Buffers de halo separados | Offsets independientes de `nx` cuando `N` no divide a `P`. |
| `bsp_put` para halos | Evita cola BSMP y `qsize`/`move` por paso. |
| OpenMP por fase | BSPlib no thread-safe; barrera implícita por fase. |
| VTK `.vti` appended raw | Un archivo por frame, binario compacto, sin dependencias. |
| Bordes periódicos | Scaling y simplicidad; paredes como futuro. |
| Esfera + cámara | Esfera baseline; cámara con gradiente explícito. |
| Neón por defecto | Monoatómico, `c_s` moderada, densidad fácil de verificar. |
| Presentación sin resultados | Es una propuesta; los resultados van en el informe. |

### 13.1 Alternativas consideradas

| Alternativa | Por qué no |
|---|---|
| LBM D3Q7 | Más realista pero no es CA estricto; flotantes. |
| Stencil PDE | Simple pero no es un CA. |
| MPI + OpenMP directo | No cumple el requisito BSP/BSPlib. |
| Descomposición 2D/3D | Más compleja; 1D basta para el objetivo. |
| Paredes reflectantes | Se priorizó scaling; quedan como opción. |
| Solo esfera | No da un gradiente de presión explícito. |
| Helio | Alternativa válida; `c_s` más alta y menor densidad. |

---

## 14. Reproducibilidad desde cero

```bash
# 1) imagen
docker build -t hpc-lab2:latest -f docker/Dockerfile .

# 2) compilar
./scripts/docker_shell.sh bash -lc './scripts/build.sh'

# 3) tests
./scripts/docker_shell.sh bash -lc \
  'mpirun -np 4 --bind-to none ./build/hello 4'
./scripts/docker_shell.sh bash -lc \
  'mpirun -np 4 --bind-to none ./build/halo_test 4'

# 4) corrida de humo
./scripts/docker_shell.sh bash -lc \
  'GEOM=chamber P_RATIO=10 RAMP=8 OMP_NUM_THREADS=4 \
   mpirun -np 4 --bind-to none ./build/gas 128 100 4 -1'

# 5) scaling (sin VTK)
./scripts/docker_shell.sh bash -lc \
  'N=256 STEPS=200 ./scripts/scaling.sh'

# 6) figuras de diseño
python3 scripts/make_comm_figures.py

# 7) figuras de simulación
#    (correr los comandos de §11.5 para generar data/figs)
python3 scripts/make_figures.py
python3 scripts/plot_evolution.py

# 8) documentos
cd docs
latexmk -pdf informe.tex
latexmk -pdf presentacion.tex
latexmk -pdf presentacion_clase.tex
```

---

## 15. Cómo continuar en una sesión nueva

1. Leer este archivo y `README.md`.
2. Revisar `git status` y no borrar `ejemplos/`.
3. Construir la imagen y correr los tests de §14.
4. Si se cambia el modelo:
   - mantener conservación de masa/momento;
   - volver a comparar `P=1` vs `P=4` byte a byte;
   - actualizar las figuras y tablas.
5. Si se cambia la comunicación:
   - respetar el quirk de `bsp_get_tag` (§5.4);
   - no llamar BSPlib desde hilos;
   - verificar que los halos se lean después del `bsp_sync`.
6. Si se cambia la presentación:
   - la que se expone en clase es `docs/presentacion_clase.tex`; editar
     **solo ese archivo** para la clase (sus slides son copias de
     `presentacion.tex` y no se sincronizan solas);
   - `presentacion.tex` es la propuesta completa: recordar que es una
     **propuesta sin resultados** y mantener las diapositivas de
     justificación, riesgos y alternativas.
7. Actualizar `README.md`, este archivo y los PDFs.

---

## 16. Referencias

- L. G. Valiant, *A Bridging Model for Parallel Computation*, CACM 33(8),
  1990.
- J. M. D. Hill et al., *BSPlib: The BSP Programming Library*, Parallel
  Computing 24, 1998.
- W. J. Suijlen, *BSPonMPI: BSPlib on top of MPI*,
  <https://sourceforge.net/projects/bsponmpi/>.
- U. Frisch, B. Hasslacher, Y. Pomeau, *Lattice-Gas Automata for the
  Navier–Stokes Equation*, PRL 56, 1986.
- J. H. Dymond y E. B. Smith, *The Virial Coefficients of Pure Gases and
  Mixtures*, Oxford University Press.
- OpenMP Architecture Review Board, <https://www.openmp.org/>.
- Documentación del curso en `ejemplos/`.
