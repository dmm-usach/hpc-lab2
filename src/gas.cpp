// ---------------------------------------------------------------------------
// hpc-lab2: expansion de un gas 3D con un automata celular de gas reticular
// (LGA, 6 direcciones, sin reposo) sobre BSPlib (BSPonMPI) + OpenMP.
//
// Modelo:
//   - Malla cubica N x N x N, estado por celda = mascara de 6 bits
//     (una particula por direccion, principio de exclusion).
//   - Paso: colision local (tabla que conserva masa y momento) y luego
//     streaming hacia los 6 vecinos.
//   - Condicion inicial: esfera densa (p=1/2 por direccion) en el centro.
//   - Bordes periodicos en x (anillo entre ranks) y en y/z (copia local).
//
// Paralelizacion:
//   - Procesos: descomposicion 1D en x. Cada rank tiene nx capas reales
//     y una capa fantasma a cada lado. Intercambio con bsp_push_reg +
//     bsp_put + bsp_sync (superpaso BSP).
//   - Hilos: OpenMP sobre el interior (collide y stream), sin tocar halos.
//
// Compilar: ver CMakeLists.txt (requiere -DOMPI_SKIP_MPICXX con OpenMPI 4).
// Ejecutar:  OMP_NUM_THREADS=4 mpirun -np 4 --bind-to none ./gas 256 1000 4 0 gas
// ---------------------------------------------------------------------------

extern "C" {
#include <bsp.h>
}

#include <omp.h>

#include <algorithm>
#include <array>
#include <cstddef>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <map>
#include <random>
#include <string>
#include <vector>

namespace {

// ---------------------------------------------------------------------------
// Configuracion
// ---------------------------------------------------------------------------
struct Config {
    int N = 256;          // lado de la malla cubica global
    int steps = 1000;     // pasos de simulacion
    int threads = 1;      // hilos OpenMP por proceso
    int out_every = 0;    // 0: solo inicio/final; >0: cada n pasos; <0: sin VTK
    std::string prefix = "gas";

    // Gas concreto para el mapeo lattice -> SI
    std::string gas = "neon";
    double dx_mm = 1.0;   // tamano de celda en mm; 0 = 1 particula = 1 molecula
    double T0 = 300.0;    // K
    double p0_atm = 1.0;  // atm

    // Geometria inicial
    std::string geom = "sphere";  // sphere | chamber
    double rho_h = 3.0;           // densidad alta (particulas/celda)
    double rho_l = 0.0;           // densidad baja / fondo
    double iface_frac = 0.25;     // interfaz de la camara (fraccion de N)
    double ramp = 0.0;            // ancho del gradiente en celdas (0 = escalon)

    // Diagnostico temporal: perfil x(t) y metricas por paso
    std::string profile_csv;      // si no esta vacio, escribe el CSV
    int profile_every = 1;        // cada cuantos pasos muestrear
};

Config g_cfg;

// Gas concreto solo para el mapeo de unidades lattice -> SI.
// La dinamica del LGA no depende de la especie: m y gamma entran en la
// conversion de escalas y en las magnitudes fisicas reportadas.
struct GasInfo {
    const char *name;
    double m_u;    // masa molecular en unidades de masa atomica
    double gamma;  // cp/cv (gas monoatomico = 5/3)
};

GasInfo select_gas(const char *requested)
{
    static const GasInfo gases[] = {
        {"neon", 20.1797, 5.0 / 3.0},
        {"helium", 4.002602, 5.0 / 3.0},
    };
    if (requested != nullptr) {
        for (const GasInfo &g : gases) {
            if (std::strcmp(requested, g.name) == 0) {
                return g;
            }
        }
    }
    return gases[0];  // neon por defecto
}

// Lee la configuracion fisica/geometrica desde el entorno. Se llama antes de
// bsp_init, de modo que todos los ranks tengan la misma configuracion.
void load_env(Config &c)
{
    if (const char *v = std::getenv("GAS")) c.gas = v;
    if (const char *v = std::getenv("DX_MM")) c.dx_mm = std::atof(v);
    if (const char *v = std::getenv("T0_K")) c.T0 = std::atof(v);
    if (const char *v = std::getenv("P0_ATM")) c.p0_atm = std::atof(v);

    if (const char *v = std::getenv("GEOM")) c.geom = v;
    if (const char *v = std::getenv("RHO_H")) c.rho_h = std::atof(v);

    bool rho_l_set = false;
    if (const char *v = std::getenv("RHO_L")) {
        c.rho_l = std::atof(v);
        rho_l_set = true;
    }
    if (const char *v = std::getenv("P_RATIO")) {
        const double r = std::atof(v);
        if (r > 0.0) {
            c.rho_l = c.rho_h / r;
            rho_l_set = true;
        }
    }
    if (!rho_l_set) {
        c.rho_l = (c.geom == "chamber") ? 0.3 : 0.0;
    }

    if (const char *v = std::getenv("IFACE_FRAC")) c.iface_frac = std::atof(v);
    if (const char *v = std::getenv("RAMP")) c.ramp = std::atof(v);
    if (const char *v = std::getenv("PROFILE_CSV")) c.profile_csv = v;
    if (const char *v = std::getenv("PROFILE_EVERY")) c.profile_every = std::atoi(v);
    if (c.profile_every < 1) c.profile_every = 1;

    c.rho_h = std::clamp(c.rho_h, 0.0, 6.0);
    c.rho_l = std::clamp(c.rho_l, 0.0, 6.0);
    c.iface_frac = std::clamp(c.iface_frac, 0.0, 1.0);
    c.ramp = std::max(0.0, c.ramp);
}

// ---------------------------------------------------------------------------
// Automata: 6 direcciones (sin particula en reposo)
// ---------------------------------------------------------------------------
constexpr int NDIR = 6;
constexpr int NSTATE = 1 << NDIR;  // 64 estados posibles
constexpr int NTAB = 8;            // tablas de colision para variar por celda
constexpr int TAG_BYTES = static_cast<int>(sizeof(int));  // tamano de tag BSMP

constexpr uint8_t B_XP = 1u << 0;
constexpr uint8_t B_XM = 1u << 1;
constexpr uint8_t B_YP = 1u << 2;
constexpr uint8_t B_YM = 1u << 3;
constexpr uint8_t B_ZP = 1u << 4;
constexpr uint8_t B_ZM = 1u << 5;

const int DX[NDIR] = {+1, -1, 0, 0, 0, 0};
const int DY[NDIR] = {0, 0, +1, -1, 0, 0};
const int DZ[NDIR] = {0, 0, 0, 0, +1, -1};

uint8_t g_table[NTAB][NSTATE];

inline uint32_t mix32(uint32_t x)
{
    x ^= x >> 16;
    x *= 0x7feb352du;
    x ^= x >> 15;
    x *= 0x846ca68bu;
    x ^= x >> 16;
    return x;
}

inline uint32_t coord_hash(int x, int y, int z, int step)
{
    uint32_t h = static_cast<uint32_t>(x) * 0x9e3779b1u
               ^ static_cast<uint32_t>(y) * 0x85ebca77u
               ^ static_cast<uint32_t>(z) * 0xc2b2ae3du
               ^ static_cast<uint32_t>(step) * 0x27d4eb2fu;
    return mix32(h);
}

// Construye tablas de colision: para cada estado, se elige al azar otro
// estado con la misma masa y el mismo momento. Como la eleccion es dentro
// del mismo grupo de invariantes, la colision conserva masa y momento.
void build_collision_table()
{
    std::map<std::array<int, 4>, std::vector<int>> groups;

    for (int s = 0; s < NSTATE; ++s) {
        int mass = 0, px = 0, py = 0, pz = 0;
        for (int d = 0; d < NDIR; ++d) {
            if (s & (1 << d)) {
                ++mass;
                px += DX[d];
                py += DY[d];
                pz += DZ[d];
            }
        }
        groups[{mass, px, py, pz}].push_back(s);
    }

    std::mt19937 rng(0x5eed1234u);
    for (const auto &kv : groups) {
        const std::vector<int> &g = kv.second;
        for (int k = 0; k < NTAB; ++k) {
            for (int s : g) {
                g_table[k][s] = static_cast<uint8_t>(g[rng() % g.size()]);
            }
        }
    }
}

// Densidad objetivo (particulas/celda) segun la geometria inicial.
//   sphere:  rho_h dentro de la esfera central, rho_l fuera.
//   chamber: rho_h en x < interfaz, rho_l en x >= interfaz, con un ramp
//            lineal opcional de `ramp` celdas centrado en la interfaz.
double target_density(int gx, int y, int z)
{
    const int N = g_cfg.N;

    if (g_cfg.geom == "chamber") {
        const double xi = g_cfg.iface_frac * N;
        const double w = g_cfg.ramp;
        double t;
        if (w <= 0.0) {
            t = (gx < xi) ? 0.0 : 1.0;
        } else {
            t = (gx - (xi - 0.5 * w)) / w;
            t = std::clamp(t, 0.0, 1.0);
        }
        return g_cfg.rho_h + (g_cfg.rho_l - g_cfg.rho_h) * t;
    }

    const double c = 0.5 * (N - 1);
    const double r = std::max(6.0, 0.125 * N);
    const double dx = gx - c, dy = y - c, dz = z - c;
    return (dx * dx + dy * dy + dz * dz <= r * r) ? g_cfg.rho_h : g_cfg.rho_l;
}

// Muestrea las 6 direcciones con probabilidad rho/6 por direccion. El hash usa
// coordenadas globales, asi que el estado inicial no depende de P.
inline uint8_t sample_state(double rho, int gx, int y, int z)
{
    rho = std::clamp(rho, 0.0, 6.0);
    if (rho <= 0.0) return 0;
    if (rho >= 6.0) return 0x3f;

    const double p = rho / 6.0;
    uint8_t v = 0;
    for (int d = 0; d < NDIR; ++d) {
        const uint32_t seed =
            0x51ed2701u + static_cast<uint32_t>(d) * 0x9e3779b9u;
        const uint32_t h = coord_hash(gx, y, z, static_cast<int>(seed));
        const double u = (h & 0x00ffffffu) / 16777216.0;
        if (u < p) v |= static_cast<uint8_t>(1u << d);
    }
    return v;
}

inline uint8_t initial_state(int gx, int y, int z)
{
    return sample_state(target_density(gx, y, z), gx, y, z);
}

// ---------------------------------------------------------------------------
// Malla local con halos
//   x en [-1, nx]   (halo real, viene del rank vecino)
//   y, z en [-1, N] (halo periodico local)
//   idx(x,y,z) = ((x+1)*LY + (y+1))*LZ + (z+1)
// ---------------------------------------------------------------------------
struct Grid {
    int nx = 0;         // capas reales en x
    int x0 = 0;         // primer indice global en x
    int N = 0;          // lado global (y, z)
    int LY = 0, LZ = 0;
    size_t layer = 0;   // LY * LZ (una capa en x)
    size_t cells = 0;   // (nx+2) * layer

    void init(int nx_, int x0_, int N_)
    {
        nx = nx_;
        x0 = x0_;
        N = N_;
        LY = N + 2;
        LZ = N + 2;
        layer = static_cast<size_t>(LY) * LZ;
        cells = static_cast<size_t>(nx + 2) * layer;
    }

    inline size_t idx(int x, int y, int z) const
    {
        return ((static_cast<size_t>(x + 1) * LY + (y + 1)) * LZ) + (z + 1);
    }

    // indice de (x,y,z=0) para recorrer z de forma contigua
    inline size_t row(int x, int y) const
    {
        return ((static_cast<size_t>(x + 1) * LY + (y + 1)) * LZ) + 1;
    }
};

// ---------------------------------------------------------------------------
// Inicializacion
// ---------------------------------------------------------------------------
void init_state(const Grid &G, std::vector<uint8_t> &state)
{
#pragma omp parallel for collapse(3) schedule(static)
    for (int x = 0; x < G.nx; ++x) {
        for (int y = 0; y < G.N; ++y) {
            for (int z = 0; z < G.N; ++z) {
                state[G.idx(x, y, z)] = initial_state(G.x0 + x, y, z);
            }
        }
    }
}

// ---------------------------------------------------------------------------
// Colision local: state -> tmp (solo celdas reales)
// ---------------------------------------------------------------------------
void collide(const Grid &G, const std::vector<uint8_t> &state,
             std::vector<uint8_t> &tmp, int step)
{
    const int nx = G.nx, N = G.N, x0 = G.x0;

#pragma omp parallel for collapse(3) schedule(static)
    for (int x = 0; x < nx; ++x) {
        for (int y = 0; y < N; ++y) {
            for (int z = 0; z < N; ++z) {
                const size_t i = G.idx(x, y, z);
                const uint32_t h = coord_hash(x0 + x, y, z, step);
                tmp[i] = g_table[h & (NTAB - 1)][state[i]];
            }
        }
    }
}

// ---------------------------------------------------------------------------
// Halos periodicos en y/z dentro de tmp (copia local, sin comunicacion)
// ---------------------------------------------------------------------------
void fill_yz_ghosts(const Grid &G, std::vector<uint8_t> &tmp)
{
    const int N = G.N;

#pragma omp parallel for schedule(static)
    for (int x = 0; x < G.nx; ++x) {
        // z = -1 y z = N
        for (int y = 0; y < N; ++y) {
            tmp[G.idx(x, y, -1)] = tmp[G.idx(x, y, N - 1)];
            tmp[G.idx(x, y, N)] = tmp[G.idx(x, y, 0)];
        }
        // y = -1 y y = N (incluye esquinas en z, ya copiadas)
        for (int z = -1; z <= N; ++z) {
            tmp[G.idx(x, -1, z)] = tmp[G.idx(x, N - 1, z)];
            tmp[G.idx(x, N, z)] = tmp[G.idx(x, 0, z)];
        }
    }
}

// ---------------------------------------------------------------------------
// Intercambio de halos en x: superpaso BSP
//   ghost_lo (x=-1) recibe la frontera alta del vecino izquierdo
//   ghost_hi (x=nx) recibe la frontera baja del vecino derecho
// ---------------------------------------------------------------------------
void exchange_x(const Grid &G, std::vector<uint8_t> &tmp,
                std::vector<uint8_t> &ghost_lo,
                std::vector<uint8_t> &ghost_hi,
                int left, int right)
{
    const size_t bytes = G.layer;
    const uint8_t *lo = &tmp[G.idx(0, -1, -1)];
    const uint8_t *hi = &tmp[G.idx(G.nx - 1, -1, -1)];

    if (bsp_nprocs() == 1) {
        // anillo de un solo proceso: el halo es la propia frontera opuesta
        std::memcpy(ghost_lo.data(), hi, bytes);
        std::memcpy(ghost_hi.data(), lo, bytes);
    } else {
        bsp_put(left, lo, ghost_hi.data(), 0, bytes);
        bsp_put(right, hi, ghost_lo.data(), 0, bytes);
    }

    bsp_sync();

    std::memcpy(&tmp[G.idx(-1, -1, -1)], ghost_lo.data(), bytes);
    std::memcpy(&tmp[G.idx(G.nx, -1, -1)], ghost_hi.data(), bytes);
}

// ---------------------------------------------------------------------------
// Streaming: tmp (post-colision) -> state (t+1), solo celdas reales
// ---------------------------------------------------------------------------
void stream(const Grid &G, const std::vector<uint8_t> &tmp,
            std::vector<uint8_t> &state)
{
    const int nx = G.nx, N = G.N;
    const size_t L = G.layer, W = static_cast<size_t>(G.LZ);

#pragma omp parallel for collapse(3) schedule(static)
    for (int x = 0; x < nx; ++x) {
        for (int y = 0; y < N; ++y) {
            for (int z = 0; z < N; ++z) {
                const size_t r = G.row(x, y);
                const size_t xm = r - L, xp = r + L;
                const size_t ym = r - W, yp = r + W;

                uint8_t v = 0;
                if (tmp[xm + z] & B_XP) v |= B_XP;  // venia de x-1 con +x
                if (tmp[xp + z] & B_XM) v |= B_XM;  // venia de x+1 con -x
                if (tmp[ym + z] & B_YP) v |= B_YP;
                if (tmp[yp + z] & B_YM) v |= B_YM;
                if (tmp[r + z - 1] & B_ZP) v |= B_ZP;
                if (tmp[r + z + 1] & B_ZM) v |= B_ZM;

                state[r + z] = v;
            }
        }
    }
}

// ---------------------------------------------------------------------------
// Metricas: masa y momento total (invariantes del modelo)
// ---------------------------------------------------------------------------
struct Stats {
    unsigned long long mass = 0;
    long long px = 0, py = 0, pz = 0;
    double t_comp = 0.0;
    double t_comm = 0.0;
};

Stats local_stats(const Grid &G, const std::vector<uint8_t> &state)
{
    unsigned long long mass = 0;
    long long px = 0, py = 0, pz = 0;

#pragma omp parallel for collapse(3) schedule(static) \
    reduction(+ : mass, px, py, pz)
    for (int x = 0; x < G.nx; ++x) {
        for (int y = 0; y < G.N; ++y) {
            for (int z = 0; z < G.N; ++z) {
                const uint8_t v = state[G.idx(x, y, z)];
                mass += static_cast<unsigned long long>(
                    __builtin_popcount(static_cast<unsigned>(v)));
                px += (v & B_XP) ? 1 : 0;
                px -= (v & B_XM) ? 1 : 0;
                py += (v & B_YP) ? 1 : 0;
                py -= (v & B_YM) ? 1 : 0;
                pz += (v & B_ZP) ? 1 : 0;
                pz -= (v & B_ZM) ? 1 : 0;
            }
        }
    }

    Stats s;
    s.mass = mass;
    s.px = px;
    s.py = py;
    s.pz = pz;
    return s;
}

// Reduce las metricas de todos los ranks en el rank 0.
Stats gather_stats(const Grid &G, const std::vector<uint8_t> &state,
                   double t_comp, double t_comm)
{
    const int pid = bsp_pid();
    Stats mine = local_stats(G, state);
    mine.t_comp = t_comp;
    mine.t_comm = t_comm;

    Stats total;
    if (pid == 0) {
        total = mine;
        bsp_sync();
        int packets = 0;
        size_t nbytes = 0;
        bsp_qsize(&packets, &nbytes);
        for (int i = 0; i < packets; ++i) {
            int size = -1, tag = -1;
            bsp_get_tag(&size, &tag);
            if (size <= 0) {
                break;
            }
            // BSPonMPI informa payload + tag en el status de bsp_get_tag.
            const int payload = size - TAG_BYTES;
            if (payload != static_cast<int>(sizeof(Stats))) {
                std::fprintf(stderr, "gather_stats: payload inesperado %d\n",
                             payload);
                continue;
            }
            Stats other{};
            bsp_move(&other, static_cast<size_t>(payload));
            total.mass += other.mass;
            total.px += other.px;
            total.py += other.py;
            total.pz += other.pz;
            total.t_comp = std::max(total.t_comp, other.t_comp);
            total.t_comm = std::max(total.t_comm, other.t_comm);
        }
    } else {
        int tag = pid;
        bsp_send(0, &tag, &mine, sizeof(Stats));
        bsp_sync();
    }
    return total;  // valido solo en pid 0
}

// ---------------------------------------------------------------------------
// Salida VTK (.vti) con la densidad (0..6) por celda
// ---------------------------------------------------------------------------
void write_vti(const Config &cfg, const std::vector<uint8_t> &sim, int step)
{
    const int N = cfg.N;
    const size_t total = static_cast<size_t>(N) * N * N;
    if (total > 0xffffffffull) {
        std::fprintf(stderr, "write_vti: malla demasiado grande para UInt32\n");
        return;
    }

    // La simulacion guarda (x,y,z) con z contiguo; VTK espera x contiguo.
    std::vector<uint8_t> vtk(total);
#pragma omp parallel for collapse(2) schedule(static)
    for (int z = 0; z < N; ++z) {
        for (int y = 0; y < N; ++y) {
            for (int x = 0; x < N; ++x) {
                vtk[static_cast<size_t>(x)
                    + static_cast<size_t>(N) * (y + static_cast<size_t>(N) * z)] =
                    sim[(static_cast<size_t>(x) * N + y) * N + z];
            }
        }
    }

    char name[1024];
    std::snprintf(name, sizeof(name), "%s_%05d.vti", cfg.prefix.c_str(), step);

    FILE *f = std::fopen(name, "wb");
    if (!f) {
        std::perror(name);
        return;
    }

    std::fprintf(f,
                 "<?xml version=\"1.0\"?>\n"
                 "<VTKFile type=\"ImageData\" version=\"0.1\" "
                 "byte_order=\"LittleEndian\" header_type=\"UInt32\">\n"
                 "  <ImageData WholeExtent=\"0 %d 0 %d 0 %d\" "
                 "Origin=\"0 0 0\" Spacing=\"1 1 1\">\n"
                 "    <Piece Extent=\"0 %d 0 %d 0 %d\">\n"
                 "      <PointData Scalars=\"density\">\n"
                 "        <DataArray type=\"UInt8\" Name=\"density\" "
                 "format=\"appended\" offset=\"0\"/>\n"
                 "      </PointData>\n"
                 "    </Piece>\n"
                 "  </ImageData>\n"
                 "  <AppendedData encoding=\"raw\">\n"
                 "   _",
                 N - 1, N - 1, N - 1, N - 1, N - 1, N - 1);

    const uint32_t nbytes = static_cast<uint32_t>(vtk.size());
    std::fwrite(&nbytes, sizeof(nbytes), 1, f);
    std::fwrite(vtk.data(), 1, vtk.size(), f);
    std::fprintf(f, "\n  </AppendedData>\n</VTKFile>\n");
    std::fclose(f);
}

void gather_and_write(const Grid &G, const std::vector<uint8_t> &state,
                      int step, const Config &cfg)
{
    const int P = bsp_nprocs();
    const int pid = bsp_pid();

    std::vector<uint8_t> local(static_cast<size_t>(G.nx) * G.N * G.N);
#pragma omp parallel for collapse(3) schedule(static)
    for (int x = 0; x < G.nx; ++x) {
        for (int y = 0; y < G.N; ++y) {
            for (int z = 0; z < G.N; ++z) {
                local[(static_cast<size_t>(x) * G.N + y) * G.N + z] =
                    static_cast<uint8_t>(__builtin_popcount(
                        static_cast<unsigned>(state[G.idx(x, y, z)])));
            }
        }
    }

    std::vector<uint8_t> full;
    if (pid == 0) {
        full.resize(static_cast<size_t>(cfg.N) * cfg.N * cfg.N);
        std::memcpy(&full[static_cast<size_t>(G.x0) * cfg.N * cfg.N],
                    local.data(), local.size());
    } else {
        int tag = pid;
        bsp_send(0, &tag, local.data(), local.size());
    }

    bsp_sync();

    if (pid == 0) {
        int packets = 0;
        size_t nbytes = 0;
        bsp_qsize(&packets, &nbytes);
        const int base = cfg.N / P;
        const int rem = cfg.N % P;

        for (int i = 0; i < packets; ++i) {
            int size = -1, tag = -1;
            bsp_get_tag(&size, &tag);
            if (size <= 0) {
                break;
            }
            // BSPonMPI informa payload + tag en el status de bsp_get_tag.
            const int payload = size - TAG_BYTES;
            const int nx = base + (tag < rem ? 1 : 0);
            const int x0 = tag * base + std::min(tag, rem);
            const size_t expect = static_cast<size_t>(nx) * cfg.N * cfg.N;
            if (static_cast<size_t>(payload) != expect) {
                std::fprintf(stderr,
                             "gather: tamano inesperado de rank %d: %d vs %zu\n",
                             tag, payload, expect);
                continue;
            }
            bsp_move(&full[static_cast<size_t>(x0) * cfg.N * cfg.N],
                     static_cast<size_t>(payload));
        }

        write_vti(cfg, full, step);
    }

    bsp_sync();  // cierra el superpaso de salida
}

// ---------------------------------------------------------------------------
// Mapeo fisico: unidades lattice -> SI para un gas concreto.
//   rho_lat = 3 (cada direccion p=1/2) se toma como p0, T0.
//   Si DX_MM > 0: celda = DX_MM milimetros (mapeo macroscopico, cada
//                 particula lattice representa un conglomerado de moleculas).
//   Si DX_MM = 0: se elige dx para que 1 particula lattice = 1 molecula.
// ---------------------------------------------------------------------------
void print_physical_mapping(const GasInfo &gas, unsigned long long n_particles)
{
    const double kB = 1.380649e-23;      // J/K
    const double u = 1.66053906660e-27;  // kg/u
    const double atm = 101325.0;         // Pa

    const double T0 = g_cfg.T0;
    const double p0 = g_cfg.p0_atm * atm;
    const double dx_mm = g_cfg.dx_mm;

    const double m = gas.m_u * u;
    const double n0 = p0 / (kB * T0);           // moleculas/m^3
    const double dx = (dx_mm > 0.0)
                          ? dx_mm * 1e-3
                          : std::cbrt(3.0 * kB * T0 / p0);

    const double cell_volume = dx * dx * dx;
    const double molecules_per_particle = (n0 * cell_volume) / 3.0;
    const double mass_per_particle = molecules_per_particle * m;
    const double rho_mass = n0 * m;
    const double cs_iso = std::sqrt(kB * T0 / m);
    const double cs_ad = std::sqrt(gas.gamma) * cs_iso;
    const double dt = (1.0 / std::sqrt(3.0)) * dx / cs_iso;
    const double v_lat = dx / dt;
    const double domain = g_cfg.N * dx;
    const double mass_total = static_cast<double>(n_particles) * mass_per_particle;

    const double p_hi = (g_cfg.rho_h / 3.0) * p0;
    const double p_lo = (g_cfg.rho_l / 3.0) * p0;
    const double iface_m = g_cfg.iface_frac * g_cfg.N * dx;

    std::printf("gas=%s m=%.6g kg gamma=%.6g\n", gas.name, m, gas.gamma);
    std::printf("T0=%.6g K p0=%.6g Pa dx=%.6g m\n", T0, p0, dx);
    std::printf("cs_iso=%.6g m/s cs_ad=%.6g m/s dt=%.6g s v_lat=%.6g m/s\n",
                cs_iso, cs_ad, dt, v_lat);
    std::printf("rho_mass=%.6g kg/m3 molecules_per_particle=%.6g\n",
                rho_mass, molecules_per_particle);
    std::printf("geom=%s rho_H=%.6g rho_L=%.6g\n", g_cfg.geom.c_str(),
                g_cfg.rho_h, g_cfg.rho_l);
    std::printf("p_H=%.6g Pa p_L=%.6g Pa ratio=", p_hi, p_lo);
    if (g_cfg.rho_l > 0.0) {
        std::printf("%.6g\n", g_cfg.rho_h / g_cfg.rho_l);
    } else {
        std::printf("inf\n");
    }
    if (g_cfg.geom == "chamber") {
        std::printf("interface=%.6g m ramp=%.6g m\n", iface_m,
                    g_cfg.ramp * dx);
    }
    std::printf("domain=%.6g m mass_total=%.6g kg\n", domain, mass_total);
    std::fflush(stdout);
}

// ---------------------------------------------------------------------------
// Escalas fisicas para el CSV de evolucion
// ---------------------------------------------------------------------------
double physical_dx_mm()
{
    if (g_cfg.dx_mm > 0.0) return g_cfg.dx_mm;
    const double kB = 1.380649e-23;
    const double atm = 101325.0;
    const GasInfo gas = select_gas(g_cfg.gas.c_str());
    (void)gas;
    const double dx = std::cbrt(3.0 * kB * g_cfg.T0 / (g_cfg.p0_atm * atm));
    return dx * 1e3;
}

double physical_dt_us()
{
    const double kB = 1.380649e-23;
    const double u = 1.66053906660e-27;
    const GasInfo gas = select_gas(g_cfg.gas.c_str());
    const double m = gas.m_u * u;
    const double cs = std::sqrt(kB * g_cfg.T0 / m);
    const double dx = physical_dx_mm() * 1e-3;
    return (1.0 / std::sqrt(3.0)) * dx / cs * 1e6;
}

// ---------------------------------------------------------------------------
// Diagnostico temporal: perfil global de densidad en x y metricas.
// Todos los ranks llaman; solo rank 0 escribe el CSV.
// ---------------------------------------------------------------------------
void sample_profile(const Grid &G, const std::vector<uint8_t> &state, int step,
                    FILE *fp)
{
    const int P = bsp_nprocs();
    const int pid = bsp_pid();

    std::vector<double> prof(G.nx, 0.0);
    double mass = 0.0, px = 0.0, py = 0.0, pz = 0.0, max_rho = 0.0;

#pragma omp parallel for schedule(static) \
    reduction(+ : mass, px, py, pz) reduction(max : max_rho)
    for (int x = 0; x < G.nx; ++x) {
        for (int y = 0; y < G.N; ++y) {
            for (int z = 0; z < G.N; ++z) {
                const uint8_t v = state[G.idx(x, y, z)];
                const int c = __builtin_popcount(static_cast<unsigned>(v));
                prof[x] += c;
                mass += c;
                px += (v & B_XP) ? 1 : 0;
                px -= (v & B_XM) ? 1 : 0;
                py += (v & B_YP) ? 1 : 0;
                py -= (v & B_YM) ? 1 : 0;
                pz += (v & B_ZP) ? 1 : 0;
                pz -= (v & B_ZM) ? 1 : 0;
                if (c > max_rho) max_rho = c;
            }
        }
    }

    std::vector<double> buf(5 + static_cast<size_t>(G.nx));
    buf[0] = mass;
    buf[1] = px;
    buf[2] = py;
    buf[3] = pz;
    buf[4] = max_rho;
    for (int x = 0; x < G.nx; ++x) buf[5 + x] = prof[x];

    std::vector<double> gprof(pid == 0 ? g_cfg.N : 0, 0.0);
    double gmass = 0.0, gpx = 0.0, gpy = 0.0, gpz = 0.0, gmax = 0.0;

    if (pid == 0) {
        gmass = buf[0];
        gpx = buf[1];
        gpy = buf[2];
        gpz = buf[3];
        gmax = buf[4];
        for (int x = 0; x < G.nx; ++x) gprof[G.x0 + x] += buf[5 + x];
    } else {
        int tag = pid;
        bsp_send(0, &tag, buf.data(), buf.size() * sizeof(double));
    }

    bsp_sync();

    if (pid == 0) {
        int packets = 0;
        size_t nbytes = 0;
        bsp_qsize(&packets, &nbytes);
        const int base = g_cfg.N / P;
        const int rem = g_cfg.N % P;
        for (int i = 0; i < packets; ++i) {
            int size = -1, tag = -1;
            bsp_get_tag(&size, &tag);
            if (size <= 0) break;
            const int payload = size - TAG_BYTES;
            const int nx = base + (tag < rem ? 1 : 0);
            const int x0 = tag * base + std::min(tag, rem);
            if (payload != static_cast<int>((5 + nx) * sizeof(double))) {
                std::fprintf(stderr, "profile: payload inesperado %d\n", payload);
                continue;
            }
            std::vector<double> rbuf(5 + static_cast<size_t>(nx));
            bsp_move(rbuf.data(), static_cast<size_t>(payload));
            gmass += rbuf[0];
            gpx += rbuf[1];
            gpy += rbuf[2];
            gpz += rbuf[3];
            gmax = std::max(gmax, rbuf[4]);
            for (int x = 0; x < nx; ++x) gprof[x0 + x] += rbuf[5 + x];
        }

        const double norm = 1.0 / (static_cast<double>(g_cfg.N) * g_cfg.N);
        for (int x = 0; x < g_cfg.N; ++x) gprof[x] *= norm;

        const double sum_prof = gmass * norm;
        double com = 0.0, rms = 0.0, right = 0.0;
        int front = 0;
        const double iface = g_cfg.iface_frac * g_cfg.N;
        const double front_thr =
            g_cfg.rho_l + 0.1 * (g_cfg.rho_h - g_cfg.rho_l);
        for (int x = 0; x < g_cfg.N; ++x) {
            com += x * gprof[x];
            if (x >= iface) right += gprof[x];
            if (gprof[x] > front_thr) front = x;
        }
        if (sum_prof > 0.0) {
            com /= sum_prof;
            for (int x = 0; x < g_cfg.N; ++x) {
                const double d = x - com;
                rms += d * d * gprof[x];
            }
            rms = std::sqrt(rms / sum_prof);
        }

        if (fp == nullptr) return;

        const double dx_mm = physical_dx_mm();
        const double dt_us = physical_dt_us();
        const double right_frac = (sum_prof > 0.0) ? (right / sum_prof) : 0.0;

        std::fprintf(fp,
                     "%d,%.6g,%.10g,%.10g,%.10g,%.10g,%.6g,%.6g,%.6g,%.6g,%.6g",
                     step, step * dt_us, gmass, gpx, gpy, gpz, gmax,
                     com * dx_mm, rms * dx_mm, front * dx_mm, right_frac);
        for (int x = 0; x < g_cfg.N; ++x) {
            std::fprintf(fp, ",%.6g", gprof[x]);
        }
        std::fprintf(fp, "\n");
    }
}

// ---------------------------------------------------------------------------
// Driver
// ---------------------------------------------------------------------------
void bsp_main()
{
    bsp_begin(bsp_nprocs());
    omp_set_num_threads(g_cfg.threads);

    const int P = bsp_nprocs();
    const int pid = bsp_pid();

    if (P > g_cfg.N) {
        if (pid == 0) {
            std::fprintf(stderr, "P=%d no puede ser mayor que N=%d\n", P,
                         g_cfg.N);
        }
        bsp_abort("P > N");
    }

    const int base = g_cfg.N / P;
    const int rem = g_cfg.N % P;
    const int nx = base + (pid < rem ? 1 : 0);
    const int x0 = pid * base + std::min(pid, rem);

    Grid G;
    G.init(nx, x0, g_cfg.N);

    std::vector<uint8_t> state(G.cells, 0), tmp(G.cells, 0);
    std::vector<uint8_t> ghost_lo(G.layer, 0), ghost_hi(G.layer, 0);

    size_t tag_bytes = static_cast<size_t>(TAG_BYTES);
    bsp_set_tagsize(&tag_bytes);
    bsp_push_reg(ghost_lo.data(), G.layer);
    bsp_push_reg(ghost_hi.data(), G.layer);
    bsp_sync();

    init_state(G, state);
    const int left = (pid - 1 + P) % P;
    const int right = (pid + 1) % P;

    Stats s0 = gather_stats(G, state, 0.0, 0.0);
    if (pid == 0) {
        std::printf("N=%d P=%d T=%d steps=%d geom=%s\n", g_cfg.N, P,
                    g_cfg.threads, g_cfg.steps, g_cfg.geom.c_str());
        std::printf("mass0=%llu momentum0=(%lld,%lld,%lld)\n", s0.mass, s0.px,
                    s0.py, s0.pz);
        std::fflush(stdout);

        print_physical_mapping(select_gas(g_cfg.gas.c_str()), s0.mass);
    }

    if (g_cfg.out_every >= 0) {
        gather_and_write(G, state, 0, g_cfg);
    }

    FILE *profile_fp = nullptr;
    if (!g_cfg.profile_csv.empty()) {
        if (pid == 0) {
            profile_fp = std::fopen(g_cfg.profile_csv.c_str(), "w");
            if (!profile_fp) {
                std::perror(g_cfg.profile_csv.c_str());
            } else {
                std::fprintf(profile_fp,
                             "# gas=%s dx_mm=%.6g dt_us=%.6g N=%d rho_h=%.6g "
                             "rho_l=%.6g iface_frac=%.6g\n",
                             g_cfg.gas.c_str(), physical_dx_mm(),
                             physical_dt_us(), g_cfg.N, g_cfg.rho_h,
                             g_cfg.rho_l, g_cfg.iface_frac);
                std::fprintf(profile_fp,
                             "step,time_us,mass,mx,my,mz,max_rho,com_x_mm,"
                             "rms_x_mm,front_x_mm,right_frac");
                for (int x = 0; x < g_cfg.N; ++x) {
                    std::fprintf(profile_fp, ",rho_%d", x);
                }
                std::fprintf(profile_fp, "\n");
            }
        }
        sample_profile(G, state, 0, profile_fp);
    }

    double t_comp = 0.0, t_comm = 0.0;
    bsp_sync();
    const double t0 = bsp_time();

    for (int step = 0; step < g_cfg.steps; ++step) {
        const double a = bsp_time();
        collide(G, state, tmp, step);
        fill_yz_ghosts(G, tmp);
        const double b = bsp_time();

        exchange_x(G, tmp, ghost_lo, ghost_hi, left, right);
        const double c = bsp_time();

        stream(G, tmp, state);
        const double d = bsp_time();

        t_comp += (b - a) + (d - c);
        t_comm += (c - b);

        if (!g_cfg.profile_csv.empty()
            && ((step + 1) % g_cfg.profile_every == 0)) {
            sample_profile(G, state, step + 1, profile_fp);
        }

        if (g_cfg.out_every > 0 && (step + 1) % g_cfg.out_every == 0
            && step + 1 < g_cfg.steps) {
            gather_and_write(G, state, step + 1, g_cfg);
        }
    }

    bsp_sync();
    const double t_loop = bsp_time() - t0;

    Stats s1 = gather_stats(G, state, t_comp, t_comm);
    if (pid == 0) {
        const unsigned long long dm = (s1.mass > s0.mass)
                                          ? (s1.mass - s0.mass)
                                          : (s0.mass - s1.mass);
        std::printf("mass1=%llu (delta=%llu)\n", s1.mass, dm);
        std::printf("momentum1=(%lld,%lld,%lld)\n", s1.px, s1.py, s1.pz);
        std::printf("time_loop=%.6f s compute_max=%.6f s comm_max=%.6f s\n",
                    t_loop, s1.t_comp, s1.t_comm);
        std::fflush(stdout);
    }

    if (g_cfg.out_every >= 0) {
        gather_and_write(G, state, g_cfg.steps, g_cfg);
    }

    if (pid == 0 && profile_fp != nullptr) {
        std::fclose(profile_fp);
    }

    bsp_pop_reg(ghost_hi.data());
    bsp_pop_reg(ghost_lo.data());
    bsp_sync();
    bsp_end();
}

void usage(const char *prog)
{
    std::printf(
        "Uso: %s [N] [pasos] [hilos] [cada_cuanto_vtk] [prefijo]\n"
        "  N              lado de la malla cubica (default 256)\n"
        "  pasos          pasos de simulacion (default 1000)\n"
        "  hilos          hilos OpenMP por proceso (default 1)\n"
        "  cada_cuanto    VTK cada n pasos; 0 = solo inicio/final;\n"
        "                 <0 = sin VTK (default 0)\n"
        "  prefijo        prefijo de los archivos .vti (default gas)\n"
        "\n"
        "Entorno: GAS=neon|helium, DX_MM, T0_K, P0_ATM,\n"
        "         GEOM=sphere|chamber, RHO_H, RHO_L, P_RATIO,\n"
        "         IFACE_FRAC, RAMP\n"
        "Ejemplo: OMP_NUM_THREADS=4 mpirun -np 4 --bind-to none %s 256 1000 4 100 gas\n"
        "Ejemplo camara: GEOM=chamber P_RATIO=10 RAMP=8 OMP_NUM_THREADS=4 \\\n"
        "                mpirun -np 4 --bind-to none %s 256 1000 4 100 gas\n",
        prog, prog, prog);
}

bool parse_args(int argc, char **argv, Config &cfg)
{
    int i = 1;
    if (i < argc && (std::strcmp(argv[i], "-h") == 0
                     || std::strcmp(argv[i], "--help") == 0)) {
        return false;
    }
    if (i < argc) cfg.N = std::atoi(argv[i++]);
    if (i < argc) cfg.steps = std::atoi(argv[i++]);
    if (i < argc) cfg.threads = std::atoi(argv[i++]);
    if (i < argc) cfg.out_every = std::atoi(argv[i++]);
    if (i < argc) cfg.prefix = argv[i++];

    return cfg.N >= 8 && cfg.N <= 1024 && cfg.steps >= 0 && cfg.threads >= 1;
}

}  // namespace

int main(int argc, char **argv)
{
    if (!parse_args(argc, argv, g_cfg)) {
        usage(argv[0]);
        return 1;
    }

    load_env(g_cfg);
    build_collision_table();

    bsp_init(bsp_main, argc, argv);
    bsp_main();
    return 0;
}
