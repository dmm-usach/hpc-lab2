// Prueba de intercambio de halos 1D con BSPlib (DRMA) + OpenMP
extern "C" {
#include <bsp.h>
}

#include <omp.h>
#include <vector>
#include <cstdio>
#include <cstdlib>
#include <cmath>

static int NY = 64, NZ = 64, NTHREADS = 4;

void bsp_main()
{
    bsp_begin(bsp_nprocs());
    omp_set_num_threads(NTHREADS);

    const int P    = bsp_nprocs();
    const int pid  = bsp_pid();
    const int left  = (pid - 1 + P) % P;
    const int right = (pid + 1) % P;
    const int face = NY * NZ;
    const int bytes = face * (int)sizeof(float);

    std::vector<float> face_lo(face), face_hi(face);
    std::vector<float> ghost_lo(face, -1.0f), ghost_hi(face, -1.0f);

    bsp_push_reg(ghost_lo.data(), bytes);
    bsp_push_reg(ghost_hi.data(), bytes);
    bsp_sync();  // los registros son efectivos tras la barrera

    // Cada rank llena sus caras con un valor identificable
    #pragma omp parallel for collapse(2)
    for (int j = 0; j < NY; ++j)
        for (int k = 0; k < NZ; ++k) {
            const int idx = j * NZ + k;
            face_lo[idx] = 1000.0f * pid + 1.0f;
            face_hi[idx] = 1000.0f * pid + 2.0f;
        }

    // face_lo va al ghost_hi del vecino izquierdo; face_hi al ghost_lo del derecho
    bsp_put(left,  face_lo.data(), ghost_hi.data(), 0, bytes);
    bsp_put(right, face_hi.data(), ghost_lo.data(), 0, bytes);
    bsp_sync();

    // Verificación: ghost_lo debe tener face_hi del vecino izquierdo,
    //               ghost_hi debe tener face_lo del vecino derecho.
    const float expect_lo = 1000.0f * left  + 2.0f;
    const float expect_hi = 1000.0f * right + 1.0f;

    long bad = 0;
    #pragma omp parallel for collapse(2) reduction(+:bad)
    for (int j = 0; j < NY; ++j)
        for (int k = 0; k < NZ; ++k) {
            const int idx = j * NZ + k;
            if (ghost_hi[idx] != expect_hi) ++bad;
            if (ghost_lo[idx] != expect_lo) ++bad;
        }

    printf("pid %d/%d: hilos=%d, errores=%ld (esperaba %g y %g)\n",
           pid, P, omp_get_max_threads(), bad, expect_lo, expect_hi);
    fflush(stdout);

    bsp_pop_reg(ghost_hi.data());
    bsp_pop_reg(ghost_lo.data());
    bsp_sync();
    bsp_end();
}

int main(int argc, char **argv)
{
    if (argc > 1) NTHREADS = atoi(argv[1]);
    bsp_init(bsp_main, argc, argv);
    bsp_main();
    return 0;
}
