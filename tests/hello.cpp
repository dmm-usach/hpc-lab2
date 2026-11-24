// Prueba de humo: BSPlib (BSPonMPI) + OpenMP
extern "C" {
#include <bsp.h>
}

#include <omp.h>
#include <cstdio>
#include <cstdlib>

static int nthreads = 1;

void bsp_main()
{
    bsp_begin(bsp_nprocs());
    omp_set_num_threads(nthreads);

    const int pid = bsp_pid();
    const int P = bsp_nprocs();

    int inside = 0;
    #pragma omp parallel
    {
        #pragma omp single
        inside = omp_get_num_threads();
    }

    printf("pid %d/%d: %d hilos OpenMP\n", pid, P, inside);
    fflush(stdout);

    bsp_end();
}

int main(int argc, char **argv)
{
    nthreads = (argc > 1) ? atoi(argv[1]) : 1;
    bsp_init(bsp_main, argc, argv);
    bsp_main();
    return 0;
}
