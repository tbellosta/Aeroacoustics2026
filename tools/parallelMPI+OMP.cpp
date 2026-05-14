//============================================================
//
//      Type:        FWH include file
//
//      Author:      Tommaso Bellosta on 06/05/26.
//                   Dipartimento di Scienze e Tecnologie Aerospaziali
//                   Politecnico di Milano
//                   Via La Masa 34, 20156 Milano, ITALY
//                   e-mail: tommaso.bellosta@polimi.it
//
//      Copyright:   2026, authors above and the FWH contributors.
//                   This software is distributed under the MIT license, see LICENSE.txt
//
//============================================================
#include <iostream>
#include <omp.h>
#include <mpi.h>

int main(int argc, char** argv) {

    MPI_Init(&argc, &argv);

    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    // Total problem size
    const int N = 10000000;

    // Domain decomposition (block)
    int chunk = N / size;
    int start = rank * chunk;


    int end   = (rank == size - 1) ? N : start + chunk;

    /** same as ternary operator above **/
    // if (rank == size - 1) end = N;
    // else end = start + chunk;

    double local_sum = 0.0;

    // ============================================
    // OpenMP parallel region inside each MPI rank
    // ============================================
#pragma omp parallel for reduction(+:local_sum)
    for (int i = start; i < end; ++i) {
        local_sum += 1.0;  // trivial workload
    }

    // Reduce across MPI processes
    double global_sum = 0.0;

    MPI_Reduce(&local_sum, &global_sum, 1,
               MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

    if (rank == 0) {
        std::cout << "Global sum = " << global_sum << std::endl;
    }

    // Debug print (optional)
    std::cout << "Rank " << rank
              << " used " << omp_get_max_threads()
              << " threads, local sum = " << local_sum << std::endl;

    MPI_Finalize();
    return 0;
}