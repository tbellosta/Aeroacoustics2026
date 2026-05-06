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
#include <mpi.h>

int main(int argc, char** argv) {

    MPI_Init(&argc, &argv);

    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    if (rank == 0) {
        std::cout << "Running with " << size << " processes\n";
    }
    MPI_Barrier(MPI_COMM_WORLD);

    // ============================================================
    // STEP 1: EACH PROCESS HAS ITS OWN MEMORY
    // ============================================================
    if (rank == 0) std::cout << "\nSTEP 1\n";
    {
        int x = 0;

        x += rank;

        std::cout << "Rank " << rank
                  << " has x = " << x  << "\n";
        MPI_Barrier(MPI_COMM_WORLD);

    }



    // ============================================================
    // STEP 2: NO AUTOMATIC SHARING
    // ============================================================
//    if (rank == 0) std::cout << "\nSTEP 2\n";
//    {
//        int x = 0;
//
//        if (rank == 0) {
//            x = 42;
//        }
//
//        // No communication here!
//
//        std::cout << "Rank " << rank
//                  << " sees x = " << x << std::endl;
//        MPI_Barrier(MPI_COMM_WORLD);
//    }



    // ============================================================
    // STEP 3: EXPLICIT COMMUNICATION (SEND/RECV)
    // ============================================================
//    if (rank == 0) std::cout << "\nSTEP 3\n";
//    {
//        int x = 0;
//
//        if (rank == 0) {
//            x = 42;
//            MPI_Send(&x, 1, MPI_INT, 1, 0, MPI_COMM_WORLD);
//            std::cout << "Rank 0 sent x = " << x << std::endl;
//        }
//
//        if (rank == 1) {
//            MPI_Recv(&x, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
//            std::cout << "Rank 1 received x = " << x << std::endl;
//        }
//        MPI_Barrier(MPI_COMM_WORLD);
//    }



    // ============================================================
    // STEP 4: BROADCAST (ONE-TO-ALL COMMUNICATION)
    // ============================================================
//    if (rank == 0) std::cout << "\nSTEP 4\n";
//    {
//        int x = 0;
//
//        if (rank == 0) {
//            x = 42;
//        }
//
//        MPI_Bcast(&x, 1, MPI_INT, 0, MPI_COMM_WORLD);
//
//        std::cout << "Rank " << rank
//                  << " now has x = " << x << std::endl;
//        MPI_Barrier(MPI_COMM_WORLD);
//    }




    // ============================================================
    // STEP 5: LOCAL WORK + GLOBAL REDUCTION
    // ============================================================
//    if (rank == 0) std::cout << "\nSTEP 5\n";
//    {
//        int local = rank + 1;  // each process has different data
//
//        int global = 0;
//
//        MPI_Reduce(&local, &global, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
//
//        std::cout << "Rank " << rank
//                  << " local = " << local << std::endl;
//
//        if (rank == 0) {
//            std::cout << "Global sum = " << global << std::endl;
//        }
//        MPI_Barrier(MPI_COMM_WORLD);
//    }



    // ============================================================
    // STEP 6: ALLREDUCE (EVERYONE GETS RESULT)
    // ============================================================
//    if (rank == 0) std::cout << "\nSTEP 6\n";
//    {
//        int local = rank + 1;
//
//        int global = 0;
//
//        MPI_Allreduce(&local, &global, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
//
//        std::cout << "Rank " << rank
//                  << " sees global sum = " << global << std::endl;
//        MPI_Barrier(MPI_COMM_WORLD);
//    }


    MPI_Finalize();
    return 0;
}

