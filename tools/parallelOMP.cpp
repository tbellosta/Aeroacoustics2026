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
#include <omp.h>
#include <iostream>

#include <iostream>
#include <vector>
#include <omp.h>

int main() {

    std::cout << "Using " << omp_get_max_threads() << " threads\n";

    // ============================================================
    // STEP 1: RACE CONDITION (shared variable)
    // ============================================================

    {
        std::cout << "\nSTEP 1: Race condition\n";

        int x = 0;

        #pragma omp parallel
        {
            for (int i = 0; i < 1000; ++i) {
                x++;  // race condition
            }
        }

        std::cout << "Final x (wrong) = " << x << std::endl;
    }


    // ============================================================
    // STEP 2: PRIVATE VARIABLES (thread-local copies)
    // ============================================================

    {
        std::cout << "\nSTEP 2: Private copies\n";

        int x = 0;

        #pragma omp parallel private(x)
        {
            x = 0; // each thread has its own copy

            for (int i = 0; i < 1000; ++i) {
                x++;
            }

            #pragma omp critical
            std::cout << "Thread " << omp_get_thread_num()
                      << " local x = " << x << std::endl;
        }

        std::cout << "Global x (unchanged) = " << x << std::endl;
    }


    // ============================================================
    // STEP 3: REDUCTION (correct parallel accumulation)
    // ============================================================

    {
        std::cout << "\nSTEP 3: Reduction\n";

        int x = 0;

        #pragma omp parallel for reduction(+:x)
        for (int i = 0; i < 1000; ++i) {
            x++;
        }

        std::cout << "Final x (correct) = " << x << std::endl;
    }


    // ============================================================
    // STEP 4: MANUAL REDUCTION (what OpenMP does internally)
    // ============================================================

    {
        std::cout << "\nSTEP 4: Manual reduction\n";

        int nthreads = omp_get_max_threads();
        std::vector<int> local(nthreads, 0);

           #pragma omp parallel for
            for (int i = 0; i < 1000; ++i) {
                int tid = omp_get_thread_num();
                local[tid]++;
            }

        int global = 0;
        for (int i = 0; i < nthreads; ++i) {
            global += local[i];
        }

        std::cout << "Final sum = " << global << std::endl;
    }


    // ============================================================
    // STEP 5: MEMORY VISIBILITY PROBLEM
    // ============================================================

    {
        std::cout << "\nSTEP 5: Memory visibility issue (may hang!)\n";

        int flag = 0;

        #pragma omp parallel num_threads(2) shared(flag)
        {

            if (omp_get_thread_num() == 0) {
                flag = 1;
            } else {
                while (flag == 0) {
                    // spin wait
                }
                std::cout << "Thread sees flag = " << flag << std::endl;
            }
        }
    }


    // ============================================================
    // STEP 6: FIX WITH FLUSH
    // ============================================================

    {
        std::cout << "\nSTEP 6: Fix with flush\n";

        int flag = 0;

        #pragma omp parallel num_threads(2) shared(flag)
        {

            if (omp_get_thread_num() == 0) {
                flag = 1;
                #pragma omp flush(flag)
            } else {
                while (true) {
                    #pragma omp flush(flag)
                    if (flag == 1) break;
                }
                std::cout << "Thread sees flag = " << flag << std::endl;
            }
        }
    }


    // ============================================================
    // STEP 7: FIX WITH BARRIER (cleaner)
    // ============================================================

    {
        std::cout << "\nSTEP 7: Barrier synchronization\n";

        int x = 0;

        #pragma omp parallel shared(x)
        {
            int tid = omp_get_thread_num();

            if (tid == 0) {
                x = 42;
            }

            #pragma omp barrier
            #pragma omp critical
            std::cout << "Thread " << tid
                      << " sees x = " << x << std::endl;
        }
    }


    // ============================================================
    // STEP 8: ATOMIC VS CRITICAL
    // ============================================================

    {
        std::cout << "\nSTEP 8: atomic vs critical\n";

        int x = 0;

        #pragma omp parallel
        {
            for (int i = 0; i < 1000; ++i) {

                // Try switching between these:

                 #pragma omp atomic
                 x++;

//                #pragma omp critical
//                x++;
            }
        }

        std::cout << "Final x = " << x << std::endl;
    }



}

