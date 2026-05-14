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

#include "meshStructure.h"

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
    // if (rank == 0) std::cout << "\nSTEP 1\n";
    // {
    //     int x = 0;
    //
    //     x += rank;
    //
    //     std::cout << "Rank " << rank
    //               << " has x = " << x  << "\n";
    //     MPI_Barrier(MPI_COMM_WORLD);
    //
    // }



    // ============================================================
    // STEP 2: NO AUTOMATIC SHARING
    // ============================================================
    // if (rank == 0) std::cout << "\nSTEP 2\n";
    // {
    //     int x = 0;
    //
    //     if (rank == 0) {
    //         x = 42;
    //     }
    //
    //     // No communication here!
    //
    //     std::cout << "Rank " << rank
    //               << " sees x = " << x << std::endl;
    //     MPI_Barrier(MPI_COMM_WORLD);
    // }



    // ============================================================
    // STEP 3: EXPLICIT COMMUNICATION (SEND/RECV)
    // ============================================================
    // if (rank == 0) std::cout << "\nSTEP 3\n";
    // {
    //     int x = 0;
    //     int y = 0;
    //
    //     if (rank == 0) {
    //         x = 42;
    //         MPI_Send(&x, 1, MPI_INT, 1, 0, MPI_COMM_WORLD);
    //         // std::cout << "Rank 0 sent x = " << x << std::endl;
    //     }
    //
    //     if (rank == 1) {
    //         MPI_Recv(&x, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    //         // std::cout << "Rank 1 received x = " << x << std::endl;
    //     }
    //
    //     std::cout << "Rank " << rank
    //         << " sees x = " << x << std::endl;
    //     MPI_Barrier(MPI_COMM_WORLD);
    // }


    // if (rank == 0) std::cout << "\nSTEP 3\n";
    // {
    //     int x0 = 0;
    //     int x1 = 0;
    //
    //     MPI_Request req_send;
    //     MPI_Request req_recv;
    //
    //     if (rank == 0) {
    //         x0 = 42;
    //         MPI_Recv(&x1, 1, MPI_INT, 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    //         std::cout << "Rank 0 received x1 = " << x1 << std::endl;
    //         MPI_Send(&x0, 1, MPI_INT, 1, 0, MPI_COMM_WORLD);
    //         std::cout << "Rank 0 sent x0 = " << x0 << std::endl;
    //     }
    //
    //     if (rank == 1) {
    //         x1 = 32;
    //         MPI_Send(&x1, 1, MPI_INT, 0, 0, MPI_COMM_WORLD);
    //         std::cout << "Rank 1 sent x1 = " << x1 << std::endl;
    //         MPI_Recv(&x0, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    //         std::cout << "Rank 1 received x0 = " << x0 << std::endl;
    //
    //     }
    //
    //     // MPI_Wait(&req,MPI_STATUS_IGNORE);
    //
    //     MPI_Barrier(MPI_COMM_WORLD);
    // }



    // ============================================================
    // STEP 4: BROADCAST (ONE-TO-ALL COMMUNICATION)
    // ============================================================
    // if (rank == 0) std::cout << "\nSTEP 4\n";
    // {
    //     int x = 0;
    //
    //     if (rank == 0) {
    //         x = 42;
    //     }
    //
    //     MPI_Bcast(&x, 1, MPI_INT, 0, MPI_COMM_WORLD);
    //
    //     std::cout << "Rank " << rank
    //               << " now has x = " << x << std::endl;
    //     MPI_Barrier(MPI_COMM_WORLD);
    // }




    // ============================================================
    // STEP 5: LOCAL WORK + GLOBAL REDUCTION (ALL-TO-ONE comm)
    // ============================================================
    // if (rank == 0) std::cout << "\nSTEP 5\n";
    // {
    //     int local = rank + 1;  // each process has different data
    //
    //     int global = 0;
    //
    //     MPI_Reduce(&local, &global, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
    //
    //     std::cout << "Rank " << rank
    //               << " local = " << local << std::endl;
    //
    //     if (rank == 0) {
    //         std::cout << "Global sum = " << global << std::endl;
    //     }
    //     MPI_Barrier(MPI_COMM_WORLD);
    // }



    // ============================================================
    // STEP 6: ALLREDUCE (EVERYONE GETS RESULT) (ALL-TO-ALL comm)
    // ============================================================
    // if (rank == 0) std::cout << "\nSTEP 6\n";
    // {
    //     int local = rank + 1;
    //
    //     int global = 0;
    //
    //     MPI_Allreduce(&local, &global, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
    //
    //     std::cout << "Rank " << rank
    //               << " sees global sum = " << global << std::endl;
    //     MPI_Barrier(MPI_COMM_WORLD);
    // }

    // ============================================================
    // STEP 7: CUSTOM DATATYPES
    // ============================================================

    /** Only run if size > 2 **/
    if (size < 2) {
        std::cout << "Not enough processes, need at least 2\n";
        MPI_Finalize();
        return 1;
    }

    /** Create the custom type.
     * Ideally this would be a function.
     **/

    MPI_Datatype nodeType;

    /** number of fields in the class we want to comm **/
    constexpr int nMembers = 5;

    /** number of elements for each field **/
    int blockLengths[nMembers] = {
        3, // x
        3, // dS
        1, // p
        3, // u
        1  // rho
    };

    /** offset in memory for each field in the class.
     * This tells MPI where to find a given field given
     * the address in memory of the beginning of a
     * Node object.
     **/
    MPI_Aint offsets[nMembers];

    offsets[0] = offsetof(Node, x);
    offsets[1] = offsetof(Node, dS);
    offsets[2] = offsetof(Node, p);
    offsets[3] = offsetof(Node, u);
    offsets[4] = offsetof(Node, rho);

    /** The types of the fields in Node **/
    MPI_Datatype types[nMembers] = {
        MPI_DOUBLE,
        MPI_DOUBLE,
        MPI_DOUBLE,
        MPI_DOUBLE,
        MPI_DOUBLE
    };

    /** Create the custom datatype **/
    MPI_Type_create_struct(
        nMembers,
        blockLengths,
        offsets,
        types,
        &nodeType
    );

    /** This is mandatory. To use a custom
     * datatype, you need to commit it.
     **/
    MPI_Type_commit(&nodeType);


    /** Now we can use the custom type we created.
     * First we are sending a single node from Rank0
     * to Rank1 **/

    Node toBeSent;  // this is only meaningful on rank 0
    Node toBeRecvd; // this is only meaningful on rank 1

    MPI_Request request;

    if (rank == 0) {
        /** do stuff with the Node you want
         * to send.
         **/
        toBeSent.x = Vect3(1.0, 2.0, 3.0);
        // toBeSent.* = ...

        /** Use a nonblocking comm. The next function
         * will return before the actual communication
         * is complete. Need MPI_Wait if you actually
         * need to perform operations that depend on
         * the result of the comm operation.
         **/
        MPI_Isend(&toBeSent,1,nodeType,1,0,MPI_COMM_WORLD,&request);

    } else if (rank ==1) {
        /** receive the node **/
        MPI_Irecv(&toBeRecvd,1,nodeType,0,0,MPI_COMM_WORLD,&request);

        /** If you want to use it, you need to wait for the
         * comm to complete.
         **/
        MPI_Wait(&request,MPI_STATUS_IGNORE);
        std::cout << "Rank 1 received a node with position:\t"
                  << toBeRecvd.x[0] << " "
                  << toBeRecvd.x[1] << " "
                  << toBeRecvd.x[2] << "\n";

    }


    /** Now we will send a vector of Nodes
     * from Rank0 to Rank1.
     **/

    std::vector<MPI_Request> requests;

    std::vector<Node> nodesToComm; // filled only on rank0
    std::vector<Node> recvNodes; // filled only on rank1

    if (rank == 0) {
        /** This is the vector of nodes we want
         * to send to Rank1. **/
        int nNodes = 1000;
        nodesToComm.resize(nNodes);
        {
            // to something with the nodes
        }

        /** before sending the whole data,
         * we need to comm the number of elements
         * we are sending over to Rank1. Rank1 needs
         * to know how many nodes we are sending to allocate
         * the memory needed to store the whole set of nodes
         * (it simply needs to initialize a vector with the right
         * number of elements).
         * **/

        requests.resize(2); // one for the size, the other for the data

        MPI_Isend(&nNodes,1,MPI_INT,1,0,MPI_COMM_WORLD,&requests[0]);
        MPI_Isend(nodesToComm.data(),nNodes,nodeType,1,0,MPI_COMM_WORLD,&requests[1]);

    } else if (rank == 1) {
        /** we only need one request (we will reuse the same
         * since we can receive the data only after we got the
         * number of nodes).
         * **/
        requests.resize(1);

        int nNodes;

        /** receive the number of nodes to setup the receiving buffer **/
        MPI_Irecv(&nNodes,1,MPI_INT,0,0,MPI_COMM_WORLD,&requests[0]);

        /** We need to wait for the comm to complete.
         * We could have used a blocking call instead.
         * **/
        MPI_Wait(&requests[0],MPI_STATUS_IGNORE);

        /** Now we can prepare the receiving vector **/
        recvNodes.resize(nNodes);
        MPI_Irecv(recvNodes.data(),nNodes,nodeType,0,0,MPI_COMM_WORLD,&requests[0]);
    }

    /** Remember that comms where non-blocking.
     * To use the communicated data, wait for all comms
     * to complete. Rank0 needs to complete 2, Rank1 only one
     * (we already made sure the first recv is complete).
     * **/
    MPI_Waitall(static_cast<int>(requests.size()),requests.data(),MPI_STATUS_IGNORE);
    if (rank == 0) {
        std::cout << "Rank0 sent 1000 nodes to Rank1\n";
    }

    MPI_Finalize();
    return 0;
}

