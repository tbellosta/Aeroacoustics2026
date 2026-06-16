//============================================================
//
//      Type:        FWH include file
//
//      Author:      Tommaso Bellosta on 29/04/26.
//                   Dipartimento di Scienze e Tecnologie Aerospaziali
//                   Politecnico di Milano
//                   Via La Masa 34, 20156 Milano, ITALY
//                   e-mail: tommaso.bellosta@polimi.it
//
//      Copyright:   2026, authors above and the FWH contributors.
//                   This software is distributed under the MIT license, see LICENSE.txt
//
//============================================================
#include "FWHSolver.h"
#include "DataProvider.h"

#include <fstream>
#include <iostream>
#include <numeric>
#include <mpi.h>

int main(int argc, char** argv) {

  MPI_Init(&argc, &argv);

  int rank, size;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);

  /** monopole characteristics (and FWH surface) **/
  double A = 1.0;
  double frq = 10.0;
  double r = 1.0;
  size_t nRefine = 2;

  double src_dt = 1e-4;
  size_t nsteps = 12000;
  double M0_x = 0.9;

  Vect3 M0 = Vect3(M0_x, 0, 0);

  /** This is the data provider **/
  AnaliticalDipoleProvider provider(A,frq,src_dt,nsteps,r,nRefine,M0_x);

  /** This stores the 3 snapshots and cycles them **/
  SnapshotBuffer buffer;
  buffer.initialize(provider);

  /** the solver object **/
  FWHSolver solver;
  solver.permeable = true; // tells the solver we want the permeable surface formulation
  solver.M0_flow  = M0; // mean convection

  /** define a microphone array **/
  double R = 10.0;
  size_t nMicsGlobal = 72;
  size_t nMicsLocal;

  /** define small angular offset so that mics are not aligned with x-axis **/
  double off_deg = 1.0;
  double off_rad = off_deg * M_PI / 180.0;

  /** MPI split the mics in between ranks **/
  size_t chunk_size = nMicsGlobal / size;
  size_t iMicStart, iMicEnd;
  iMicStart = chunk_size*rank;
  iMicEnd = (rank == size - 1) ? nMicsGlobal : iMicStart + chunk_size;

  /** initializes the mics positions and store the observer into the solver object **/
  for (size_t m = iMicStart; m < iMicEnd; m++) {
    double theta = off_rad + 2.0 * M_PI * m / nMicsGlobal;
    Observer obs;
    obs.position = Vect3(R * std::cos(theta),
                         R * std::sin(theta),
                         0.0);
    obs.dt = 10 * src_dt;

    solver.observers.push_back(obs);

  }
  nMicsLocal = solver.observers.size();

  double t_source_fist = 0.0;
  double t_source_end = (nsteps - 1) * src_dt;

  solver.initialize_observers(buffer.curr,t_source_fist,t_source_end);

  /** execute the FWH solver **/

  do {

    solver.process(buffer.prev,
                   buffer.curr,
                   buffer.next,
                   buffer.t_curr,
                   src_dt);


  } while(buffer.advance(provider));


  /** postprocessing the result **/
    /** save both the time histories at the mic positions as well as the directivity **/
  double p_ref = 2.0e-5;

  /** need to communicate data needed for printing the directivity file.
   * MPI_Gather **/

  size_t nFields = 5;
  std::vector<double> data_local(nMicsLocal*nFields);


  for (size_t iMic = 0; iMic < nMicsLocal; iMic++) {
    /** get the time history at the mic and compute the rms **/
   const auto& obs = solver.observers[iMic];
    std::vector<double> pprime = obs.signalL;
    for (size_t iSampl = 0; iSampl < pprime.size(); iSampl++) {
      pprime[iSampl] += obs.signalT[iSampl];
    }

   double rms_fwh = 0, mean_fwh = 0;
   size_t len = 0;

   /** compute the signal mean **/
   for (size_t i = obs.trim_start; i < obs.trim_end; i++) {
     mean_fwh += pprime[i];
     ++len;
   }
   mean_fwh /= len;

   for (size_t i = obs.trim_start; i < obs.trim_end; i++)
    rms_fwh += std::pow(pprime[i] - mean_fwh,2);

   rms_fwh = std::sqrt(rms_fwh /= len);
   double spl_fwh = 20.0 * std::log10(std::max(rms_fwh,1e-30) / p_ref);

   /** now compute the exact SPL and rms **/
   double sum2_exact = 0.0, mean_exact = 0.0;
   for (size_t i = obs.trim_start; i < obs.trim_end; i++) {
     double t = obs.t0 + i * obs.dt;
     double p = provider.exact_presure_at(obs.position[0],
                                          obs.position[1],
                                          obs.position[2],
                                          t);
     mean_exact += p;
   }
   mean_exact /= len;

    for (size_t i = obs.trim_start; i < obs.trim_end; i++) {
      double t = obs.t0 + i * obs.dt;
      double p = provider.exact_presure_at(obs.position[0],
                                           obs.position[1],
                                           obs.position[2],
                                           t);
      sum2_exact += std::pow(p - mean_exact,2);
    }
    sum2_exact /= len;

    double rms_exact = std::sqrt(sum2_exact);
    double spl_exact = 20.0 * std::log10(std::max(rms_exact,1e-30) / p_ref);

    /** compute mic position **/
    size_t iMicGlobal = iMicStart + iMic;
    double theta_deg = off_deg + 360.0 * iMicGlobal / nMicsGlobal;

    /** store for later comm **/
    data_local[iMic*nFields + 0] = theta_deg;
    data_local[iMic*nFields + 1] = spl_fwh;
    data_local[iMic*nFields + 2] = rms_fwh;
    data_local[iMic*nFields + 3] = spl_exact;
    data_local[iMic*nFields + 4] = rms_exact;

  }

  /** now gather all data to rank0 for printing **/
  std::vector<double> data_global;

  if (rank == 0) {
    data_global.resize(nMicsGlobal*nFields);
  }

  /** MPI_Gather expects an equal amount of data
   * from each rank. In our program, the last rank
   * may have extra elements. We could either use the
   * function MPI_Gatherv which works for variable elements
   * from all ranks, or simply use MPI_Gather and then
   * send the possible extra mics from rank(size-1) to
   * rank0. We do the latter.
   * **/

  /** Those are the extra mics in the last rank **/
  const int remainder = nMicsGlobal % size;

    MPI_Gather(
          data_local.data(),
          static_cast<int>(chunk_size*nFields),
          MPI_DOUBLE,
          data_global.data(),
          static_cast<int>(chunk_size*nFields),
          MPI_DOUBLE,
          0,
          MPI_COMM_WORLD);


  if (remainder > 0) {
    if (rank == size - 1)
        MPI_Send(data_local.data()+chunk_size*nFields,
             remainder*nFields,
                   MPI_DOUBLE,
             0,
             0,
             MPI_COMM_WORLD);

    if (rank == 0)
        MPI_Recv(data_global.data()+(size*chunk_size*nFields),
        remainder*nFields,
        MPI_DOUBLE,
        size-1,
        0,
        MPI_COMM_WORLD,
        MPI_STATUS_IGNORE);

  }

  if (rank == 0) {
    std::ofstream out("directivity.dat");
    out << "# theta[deg] SPL_FWH[dB] p_rms_FWH[Pa] SPL_exact[dB] p_rms_exact[Pa]\n";

    for (int iMic = 0; iMic < nMicsGlobal; ++iMic) {
      out << data_global[iMic*nFields + 0] << " "
          << data_global[iMic*nFields + 1] << " "
          << data_global[iMic*nFields + 2] << " "
          << data_global[iMic*nFields + 3] << " "
          << data_global[iMic*nFields + 4] << "\n";
    }

    /** replicate fist mic to close directivity polar plot **/
    out << data_global[0*nFields + 0] << " "
        << data_global[0*nFields + 1] << " "
        << data_global[0*nFields + 2] << " "
        << data_global[0*nFields + 3] << " "
        << data_global[0*nFields + 4] << "\n";

    out.close();
  }



  for (size_t iMic = 0; iMic < nMicsLocal; iMic++) {

    const Observer& obs = solver.observers[iMic];

    /** get the global mic ID **/
    size_t iMicGlobal = iMicStart + iMic;

    /** write time history **/
    int width = std::to_string(nMicsGlobal).size();
    std::string num = std::to_string(iMicGlobal);
    std::string padded = std::string(width - num.length(),'0') + num;
    std::string fname = "signal_" + padded + ".dat";
    std::ofstream signal(fname);
    signal << "# t[s] p_FWH_l[Pa] p_FWH_t[Pa] p_exact[Pa]\n";

    for (size_t i = obs.trim_start; i < obs.trim_end; i++) {
      double t = obs.t0 + i * obs.dt;
      double p_fwh_l = obs.signalL[i];
      double p_fwh_t = obs.signalT[i];
      double p_exact = provider.exact_presure_at(obs.position[0],
                                                 obs.position[1],
                                                 obs.position[2],
                                                 t);
      signal << t << " " << p_fwh_l << " " << p_fwh_t << " " << p_exact << "\n";
    }
    signal.close();

  }

  MPI_Finalize();

}