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
#include <regex>
#include <sstream>
#include <iomanip>
#include <mpi.h>

/** generate a structured half-sphere to map
* half plane directivity **/
std::vector<Vect3> generateHemisphere(
    double radius,
    int nTheta,
    int nPhi);

/** generate vector of filenames for all input
* FWH surface snapshots **/
std::vector<std::string> generateFileNames(
    const std::string& firstFileName,
    size_t delta, size_t nFiles
);

int main(int argc, char** argv) {

  MPI_Init(&argc, &argv);

  int rank, size;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);

  /** input data characteristics **/
  double src_dt = 0.00006309165952;
  size_t nsteps = 2520;
  double M0_x = 0.0;

  Vect3 M0 = Vect3(M0_x, 0, 0);

  std::string fileBase = "/home/bellosta/Documents/AEROACOUSTICS/HOVER/FWH/bin/surface_flow_01080.vtu";

  /** create a vector with the path to all snapshots **/
  const auto filenames = generateFileNames(fileBase,1,nsteps);

  /** This is the data provider **/
  VTUProvider provider(src_dt, filenames);


  /** This stores the 3 snapshots and cycles them **/
  SnapshotBuffer buffer;
  buffer.initialize(provider);

  /** the solver object **/
  FWHSolver solver;
  solver.permeable = false; // tells the solver we want the solid surface formulation
  solver.M0_flow  = M0; // no mean convection

  /** define a microphone array **/
  double R = 10.0;
  /** full half-sphere directivity **/
//  auto pnts = generateHemisphere(R, 18, 72);
//  size_t nMicsGlobal = pnts.size();
  /** single arch directivity **/
  size_t nMicsGlobal = 72;
  size_t nMicsLocal;

  /** define small angular offset so that mics are not aligned with x-axis **/
  double off_deg = 0.0;
  double off_rad = off_deg * M_PI / 180.0;

  /** MPI split the mics in between ranks **/
  size_t chunk_size = nMicsGlobal / size;
  size_t iMicStart, iMicEnd;
  iMicStart = chunk_size*rank;
  iMicEnd = (rank == size - 1) ? nMicsGlobal : iMicStart + chunk_size;

  /** initializes the mics positions and store the observer into the solver object **/
  for (size_t m = iMicStart; m < iMicEnd; m++) {
    Observer obs;

    /** uncomment for simple microphone arch in xz plane **/
    double theta = off_rad + 2.0 * M_PI * m / nMicsGlobal;
    obs.position = Vect3(R * std::cos(theta),
                         0.0,
                         R * std::sin(theta));

    /** This computes the whole half sphere directivity **/
//    obs.position = pnts[m];

    obs.dt = 1 * src_dt;
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
  /** save both the time histories at the mic positions
   *  as well as the directivity **/
  double p_ref = 2.0e-5;

  /** need to communicate data needed for printing the directivity file.
   * we will use MPI_Gather **/

  /** flatten the data to be communicated **/
  size_t nFields = 3;
  std::vector<double> data_local(nMicsLocal*nFields);


  for (size_t iMic = 0; iMic < nMicsLocal; iMic++) {
   /** get the time history at the mic and compute the rms **/
   const auto& obs = solver.observers[iMic];

   double rms_fwh = 0, mean_fwh = 0;
   size_t len = 0;

   /** compute the signal mean **/
   for (size_t i = obs.trim_start; i < obs.trim_end; i++) {
     mean_fwh += obs.signal[i];
     ++len;
   }
   mean_fwh /= len;

   for (size_t i = obs.trim_start; i < obs.trim_end; i++)
    rms_fwh += std::pow(obs.signal[i] - mean_fwh,2);

   rms_fwh = std::sqrt(rms_fwh /= len);
   double spl_fwh = 20.0 * std::log10(std::max(rms_fwh,1e-30) / p_ref);

    /** compute mic position **/
    size_t iMicGlobal = iMicStart + iMic;
    double theta_deg = off_deg + 360.0 * iMicGlobal / nMicsGlobal;

    /** store for later comm **/
    data_local[iMic*nFields + 0] = theta_deg;
    data_local[iMic*nFields + 1] = spl_fwh;
    data_local[iMic*nFields + 2] = rms_fwh;

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
    out << "# theta[deg] SPL_FWH[dB] p_rms_FWH[Pa]\n";

    for (int iMic = 0; iMic < nMicsGlobal; ++iMic) {
      out << data_global[iMic*nFields + 0] << " "
          << data_global[iMic*nFields + 1] << " "
          << data_global[iMic*nFields + 2] << "\n";
    }

    /** replicate fist mic to close directivity polar plot **/
    out << data_global[0*nFields + 0] << " "
        << data_global[0*nFields + 1] << " "
        << data_global[0*nFields + 2] << "\n";

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
    signal << "# t[s] p_FWH[Pa]\n";

    for (size_t i = obs.trim_start; i < obs.trim_end; i++) {
      double t = obs.t0 + i * obs.dt;
      double p_fwh = obs.signal[i];

      signal << t << " " << p_fwh << "\n";
    }
    signal.close();

  }

  MPI_Finalize();

}


std::vector<Vect3> generateHemisphere(
    double radius,
    int nTheta,
    int nPhi)
{
  std::vector<Vect3> points;

  const double pi = std::acos(-1.0);

  for (int i = 0; i < nTheta; ++i) {
    // theta in [0, pi/2]
    double theta =
        (static_cast<double>(i) / (nTheta - 1))
        * (pi / 2.0);

    for (int j = 0; j < nPhi; ++j) {
      // phi in [0, 2pi)
      double phi =
          (static_cast<double>(j) / nPhi)
          * (2.0 * pi);

      Vect3 p;
      p[0] = radius * std::sin(theta) * std::cos(phi);
      p[1] = radius * std::sin(theta) * std::sin(phi);
      p[2] = radius * std::cos(theta);

      points.push_back(p);
    }
  }

  return points;
}


std::vector<std::string> generateFileNames(
    const std::string& firstFileName,
    size_t delta, size_t nFiles)
{

    std::regex pattern(R"((.*?)(\d+)(\.[^.]+)$)");
    std::smatch match;

    std::vector<std::string> filenames;

    if (std::regex_match(firstFileName, match, pattern)) {

        std::string prefix = match[1];      // "/home/user/.../surface_flow_"
        std::string numberStr = match[2];   // "01080"
        std::string suffix = match[3];      // ".vtu"

        int start = std::stoi(numberStr);
        int width = numberStr.size(); // preserve leading zeros

        for (int i = 0; i < nFiles; ++i) {
            std::ostringstream oss;
            oss << prefix
                << std::setw(width)
                << std::setfill('0')
                << (start + i * delta)
                << suffix;

            filenames.push_back(oss.str());
        }
    }

    return filenames;
  }
