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
#include <numeric>

int main() {

  /** monopole characteristics (and FWH surface) **/
  double A = 1.0;
  double frq = 10.0;
  double r = 1.0;
  size_t nRefine = 2;

  double src_dt = 1e-4;
  size_t nsteps = 12000;

  /** This is the data provider **/
  AnaliticalDipoleProvider provider(A,frq,src_dt,nsteps,r,nRefine);

  /** This stores the 3 snapshots and cycles them **/
  SnapshotBuffer buffer;
  buffer.initialize(provider);

  /** the solver object **/
  FWHSolver solver;
  solver.permeable = true; // tells the solver we want the permeable surface formulation
  solver.M0_flow  = Vect3(0.0, 0.0, 0.0); // no mean convection

  /** define a microphone array **/
  double R = 2.0;
  size_t nMics = 72;

  /** define small angular offset so that mics are not aligned with x-axis **/
  double off_deg = 1.0;
  double off_rad = off_deg * M_PI / 180.0;

  /** initializes the mics positions and store the observer into the solver object **/
  for (size_t m = 0; m < nMics; m++) {
    double theta = off_rad + 2.0 * M_PI * m / nMics;
    Observer obs;
    obs.position = Vect3(R * std::cos(theta),
                         R * std::sin(theta),
                         0.0);
    obs.dt = 10 * src_dt;

    solver.observers.push_back(obs);

  }

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

  std::ofstream out("directivity.dat");
  out << "# theta[deg] SPL_FWH[dB] p_rms_FWH[Pa] SPL_exact[dB] p_rms_exact[Pa]\n";

  /** create structure to circle over all mics plus tmic 0 again **/
  std::vector<size_t> loop(nMics+1,0);
  std::iota(loop.begin(),loop.end(),0);
  loop[nMics] = 0;



  for (size_t iMic : loop) {
    /** get the time history at the mic and compute the rms **/
   auto& obs = solver.observers[iMic];

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
   double SPL_fwh = 20.0 * std::log10(std::max(rms_fwh,1e-30) / p_ref);

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
    double SPL_exact = 20.0 * std::log10(std::max(rms_exact,1e-30) / p_ref);

    /** compute mic position **/
    double theta_deg = off_deg + 360.0 * iMic / nMics;

    /** print to file **/
    out << theta_deg << " "
        << SPL_fwh   << " " << rms_fwh   << " "
        << SPL_exact << " " << rms_exact << "\n";

    /** write time history **/
    int width = std::to_string(nMics).size();
    std::string num = std::to_string(iMic);
    std::string padded = std::string(width - num.length(),'0') + num;
    std::string fname = "signal_" + padded + ".dat";
    std::ofstream signal(fname);
    signal << "# t[s] p_FWH[Pa] p_exact[Pa]\n";

    for (size_t i = obs.trim_start; i < obs.trim_end; i++) {
      double t = obs.t0 + i * obs.dt;
      double p_fwh = obs.signal[i];
      double p_exact = provider.exact_presure_at(obs.position[0],
                                                 obs.position[1],
                                                 obs.position[2],
                                                 t);
      signal << t << " " << p_fwh << " " << p_exact << "\n";
    }
    signal.close();

  }

  out.close();

}