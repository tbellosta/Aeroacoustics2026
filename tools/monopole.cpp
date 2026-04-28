//============================================================
//
//      Type:        FWH include file
//
//      Author:      Tommaso Bellosta on 27/04/26.
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

int main() {

  /** monopole characteristics (and FWH surface) **/
  double A = 1.0;
  double frq = 10.0;
  double r = 1.0;
  size_t nRefine = 2;

  double src_dt = 1e-4;
  size_t nsteps = 12000;

  /** This is the data provider **/
  AnaliticalMonopoleProvider provider(A,frq,src_dt,nsteps,r,nRefine);

  /** This stores the 3 snapshots and cycles them **/
  SnapshotBuffer buffer;
  buffer.initialize(provider);

  /** the solver object **/
  FWHSolver solver;
  solver.permeable = true; // tells the solver we want the permeable surface formulation
  solver.M0_flow  = Vect3(0.0, 0.0, 0.0); // no mean convection

  /** define a microphone array **/
  double R = 50.0;
  size_t nMics = 72;

  /** initializes the mics positions and store the observer into the solver object **/
  for (size_t m = 0; m < nMics; m++) {
    double theta = 2.0 * M_PI * m / nMics;
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

}