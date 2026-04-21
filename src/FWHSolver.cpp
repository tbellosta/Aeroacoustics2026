//============================================================
//
//      Type:        FWH include file
//
//      Author:      Tommaso Bellosta on 20/04/26.
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


void FWHSolver::initialize_observers(const SurfaceData& firstSnapshot,
                             double tSource_first,
                             double tSource_end) {

  /** For each observer do: **/
  for (auto& observer : observers) {

    /** compute the max and min delay between
     *  the emitted sound at the sourface and the time
        at which the observer receives it
     **/
    double min_delay = 1e200;
    double max_delay = -1e200;

    for (auto& node : firstSnapshot.nodes) {

      Vect3 R = observer.position - node.x;
      double delay = norm(R) / c0;

      min_delay = std::min(min_delay, delay);
      max_delay = std::max(max_delay, delay);

    }

    /** time at which the observer starts receiving
      * the signal **/
    double t0   = tSource_first + min_delay;

    /** time at which the observer stops receiving
      * the signal. **/
    double tEnd = tSource_end   + max_delay;

    /** compute the number af time steps given t0,tEnd,dt **/
    int Nt = std::ceil((tEnd - t0) / observer.dt) + 1;
    /** since (tEnd-t0)/dt is not an integrar number, we
      * must choose between matching the original dt or tEnd.
      * We match dt and stretch tEnd. **/
    tEnd = t0 + (Nt-1) * observer.dt;

    observer.initialize(t0, tEnd, Nt);

    /** We can also take care of the "transient"
      * initial and final samples of the observer history.
      * Depending on the compactness of the source and resolution
      * of the observer signal, the initial/final time samples may only contain
      * signal contributions from part of the emission surface.
      * We chose to indicate inside the Observer class
      * the indexes in the signal vector corresponging to the
      * beginning and and of the "good" part of the signal.
      * This requires defining two new int in Observer. **/

    double t_clean_start = tSource_first + max_delay;
    double t_clean_end = tSource_end + min_delay;

    /** @TODO need to define those first **/
//    observer.trim_start = std::ceil((t_clean_start - t0) / observer.dt);
//    observer.trim_end   = std::floor((t_clean_end - t0) / observer.dt);


  }

 }