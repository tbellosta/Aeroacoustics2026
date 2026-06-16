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
#ifndef FWHSOLVER_H
#define FWHSOLVER_H

#include <vector>
#include <utility>

#include "geometry.h"
#include "meshStructure.h"


class FWHSolver {
 public:

   double c0   = 340.0;
   double p0   = 101325.0;
   double rho0 = 1.2;

   Vect3 M0_flow;
   bool permeable = false;

   std::vector<Observer> observers;


   void initialize_observers(const SurfaceData& firstSnapshot,
                             double tSource_first,
                             double tSource_end);

   void process(const SurfaceData& prev,
                const SurfaceData& curr,
                const SurfaceData& next,
                double tCurr,
                double dt);

private:

    std::pair<double,double> compute_FWH_F1A(const Node& n_prev,
                           const Node& n_curr,
                           const Node& n_next,
                           const Vect3& obs,
                           double dt) const;

};

#endif //FWHSOLVER_H
