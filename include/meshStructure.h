//============================================================
//
//      Type:        FWH include file
//
//      Author:      Tommaso Bellosta on 15/04/26.
//                   Dipartimento di Scienze e Tecnologie Aerospaziali
//                   Politecnico di Milano
//                   Via La Masa 34, 20156 Milano, ITALY
//                   e-mail: tommaso.bellosta@polimi.it
//
//      Copyright:   2026, authors above and the FWH contributors.
//                   This software is distributed under the MIT license, see LICENSE.txt
//
//============================================================
#ifndef MESHSTRUCTURE_H
#define MESHSTRUCTURE_H

#include <vector>

#include "geometry.h"

typedef std::vector<std::vector<size_t>> connectivity;

struct Node {
  Vect3 x; // position
  Vect3 dS; // dual area vector

  double p;
  Vect3 u;
  double rho;
};

class SurfaceData {
public:
  connectivity elements;
  std::vector<Node> nodes;

  void compute_dual_geometry();
};


class Observer {
public:
  Vect3 position;

  double t0, dt;
  std::vector<double> signal;

  void initialize(double t0_, double dt_, size_t Nt_);
  void add(double value, double t_arrival);

};




#endif //MESHSTRUCTURE_H
