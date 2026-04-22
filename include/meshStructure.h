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

  size_t trim_start, trim_end;

  void initialize(double t0_, double dt_, size_t Nt_);
  void add(double value, double t_arrival);

};

// DataProvider
// DataProvider::load_next(SurfaceData sData, double time)
// bool DataProvider::has_next()

class SnapshotBuffer {
public:
  SurfaceData prev, curr, next;
  double t_prev, t_curr, t_next;

  template<typename Provider>
  void initialize(Provider& p) {
    p.load_next(prev, t_prev);
    p.load_next(curr, t_curr);
    p.load_next(next, t_next);
  }

  template<typename Provider>
  bool advance(Provider& p) {
    prev = curr;
    curr = next;

    if (!p.has_next()) return false;

    t_prev = t_curr;
    t_curr = t_next;

    p.load_next(next, t_next);

    return true;
  }


};




#endif //MESHSTRUCTURE_H
