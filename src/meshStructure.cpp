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
#include "meshStructure.h"

void SurfaceData::compute_dual_geometry() {

  const size_t nElements = elements.size();

  for (const auto& element : elements) {
    const size_t i = element[0];
    const size_t j = element[1];
    const size_t k = element[2];

    const Vect3& xi = nodes[i].x;
    const Vect3& xj = nodes[j].x;
    const Vect3& xk = nodes[k].x;

    Vect3 v1,v2;
    v1 = xj - xi;
    v2 = xk - xi;

    /** This is an approximation (works for regular trianguler grid).
     * @TODO Consider implementing proper median dual geometry. **/
    Vect3 normal = cross(v1,v2);
    normal = normal * (0.5 / 3.0);

    nodes[i].dS += normal;
    nodes[j].dS += normal;
    nodes[k].dS += normal;

  }

}



void Observer::add(const double value, const double t_arrival) {

  const double idx = (t_arrival - t0) / dt;
  const int i = static_cast<int>(idx);

  if (i < 0 || i+1 >= signal.size()) return;

  const double alpha = idx - i;

  signal[i]   += value * (1.0 - alpha);
  signal[i+1] += value * alpha;
}

void Observer::initialize(double t0_, double dt_, size_t Nt_) {
  t0 = t0_;
  dt = dt_;

  /** All three approaches are functionally
   * equivalent **/

  // std::vector<double> zero(Nt_,0.0);
  // signal = zero;

  signal = std::vector<double>(Nt_,0.0);

  // signal.resize(Nt_);
  // for (auto& s : signal) s = 0.0;

}