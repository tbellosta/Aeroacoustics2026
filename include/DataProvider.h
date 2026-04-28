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
#ifndef DATAPROVIDER_H
#define DATAPROVIDER_H

#include "meshStructure.h"

/** Classes that feed data into out SnapshotBuffer objects **/

class AnaliticalMonopoleProvider {
private:

  connectivity sphere_elements;
  std::vector<Vect3> sphere_nodes;

public:

  double A;
  double f;
  double dt;
  size_t nsteps;
  size_t step;

  double c0, rho0, p0;
  double M0; // this is the mean-flow mach number in x-direction;

  double R;
  size_t n_refine;

  AnaliticalMonopoleProvider(double ampl_, double freq_, double dt_,
                             size_t nsteps_, double r_ = 1, size_t nrefine_ = 2,
                             double c0_ = 340.0, double rho0_ = 1.2, double p0_ = 101325.0,
                             double M0_ = 0.0);

  bool has_next() const;
  void load_next(SurfaceData& sData, double& time);


  double exact_presure_at(double x0, double y0, double z0, double t0) const;

};

#endif //DATAPROVIDER_H
