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
    const size_t k = element[3];

    const Vect3& xi = nodes[i].x;
    const Vect3& xj = nodes[j].x;
    const Vect3& xk = nodes[k].x;

    Vect3 v1,v2;
    v1 = xj - xi;
    v2 = xk - xi;

    /** This is an approximation (works for regular trianguler grid).
     * @TODO Consider implementing proper median dual geometry. **/
    Vect3 normal = cross(v1,v2);
    normal = normal * (0.5 / 3);

    nodes[i].dS += normal;
    nodes[j].dS += normal;
    nodes[k].dS += normal;

  }

}