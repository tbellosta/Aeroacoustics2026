//============================================================
//
//      Type:        FWH source file
//
//      Author:      Tommaso Bellosta on 14/04/26.
//                   Dipartimento di Scienze e Tecnologie Aerospaziali
//                   Politecnico di Milano
//                   Via La Masa 34, 20156 Milano, ITALY
//                   e-mail: tommaso.bellosta@polimi.it
//
//      Copyright:   2026, authors above and the FWH contributors.
//                   This software is distributed under the MIT license, see LICENSE.txt
//
//============================================================
#include <cmath>

#include "geometry.h"

double dot(const Vect3& v1, const Vect3& v2) {
    double out = 0;
    for (int i = 0; i < 3; i++) out += v1[i] * v2[i];
    return out;
}


Vect3 cross(const Vect3& v1, const Vect3& v2) {
    Vect3 out;
    out[0] = v1[1]*v2[2] - v1[2]*v2[1];
    out[1] = v1[2]*v2[0] - v1[0]*v2[2];
    out[2] = v1[0]*v2[1] - v1[1]*v2[0];
    return out;
}

double norm(const Vect3& v) {
    double out = 0;
    for (int i = 0; i < 3; i++) out += v[i] * v[i];
    return std::sqrt(out);
}
