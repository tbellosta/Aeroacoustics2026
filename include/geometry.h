//============================================================
//
//      Type:        FWH include file
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

#ifndef GEOMETRY_H
#define GEOMETRY_H

class Vect3 {
private:
  double x[3]{};

public:

  /** contructors are special member functions
   * that get called at object creation **/
  Vect3() {x[0] = 0; x[1] = 0; x[2] = 0;}
  Vect3(double _x, double _y, double _z) {x[0] = _x; x[1] = _y; x[2] = _z;}

  /** The desctructor gets called when an object is deleted/goes
   * out of scope **/
  ~Vect3() = default;

  double& operator[](int i) {
    return x[i];
  }
  const double& operator[](const int i) const {
    return x[i];
  }

  Vect3 operator+(const Vect3& v) const {
    Vect3 out;
    for (int i = 0; i < 3; i++) {
      out[i] = x[i] + v[i];
    }
    return out;
  }


};




#endif