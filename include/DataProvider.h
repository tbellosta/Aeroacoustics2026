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

#include <string>
#include <sstream>

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
                             size_t nsteps_, double r_ = 1, size_t nrefine_ = 2, double M0_ = 0.0,
                             double c0_ = 340.0, double rho0_ = 1.2, double p0_ = 101325.0);

  bool has_next() const;
  void load_next(SurfaceData& sData, double& time);


  double exact_presure_at(double x0, double y0, double z0, double t0) const;

};


class AnaliticalDipoleProvider {
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

    AnaliticalDipoleProvider(double ampl_, double freq_, double dt_,
                               size_t nsteps_, double r_ = 1, size_t nrefine_ = 2, double M0_ = 0.0,
                               double c0_ = 340.0, double rho0_ = 1.2, double p0_ = 101325.0);

    bool has_next() const;
    void load_next(SurfaceData& sData, double& time);


    double exact_presure_at(double x0, double y0, double z0, double t0) const;

};


class MovingMonopoleProvider {
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
    double Ms; // this is the source speed in the x-direction

    double R;
    size_t n_refine;

    MovingMonopoleProvider(double ampl_, double freq_, double Ms_, double dt_,
                               size_t nsteps_, double r_ = 1, size_t nrefine_ = 2,
                               double c0_ = 340.0, double rho0_ = 1.2, double p0_ = 101325.0);

    bool has_next() const;
    void load_next(SurfaceData& sData, double& time);



};

class VTUProvider {
private:
    double dt;
public:
    std::vector<std::string> filenames;
    size_t step = 0;

    VTUProvider(double dt_, const std::vector<std::string>& filenames_);

    bool has_next() const;
    void load_next(SurfaceData& data, double& time);
};

// ====================================================================
//  PanelCSVProvider
//
//  Reads aeroacoustic surface data exported by a vortex-lattice (VL)
//  panel code. One file = one time snapshot.
//
//  The VL surface is a thin mean sheet: pressure (the load dp across
//  the sheet) lives on the mean panel, but a thin sheet carries no
//  thickness. To recover thickness we reconstruct the actual *closed*
//  wetted surface from the upper/lower sub-panel geometry the exporter
//  provides, and emit TWO nodes per element (upper face + lower face).
//
//  On that closed surface (permeable = false):
//    - thickness  -> from the closed geometry + flight-frame convection
//    - loading    -> surface-pressure dipole F = (p - p0) n, with the
//                    mean-panel load dp split symmetrically onto the
//                    two faces:  p_u = p0 - dp/2,  p_l = p0 + dp/2.
//
//  dS is set directly (n * area); there is no connectivity, so
//  compute_dual_areas() is intentionally NOT called.
//
//  Node ordering (upper, lower) per element, in file order, is the
//  same in every snapshot, so the solver's per-node prev/curr/next
//  finite differences stay consistent.
//
//  File format (comma separated, Fortran "E+000" exponents):
//    # banner
//    # Element,Time,Pressure,Density,Sound speed,Dynamic viscosity,Flow velocity
//    nElem, time, p_ref, rho_ref, c_ref, mu, Ux, Uy, Uz
//    # cx,cy,cz,nx,ny,nz,area,rho,pressure,rhoux,rhouy,rhouz,svx,svy,svz,
//    #   cx_u,cy_u,cz_u,cx_l,cy_l,cz_l,nx_u,ny_u,nz_u,nx_l,ny_l,nz_l,
//    #   area_u,area_l,svx_u,svy_u,svz_u,svx_l,svy_l,svz_l   (35 fields)
//    <element rows...>
// ====================================================================
class PanelCSVProvider {
public:
    std::vector<std::string> filenames;
    int current = 0;
    double dt;

    // Reference state, read from each file's metadata row.
    // Pull these onto the solver after the first load_next().
    double c0   = 340.0;
    double rho0 = 1.0;
    double p0   = 0.0;
    Vect3   U_inf;            // freestream velocity vector

    PanelCSVProvider(double dt_, const std::vector<std::string>& filenames_) : dt(dt_), filenames(filenames_) {};

    bool has_next() const { return current < (int)filenames.size(); }
    void load_next(SurfaceData& data, double& time);

private:

    static std::vector<double> parse_csv(const std::string& line) {
        std::vector<double> vals;
        std::stringstream ss(line);
        std::string tok;
        while (std::getline(ss, tok, ',')) {
            size_t a = tok.find_first_not_of(" \t\r\n");
            if (a == std::string::npos) continue;
            vals.push_back(std::stod(tok.substr(a)));
        }
        return vals;
    }
};

#endif //DATAPROVIDER_H
