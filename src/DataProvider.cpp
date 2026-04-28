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

#include <complex>

#include "DataProvider.h"

AnaliticalMonopoleProvider::AnaliticalMonopoleProvider(double ampl_, double freq_, double dt_,
                             size_t nsteps_, double r_, size_t nrefine_,
                             double c0_, double rho0_, double p0_,
                             double M0_) {

  A = ampl_;
  f = freq_;
  dt = dt_;
  nsteps = nsteps_;
  R = r_;
  n_refine = nrefine_;
  c0 = c0_;
  rho0 = rho0_;
  p0 = p0_;
  M0 = M0_;


  step = 0;

  generate_icosphere(sphere_nodes, sphere_elements, n_refine,R);

 }


bool AnaliticalMonopoleProvider::has_next() const {
    return step < nsteps;
}


void AnaliticalMonopoleProvider::load_next(SurfaceData& sData, double& time) {

    time = step * dt;

    /** set the geometry into out output data (sData) **/
    size_t nNodes = sphere_nodes.size();
    sData.nodes.resize(nNodes);

    sData.elements = sphere_elements;

    for (int iNode = 0; iNode < nNodes; iNode++)
        sData.nodes[iNode].x = sphere_nodes[iNode];

    /** Set the analytical monopole solution at the reference surface **/
    const double omega = 2 * M_PI * f;
    const double beta2 = 1 - M0*M0;
    const double U0 = M0 * c0;

    /** imaginary unit **/
    const std::complex<double> j(0.0,1.0);

    for (size_t iNode = 0; iNode < nNodes; iNode++) {

        double x = sphere_nodes[iNode][0];
        double y = sphere_nodes[iNode][1];
        double z = sphere_nodes[iNode][2];

        double Rstar  = std::sqrt(x*x + beta2*y*y + beta2*z*z);
        double Rtilde = (Rstar - M0 * x) / beta2;
        double psi    = time - Rtilde / c0;

        std::complex<double> phi = A / (4.0 * M_PI * Rstar)
                                 * std::exp(j * omega * psi);

        std::complex<double> dphi_dt = j * omega * phi;

        double dRs_dx = x / Rstar;
        double dRs_dy = beta2 * y / Rstar;
        double dRs_dz = beta2 * z / Rstar;

        double dRt_dx = (x / Rstar - M0) / beta2;
        double dRt_dy = y / Rstar;
        double dRt_dz = z / Rstar;

        std::complex<double> dphi_dx = phi * (-dRs_dx/Rstar - j*omega/c0 * dRt_dx);
        std::complex<double> dphi_dy = phi * (-dRs_dy/Rstar - j*omega/c0 * dRt_dy);
        std::complex<double> dphi_dz = phi * (-dRs_dz/Rstar - j*omega/c0 * dRt_dz);

        double p_prime = std::real(-rho0 * (dphi_dt + U0 * dphi_dx));

        sData.nodes[iNode].p = p0 + p_prime;
        sData.nodes[iNode].u = Vect3(
            U0 + std::real(dphi_dx),
                 std::real(dphi_dy),
                 std::real(dphi_dz)
        );
        sData.nodes[iNode].rho = rho0 + p_prime / (c0 * c0);

    }

    sData.compute_dual_geometry();

    step++;

}


double AnaliticalMonopoleProvider::exact_presure_at(double x, double y, double z, double t0) const  {

    /** Set the analytical monopole solution at the reference surface **/
    const double omega = 2 * M_PI * f;
    const double beta2 = 1 - M0*M0;
    const double U0 = M0 * c0;

    /** imaginary unit **/
    const std::complex<double> j(0.0,1.0);

    double Rstar  = std::sqrt(x*x + beta2*y*y + beta2*z*z);
    double Rtilde = (Rstar - M0 * x) / beta2;
    double psi    = t0 - Rtilde / c0;

    std::complex<double> phi = A / (4.0 * M_PI * Rstar)
                             * std::exp(j * omega * psi);

    std::complex<double> dphi_dt = j * omega * phi;

    double dRs_dx = x / Rstar;
    double dRs_dy = beta2 * y / Rstar;
    double dRs_dz = beta2 * z / Rstar;

    double dRt_dx = (x / Rstar - M0) / beta2;
    double dRt_dy = y / Rstar;
    double dRt_dz = z / Rstar;

    std::complex<double> dphi_dx = phi * (-dRs_dx/Rstar - j*omega/c0 * dRt_dx);
    std::complex<double> dphi_dy = phi * (-dRs_dy/Rstar - j*omega/c0 * dRt_dy);
    std::complex<double> dphi_dz = phi * (-dRs_dz/Rstar - j*omega/c0 * dRt_dz);

    return std::real(-rho0 * (dphi_dt + U0 * dphi_dx));

}
