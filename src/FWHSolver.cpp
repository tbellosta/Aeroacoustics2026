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
#include "FWHSolver.h"


void FWHSolver::initialize_observers(const SurfaceData& firstSnapshot,
                             double tSource_first,
                             double tSource_end) {

  /** For each observer do: **/
  for (auto& observer : observers) {

    /** compute the max and min delay between
     *  the emitted sound at the sourface and the time
        at which the observer receives it
     **/
    double min_delay = 1e200;
    double max_delay = -1e200;

    for (const auto& node : firstSnapshot.nodes) {

      Vect3 Delta = observer.position - node.x;
      double beta2 = 1.0 - dot(M0_flow,M0_flow);
      double Mr = dot(M0_flow,Delta);
      double D2 = dot(Delta,Delta);

      // double delay = norm(R) / c0;
      double delay = (std::sqrt(D2*beta2 + Mr*Mr) - Mr) / (beta2 * c0);

      min_delay = std::min(min_delay, delay);
      max_delay = std::max(max_delay, delay);

    }

    /** time at which the observer starts receiving
      * the signal **/
    double t0   = tSource_first + min_delay;

    /** time at which the observer stops receiving
      * the signal. **/
    double tEnd = tSource_end   + max_delay;

    /** compute the number af time steps given t0,tEnd,dt **/
    int Nt = std::ceil((tEnd - t0) / observer.dt) + 1;
    /** since (tEnd-t0)/dt is not an integrar number, we
      * must choose between matching the original dt or tEnd.
      * We match dt and stretch tEnd. **/
    tEnd = t0 + (Nt-1) * observer.dt;

    observer.initialize(t0, observer.dt, Nt);

    /** We can also take care of the "transient"
      * initial and final samples of the observer history.
      * Depending on the compactness of the source and resolution
      * of the observer signal, the initial/final time samples may only contain
      * signal contributions from part of the emission surface.
      * We chose to indicate inside the Observer class
      * the indexes in the signal vector corresponging to the
      * beginning and and of the "good" part of the signal.
      * This requires defining two new int in Observer. **/

    double t_clean_start = tSource_first + max_delay;
    double t_clean_end = tSource_end + min_delay;

    /** @TODO need to define those first **/
    observer.trim_start = std::ceil((t_clean_start - t0) / observer.dt);
    observer.trim_end   = std::floor((t_clean_end - t0) / observer.dt);


  }

 }


void FWHSolver::process(const SurfaceData& prev,
                const SurfaceData& curr,
                const SurfaceData& next,
                double tCurr,
                double dt) {

  /** Main loop in the FWH computation **/
  size_t nNodes = curr.nodes.size();
  for (size_t iNode = 0; iNode < nNodes; iNode++) {

    const Node& n_prev = prev.nodes[iNode];
    const Node& n_curr = curr.nodes[iNode];
    const Node& n_next = next.nodes[iNode];

    for (auto& obs : observers) {

      /** compute the FWH integral and then assign the computed value
       * at the correct time in observer time history.
       * FWH integral computation is computed at source/retarded time
       * whereas we want the signal at the observer time. **/

      double pprime = compute_FWH_F1A(n_prev,n_curr,n_next,
                                      obs.position,dt);


      /** Assign the computed pressure at the observer time **/
      Vect3 Delta = obs.position - n_curr.x;
      double beta2 = 1.0 - dot(M0_flow,M0_flow);
      double Mr = dot(M0_flow,Delta);
      double D2 = dot(Delta,Delta);

      // double delay = norm(R) / c0;
      double delay = (std::sqrt(D2*beta2 + Mr*Mr) - Mr) / (beta2 * c0);

      double tArrival = tCurr + delay;
      obs.add(pprime*dt/obs.dt, tArrival);

    }



  }


}


double FWHSolver::compute_FWH_F1A(const Node& n_prev,
                           const Node& n_curr,
                           const Node& n_next,
                           const Vect3& obs,
                           double dt) const {

  /** Computes the FWH contribution from a single surface location to a
   * single observer location **/

  const Vect3 Uinf = M0_flow * c0;

  /** compute geometry **/
  const Vect3 Delta = obs - n_curr.x;
  double beta2 = 1.0 - dot(M0_flow,M0_flow);
  double D2 = dot(Delta,Delta);
  double M0dotDelta = dot(M0_flow,Delta);
  const double R = (std::sqrt(D2*beta2 + M0dotDelta*M0dotDelta) - M0dotDelta) / (beta2);

  Vect3 r_hat = (Delta - M0_flow * R);
  r_hat = r_hat * (1.0 / R);

  const double S = norm(n_curr.dS);
  const Vect3 n = n_curr.dS * (1.0 / S);

  const double S_prev = norm(n_prev.dS);
  const Vect3 n_m = n_prev.dS * (1.0 / S_prev);

  const double S_next = norm(n_next.dS);
  const Vect3 n_p = n_next.dS * (1.0 / S_next);

  /** Compute the velocity of the moving surface **/
  const Vect3 v_s = (n_next.x - n_prev.x) * (0.5 / dt) - Uinf;
  const Vect3 v_s_next = (n_next.x - n_curr.x) * (1.0 / dt) - Uinf;
  const Vect3 v_s_prev = (n_curr.x - n_prev.x) * (1.0 / dt) - Uinf;


  /** Compute projections in the direction of the observer **/
  const double Mr = dot(v_s, r_hat) / c0;
  const Vect3 Mdot = (v_s_next - v_s_prev) * (1.0 / (dt*c0));
  const double Mdotr = dot(Mdot,r_hat);
  const double Mr1 = 1.0 - Mr;

  /** Compute the Thickness quantities **/

  Vect3 Q, Q_prev, Q_next;

  if (!permeable) {
    Q      = v_s;
    Q_next = v_s_next;
    Q_prev = v_s_prev;
  } else {
    const Vect3 uRel_curr = n_curr.u - Uinf;
    const Vect3 uRel_prev = n_prev.u - Uinf;
    const Vect3 uRel_next = n_next.u - Uinf;
    Q = uRel_curr +  (uRel_curr - v_s) * (n_curr.rho/rho0 - 1.0);
    /** v_s_* and n_*.u(rho) are not defined at the same time instant.
     * @TODO interpolate n_*.solution at the staggered location **/
    Q_prev = uRel_prev +  (uRel_prev - v_s_prev) * (n_prev.rho/rho0 - 1.0);
    Q_next = uRel_next +  (uRel_next - v_s_next) * (n_next.rho/rho0 - 1.0);
  }

  double Qn = dot(Q,n);
  /** @TODO if solution is interpolated at the staggered time,
   * then divide by dt and not 2*dt **/
  Vect3 Qdot = (Q_next - Q_prev) * (0.5 / dt);
  double Qdotn = dot(Qdot,n);

  /** Compute the loading quantities.
   * Lij = pij + rho*ui*(uj - vj)
   * Fi = Lij nj **/

  Vect3 F, F_prev, F_next;

  if (!permeable) {
    F      = n   * (n_curr.p - p0);
    F_prev = n_m * (n_prev.p - p0);
    F_next = n_p * (n_next.p - p0);
  } else {

    const Vect3 uRel_curr = n_curr.u - Uinf;
    const Vect3 uRel_prev = n_prev.u - Uinf;
    const Vect3 uRel_next = n_next.u - Uinf;

    F = n * (n_curr.p - p0) + uRel_curr * dot(uRel_curr - v_s,n) * n_curr.rho;
    /** same comment as for the thickness vector regarding the staggering of the
     * surface velocity **/
    F_prev = n_m * (n_prev.p - p0) + uRel_prev * dot(uRel_prev - v_s_prev,n_m) * n_prev.rho;
    F_next = n_p * (n_next.p - p0) + uRel_next * dot(uRel_next - v_s_next,n_p) * n_next.rho;

  }

  /** normal time derivative **/
  Vect3 ndot = (n_p - n_m) * (0.5 / dt);
  double ndotQ = dot(Q,ndot);

  double Fr = dot(F,r_hat);

  Vect3 Fdot = (F_next - F_prev) * (0.5 / dt);
  double Fdotr = dot(Fdot,r_hat);
  double FM = dot(F,v_s) / c0;

  double K = R*Mdotr + Mr*c0 - dot(v_s,v_s) / c0;

  /** assemble the integrals **/
  double denom = std::max(1e-12,Mr1);


  const double T1 = rho0 * (Qdotn + ndotQ) / (R*std::pow(denom,2));
  const double T2 = rho0 * Qn * K / (R*R*std::pow(denom,3));

  const double L1 = Fdotr / (R*std::pow(denom,2)*c0);
  const double L2 = (Fr - FM) / (R*R*denom*denom);
  const double L3 = Fr * K / (R*R * std::pow(denom,3)*c0);

  return (T1+T2+L1+L2+L3) * S / (4*M_PI);


}
