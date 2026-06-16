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

#include <mpi.h>

#include <vtkXMLUnstructuredGridReader.h>
#include <vtkUnstructuredGrid.h>
#include <vtkPointData.h>

#include <complex>

#include "DataProvider.h"

AnaliticalMonopoleProvider::AnaliticalMonopoleProvider(double ampl_, double freq_, double dt_,
                             size_t nsteps_, double r_, size_t nrefine_, double M0_,
                             double c0_, double rho0_, double p0_) {

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


AnaliticalDipoleProvider::AnaliticalDipoleProvider(double ampl_, double freq_, double dt_,
                               size_t nsteps_, double r_, size_t nrefine_, double M0_,
                               double c0_, double rho0_, double p0_) {

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


bool AnaliticalDipoleProvider::has_next() const {
    return step < nsteps;
}


void AnaliticalDipoleProvider::load_next(SurfaceData& sData, double& time) {

    time = step * dt;

    /** set the geometry into out output data (sData) **/
    size_t nNodes = sphere_nodes.size();
    sData.nodes.resize(nNodes);

    sData.elements = sphere_elements;

    for (int iNode = 0; iNode < nNodes; iNode++)
        sData.nodes[iNode].x = sphere_nodes[iNode];

    /** Set the analytical monopole solution at the reference surface **/
    const double omega = 2 * M_PI * f;
    const double k = omega / c0; // wavenumber
    const double beta2 = 1 - M0*M0;
    const double U0 = M0 * c0;

    /** imaginary unit **/
    const std::complex<double> j(0.0,1.0);

    for (size_t iNode = 0; iNode < nNodes; iNode++) {

        double x = sphere_nodes[iNode][0];
        double y = sphere_nodes[iNode][1];
        double z = sphere_nodes[iNode][2];

        double Rstar = std::sqrt(x*x + beta2*y*y + beta2*z*z);
        double Rs2   = Rstar * Rstar;
        double Rs3   = Rs2 * Rstar;
        double Rs4   = Rs2 * Rs2;

        double Rtilde = (Rstar - M0*x) / beta2;
        double psi    = time - Rtilde / c0;

        // --- Monopole potential ---
        std::complex<double> phi = A / (4.0*M_PI*Rstar)
                                 * std::exp(j * omega * psi);

        // --- B_i factors:  d phi_mono / d x_i = phi * B_i ---
        std::complex<double> B_x = -x/Rs2 - j*k*(x/Rstar - M0)/beta2;
        std::complex<double> B_y = -beta2*y/Rs2 - j*k*y/Rstar;
        std::complex<double> B_z = -beta2*z/Rs2 - j*k*z/Rstar;

        // --- dB_i/dy  (using  dR*/dy = beta^2 y / R*) ---
        //
        //  dB_x/dy = 2 beta^2 x y / R*^4  + i k x y / R*^3
        //
        //  dB_y/dy = -beta^2/R*^2 + 2 beta^4 y^2/R*^4
        //          - i k / R*     + i k beta^2 y^2/R*^3
        //
        //  dB_z/dy = 2 beta^4 y z / R*^4  + i k beta^2 y z / R*^3

        std::complex<double> dBx_dy = 2.0*beta2*x*y/Rs4
                                    + j*k*x*y/Rs3;

        std::complex<double> dBy_dy = -beta2/Rs2 + 2.0*beta2*beta2*y*y/Rs4
                                    - j*k/Rstar  + j*k*beta2*y*y/Rs3;

        std::complex<double> dBz_dy = 2.0*beta2*beta2*y*z/Rs4
                                    + j*k*beta2*y*z/Rs3;

        // --- Dipole potential and its derivatives ---
        std::complex<double> phi_dip    = phi * B_y;
        std::complex<double> dphi_dip_dt = j * omega * phi_dip;

        std::complex<double> dphi_dip_dx = phi * (B_y * B_x + dBx_dy);
        std::complex<double> dphi_dip_dy = phi * (B_y * B_y + dBy_dy);
        std::complex<double> dphi_dip_dz = phi * (B_y * B_z + dBz_dy);

        // --- Physical fields ---
        double p_prime = std::real(-rho0 * (dphi_dip_dt + U0 * dphi_dip_dx));

        sData.nodes[iNode].p   = p0 + p_prime;
        sData.nodes[iNode].u   = Vect3(U0 + std::real(dphi_dip_dx),
                                      std::real(dphi_dip_dy),
                                      std::real(dphi_dip_dz));
        sData.nodes[iNode].rho = rho0 + p_prime / (c0 * c0);


    }

    sData.compute_dual_geometry();

    step++;

}


double AnaliticalDipoleProvider::exact_presure_at(double x, double y, double z, double t) const {

    /** Set the analytical monopole solution at the reference surface **/
    const double omega = 2 * M_PI * f;
    const double k = omega / c0; // wavenumber
    const double beta2 = 1 - M0*M0;
    const double U0 = M0 * c0;

    /** imaginary unit **/
    const std::complex<double> j(0.0,1.0);

    double Rstar = std::sqrt(x*x + beta2*y*y + beta2*z*z);
    double Rs2   = Rstar * Rstar;
    double Rs3   = Rs2 * Rstar;
    double Rs4   = Rs2 * Rs2;

    double Rtilde = (Rstar - M0*x) / beta2;
    double psi    = t - Rtilde / c0;

    // --- Monopole potential ---
    std::complex<double> phi = A / (4.0*M_PI*Rstar)
                             * std::exp(j * omega * psi);

    // --- B_i factors:  d phi_mono / d x_i = phi * B_i ---
    std::complex<double> B_x = -x/Rs2 - j*k*(x/Rstar - M0)/beta2;
    std::complex<double> B_y = -beta2*y/Rs2 - j*k*y/Rstar;
    std::complex<double> B_z = -beta2*z/Rs2 - j*k*z/Rstar;

    // --- dB_i/dy  (using  dR*/dy = beta^2 y / R*) ---
    //
    //  dB_x/dy = 2 beta^2 x y / R*^4  + i k x y / R*^3
    //
    //  dB_y/dy = -beta^2/R*^2 + 2 beta^4 y^2/R*^4
    //          - i k / R*     + i k beta^2 y^2/R*^3
    //
    //  dB_z/dy = 2 beta^4 y z / R*^4  + i k beta^2 y z / R*^3

    std::complex<double> dBx_dy = 2.0*beta2*x*y/Rs4
                                + j*k*x*y/Rs3;

    std::complex<double> dBy_dy = -beta2/Rs2 + 2.0*beta2*beta2*y*y/Rs4
                                - j*k/Rstar  + j*k*beta2*y*y/Rs3;

    std::complex<double> dBz_dy = 2.0*beta2*beta2*y*z/Rs4
                                + j*k*beta2*y*z/Rs3;

    // --- Dipole potential and its derivatives ---
    std::complex<double> phi_dip    = phi * B_y;
    std::complex<double> dphi_dip_dt = j * omega * phi_dip;

    std::complex<double> dphi_dip_dx = phi * (B_y * B_x + dBx_dy);
    std::complex<double> dphi_dip_dy = phi * (B_y * B_y + dBy_dy);
    std::complex<double> dphi_dip_dz = phi * (B_y * B_z + dBz_dy);

    // --- Physical fields ---
    return std::real(-rho0 * (dphi_dip_dt + U0 * dphi_dip_dx));

}





MovingMonopoleProvider::MovingMonopoleProvider(double ampl_, double freq_, double Ms_, double dt_,
                               size_t nsteps_, double r_, size_t nrefine_,
                               double c0_, double rho0_, double p0_) {

    A = ampl_;
    f = freq_;
    dt = dt_;
    nsteps = nsteps_;
    R = r_;
    n_refine = nrefine_;
    c0 = c0_;
    rho0 = rho0_;
    p0 = p0_;
    Ms = Ms_;


    step = 0;

    generate_icosphere(sphere_nodes, sphere_elements, n_refine,R);

}


bool MovingMonopoleProvider::has_next() const {
    return step < nsteps;
}


void MovingMonopoleProvider::load_next(SurfaceData& sData, double& time) {

    time = step * dt;

    /** set the geometry into out output data (sData) **/
    size_t nNodes = sphere_nodes.size();
    sData.nodes.resize(nNodes);

    sData.elements = sphere_elements;

    for (int iNode = 0; iNode < nNodes; iNode++)
        sData.nodes[iNode].x = sphere_nodes[iNode];

    /** Set the analytical monopole solution at the reference surface **/
    const double omega = 2 * M_PI * f;
    const double M0 = 0.0;
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

    /** rigidly diplace the sphere points according to xp = x0 + c0 * (M0, 0, 0) * (t - t0) **/

    for (size_t iNode = 0; iNode < nNodes; iNode++) {
        auto& xp = sData.nodes[iNode].x;
        const auto& x0 = sphere_nodes[iNode];
        const Vect3 vel(c0*Ms,0,0);

        xp = x0 + vel * time;

    }

    sData.compute_dual_geometry();

    step++;

}


VTUProvider::VTUProvider(double dt_, const std::vector<std::string> &filenames_) {

    dt = dt_;
    filenames = filenames_;

}

bool VTUProvider::has_next() const {
    return step < filenames.size();
}

void VTUProvider::load_next(SurfaceData& data, double& time) {

    std::string file = filenames[step];
    time = step * dt;

    vtkNew<vtkXMLUnstructuredGridReader> reader;
    reader->SetFileName(file.c_str());
    reader->Update();

    auto grid = reader->GetOutput();

    auto N = grid->GetNumberOfPoints();
    data.nodes.resize(N);

    for (int i = 0; i < N; ++i) {
        double x[3];
        grid->GetPoint(i, x);
        for (int d = 0; d < 3; ++d)
            data.nodes[i].x[d] = x[d];
    }

    auto p_array = grid->GetPointData()->GetArray("Pressure");
    for (int i = 0; i < N; ++i)
        data.nodes[i].p = p_array->GetTuple1(i);

    auto rho_array = grid->GetPointData()->GetArray("Density");
    for (int i = 0; i < N; ++i)
        data.nodes[i].rho = rho_array->GetTuple1(i);

    auto u_array = grid->GetPointData()->GetArray("Velocity");
    for (int i = 0; i < N; ++i) {
        double u[3];
        u_array->GetTuple(i, u);
        for (int d = 0; d < 3; ++d)
            data.nodes[i].u[d] = u[d];
    }

    auto nCells = grid->GetNumberOfCells();
    data.elements.resize(nCells);

    for (vtkIdType c = 0; c < grid->GetNumberOfCells(); ++c) {
        vtkIdType npts;
        const vtkIdType* pts;
        grid->GetCellPoints(c, npts, pts);

        std::vector<size_t> cell(npts);
        for (vtkIdType i = 0; i < npts; ++i) cell[i] = pts[i];
        data.elements[c] = cell;
    }

    data.compute_dual_geometry();
    step++;
}

void PanelCSVProvider::load_next(SurfaceData& data, double& time) {

        int rank;
        MPI_Comm_rank(MPI_COMM_WORLD, &rank);

        time = current * dt;
        const std::string& fname = filenames[current++];

        std::ifstream in(fname);
        if (!in)
            throw std::runtime_error("PanelCSVProvider: cannot open " + fname);

        data.nodes.clear();
        data.elements.clear();   // dS set directly; no connectivity needed

        std::string line;
        bool have_meta = false;

        while (std::getline(in, line)) {

            size_t s = line.find_first_not_of(" \t\r\n");
            if (s == std::string::npos) continue;   // blank line
            if (line[s] == '#')         continue;   // banner / column header

            std::vector<double> v = parse_csv(line);

            // --- first numeric line is the metadata row ---
            if (!have_meta) {
                if (v.size() < 9)
                    throw std::runtime_error("PanelCSVProvider: bad metadata in " + fname);

                // time  = v[1];
                p0    = v[2];
                rho0  = v[3];
                c0    = v[4];
                U_inf = Vect3(v[6], v[7], v[8]);

                have_meta = true;
                continue;
            }

            // --- element row (35 fields) ---
            if (v.size() < 35)
                throw std::runtime_error("PanelCSVProvider: short element row in " + fname);

            double rho = v[7];
            double dp  = v[8];          // load across the sheet

            // --- mean panel: LOADING only ---
            Node mean;
            mean.x    = Vect3(v[0], v[1], v[2]);
            mean.dS   = Vect3(v[3], v[4], v[5]) * v[6];     // n_mean * area_mean
            mean.p    = p0 - dp;
            mean.u    = Vect3(v[12], v[13], v[14]);
            mean.rho  = rho;
            mean.role = NodeRole::LoadingOnly;

            // --- upper face: THICKNESS only ---
            Node up;
            up.x    = Vect3(v[15], v[16], v[17]);
            up.dS   = Vect3(v[21], v[22], v[23]) * v[27];   // n_u * area_u
            up.p    = p0;                                  // no loading
            up.u    = Vect3(v[29], v[30], v[31]);
            up.rho  = rho;
            up.role = NodeRole::ThicknessOnly;

            // --- lower face: THICKNESS only ---
            Node lo;
            lo.x    = Vect3(v[18], v[19], v[20]);
            lo.dS   = Vect3(v[24], v[25], v[26]) * v[28];   // n_l * area_l
            lo.p    = p0;
            lo.u    = Vect3(v[32], v[33], v[34]);
            lo.rho  = rho;
            lo.role = NodeRole::ThicknessOnly;

            data.nodes.push_back(mean);
            data.nodes.push_back(up);
            data.nodes.push_back(lo);
        }

        if (!have_meta)
            throw std::runtime_error("PanelCSVProvider: no data in " + fname);

        // --- progress (updates in place on stderr) ---
        if (rank == 0)
        {
            int total = (int)filenames.size();
            std::cerr << "\rPanelCSV progress: " << (100 * current / total) << "% "
                    << "(" << current << "/" << total << ")" << std::flush;
            if (current == total)
                std::cerr << std::endl;
        }

    }
