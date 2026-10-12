#include "es_PIC.h"

void evalAzimuthal(
    int nspec,
    const long *nparts,
    const int ngpt[3],
    const double (*part)[MAXNPARTPROPS],
    const double *particleWeight,
    std::vector<double> &number_density,
    std::vector<double> &vel_azimuthal,
    std::vector<double> &momentum_azimuthal,                // M_{rr}, M_{rz}, M_{theta,theta} (convective)
    std::vector<double> &stress_azimuthal,                  // S_{rr}, S_{rz}, S_{theta,theta} (thermal)
    std::vector<double> &ExB_component,                     // e n_e * (E_r/B_z) or -e n_e * (E_z/B_r)
    std::vector<double> &diamagnetic_component,             // (1/B) * \partial p_e
    std::vector<double> &inertial_component,                // (1/B) * (\nabla\cdot M_e terms)
    std::vector<double> &stressdiv_component,               // -(1/B) * (\nabla\cdot \pi_e terms)
    std::vector<double> &Je_theta_total_fluid_fullmomentum, // sum of the four
    std::vector<double> &Je_theta_PIC,                      // from particle data directly
    const std::vector<double> &histVelField,                // 6*nspec*nxyz; using u_z,u_r,u_theta
    double (*E)[3],
    double (*B)[3],
    int eqSelect, // 1: Eq.20; 2: Eq.21
    int volSelect // 1: node volume, volinv; 2: cell volume, calculated locally
)
{
    const int nx = ngpt[0], ny = ngpt[1], nz = ngpt[2];
    const int nxy = nx * ny, nxyz = nxy * nz;
    const double dz = dx[0], dr = dx[1];

    // safety check
    if (nz != 1)
    {
        esPIC_exit("evalAzimuthal: nz must be 1 for azimuthal evaluation in this version.");
    }

    // --- storage (internal_*) ---------------------------------
    double *internal_ndensity = new double[nspec * nxyz * 1];           // number density
    double *internal_vel_azimuthal = new double[nspec * nxyz * 1];      // <v_theta>
    double *internal_momentum_azimuthal = new double[nspec * nxyz * 4]; // M_rr, M_rz, M_theta,theta, M_zz (convective)
    double *internal_stress_tensor = new double[nspec * nxyz * 4];      // S_rr, S_rz, S_thetatheta, S_zz (thermal)
    double *internal_p = new double[nspec * nxyz * 1];                  // scalar pressure p
    double *internal_pi_tensor = new double[nspec * nxyz * 4];          // \pi_rr, \pi_rz, \pi_theta,theta. \pi_zz
    double *internal_ExB_comp = new double[nspec * nxyz * 1];           // E×B term
    double *internal_diamagnetic_comp = new double[nspec * nxyz * 1];   // \partial p_e term
    double *internal_inertial_comp = new double[nspec * nxyz * 1];      // \nabla\cdot M_e term
    double *internal_stress_comp = new double[nspec * nxyz * 1];        // -\nabla\cdot \pi_e term
    double *internal_J_theta_total = new double[nspec * nxyz * 1];      // total J_theta,theta from fluid eq.
    double *internal_J_theta_PIC = new double[nspec * nxyz * 1];        // J_theta from particle data directly

    memset(internal_ndensity, 0, nspec * nxyz * sizeof(double));
    memset(internal_vel_azimuthal, 0, nspec * nxyz * sizeof(double));
    memset(internal_momentum_azimuthal, 0, 4 * nspec * nxyz * sizeof(double));
    memset(internal_stress_tensor, 0, 4 * nspec * nxyz * sizeof(double));
    memset(internal_p, 0, nspec * nxyz * sizeof(double));
    memset(internal_pi_tensor, 0, 4 * nspec * nxyz * sizeof(double));
    memset(internal_ExB_comp, 0, nspec * nxyz * sizeof(double));
    memset(internal_diamagnetic_comp, 0, nspec * nxyz * sizeof(double));
    memset(internal_inertial_comp, 0, nspec * nxyz * sizeof(double));
    memset(internal_stress_comp, 0, nspec * nxyz * sizeof(double));
    memset(internal_J_theta_total, 0, nspec * nxyz * sizeof(double));
    memset(internal_J_theta_PIC, 0, nspec * nxyz * sizeof(double));

    long npartsTot = nparts[MAXNSPECTYPES];

#define igpt(i, j, k) ((k) * (nxy) + (j) * (nx) + (i))
#define idx1(nid, is) ((nid) * nspec + (is))       // scalar per cell,species
#define idx3(nid, is) (((nid) * nspec + (is)) * 3) // 3 comps per cell,species
#define idx4(nid, is) (((nid) * nspec + (is)) * 4) // 4 comps per cell,species
#define idx6(nid, is) (((nid) * nspec + (is)) * 6) // histVelField

    double partX[3], partV[3];
    double xpm, ypm;
    double ddx0, ddy0, ddx1, ddy1, w;
    double r2, rj2, rjp2, vinv; // axisymmetric
    double uz, ur, uth;
    int nid, ienty1, ienty4, ienty6;
    int ix, iy, iz = 0;

    // -------------------------------------------------------
    // PARTICLE LOOP: deposit n, <v_theta>, J_theta and thermal stresses S
    // -------------------------------------------------------

    for (long ipart = 0; ipart < npartsTot; ipart++)
    {
        int itype = (int)(part[ipart][ptype] + smallValue);
        double charge_number = spec[itype][specchg];
        for (int ip = 0; ip < 3; ip++)
        {
            partX[ip] = part[ipart][xp + ip];
            partV[ip] = part[ipart][up + ip];
        }

        /* in case very few particles move more than one cell per step. */
        /* in this case the time step dt should be reduced. */
        if (partX[0] >= domain[1])
        {
            partX[0] = domain[1] - smallValue;
        }
        if (partX[1] >= domain[3])
        {
            partX[1] = domain[3] - smallValue;
        }
        if (partX[2] >= domain[5])
        {
            partX[2] = domain[5] - smallValue;
        }
        if (partX[0] <= domain[0])
        {
            partX[0] = domain[0] + smallValue;
        }
        if (partX[1] <= domain[2])
        {
            partX[1] = domain[2] + smallValue;
        }
        if (partX[2] <= domain[4])
        {
            partX[2] = domain[4] + smallValue;
        }

        xpm = (partX[0] - domain[0]) * dxinv[0]; // z
        ypm = (partX[1] - domain[2]) * dxinv[1]; // r
        ix = (int)xpm;
        iy = (int)ypm;

        r2 = SQR(partX[1]);
        rj2 = iy * dx[1] + domain[2];
        rjp2 = rj2 + dx[1];
        rj2 *= rj2;
        rjp2 *= rjp2;
        vinv = PI * (rjp2 - rj2) * dx[0];
        vinv = 1. / vinv;

        ddx0 = xpm - ix;
        ddx1 = 1.0 - ddx0;
        ddx0 *= dx[0];
        ddx1 *= dx[0];
        ddy0 = PI * (r2 - rj2);
        ddy1 = PI * (rjp2 - r2);

        struct Node
        {
            int i, j;
            double w;
        };
        Node node[4] = {
            {ix, iy, ddx1 * ddy1 * vinv},        // i, j
            {ix + 1, iy, ddx0 * ddy1 * vinv},    // i+1, j
            {ix, iy + 1, ddx1 * ddy0 * vinv},    // i, j+1
            {ix + 1, iy + 1, ddx0 * ddy0 * vinv} // i+1, j+1
        };

        for (int nit = 0; nit < 4; nit++)
        {
            const int ii = node[nit].i;
            const int jj = node[nit].j;
            w = node[nit].w;

            nid = igpt(ii, jj, iz);
            // ienty = 1 * nspec * nid + 1 * itype;
            ienty1 = idx1(nid, itype);
            // ienty3 = idx3(nid, itype);
            ienty4 = idx4(nid, itype);
            // ienty2 = 6 * nspec * nid + 6 * itype;
            ienty6 = idx6(nid, itype);

            // bulk velocities from histVelField, already divided by ns
            uz = histVelField[ienty6 + 0];  // uz
            ur = histVelField[ienty6 + 1];  // ur
            uth = histVelField[ienty6 + 2]; // uth

            if (volSelect == 2)
            {
                // number density, need particle weight, do it later
                internal_ndensity[ienty1] += w * vinv;
                // azimuthal velocity, need particle weight, do it later
                internal_vel_azimuthal[ienty1] += w * partV[2] * vinv; // v_theta
                // PIC J_theta, need particle weight, do it later
                internal_J_theta_PIC[ienty1] += w * charge_number * partV[2] * vinv;

                // thermal stress S (no mass here; add me later when forming p, \pi)
                internal_stress_tensor[ienty4 + 0] += w * (partV[1] - ur) * (partV[1] - ur) * vinv;   // S_rr
                internal_stress_tensor[ienty4 + 1] += w * (partV[1] - ur) * (partV[0] - uz) * vinv;   // S_rz
                internal_stress_tensor[ienty4 + 2] += w * (partV[2] - uth) * (partV[2] - uth) * vinv; // S_theta,theta
                internal_stress_tensor[ienty4 + 3] += w * (partV[0] - uz) * (partV[0] - uz) * vinv;   // S_zz
            }
            else
            {
                // number density, need particle weight, do it later
                internal_ndensity[ienty1] += w;
                // azimuthal velocity, need particle weight, do it later
                internal_vel_azimuthal[ienty1] += w * partV[2]; // v_theta
                // PIC J_theta, need particle weight, do it later
                internal_J_theta_PIC[ienty1] += w * charge_number * partV[2];

                // thermal stress S (no mass here; add me later when forming p, \pi)
                internal_stress_tensor[ienty4 + 0] += w * (partV[1] - ur) * (partV[1] - ur);   // S_rr
                internal_stress_tensor[ienty4 + 1] += w * (partV[1] - ur) * (partV[0] - uz);   // S_rz
                internal_stress_tensor[ienty4 + 2] += w * (partV[2] - uth) * (partV[2] - uth); // S_theta,theta
                internal_stress_tensor[ienty4 + 3] += w * (partV[0] - uz) * (partV[0] - uz);   // S_zz
            }
        }
    }
    // END OF PARTICLE LOOP - EXB, DIAMAGNETIC COMPONENTS ARE FIELD-LEVEL CALCULATIONS

    // -------------------------------------------------------
    // CELL LOOP, build M, p, \pi; then derivatives
    // -------------------------------------------------------
    const double epsB = 1e-12;

    for (int nid = 0; nid < nxyz; nid++)
    {
        // get the <vtheta> by density
        for (int s = 0; s < nspec; s++)
        {
            if (internal_ndensity[idx1(nid, s)] > 0.0)
            {
                internal_vel_azimuthal[idx1(nid, s)] /= internal_ndensity[idx1(nid, s)]; // already implicitly weighted by particleWeight in deposition
            }
            internal_ndensity[idx1(nid, s)] *= particleWeight[s];
            internal_J_theta_PIC[idx1(nid, s)] *= particleWeight[s];
            internal_stress_tensor[idx4(nid, s) + 0] *= particleWeight[s];
            internal_stress_tensor[idx4(nid, s) + 1] *= particleWeight[s];
            internal_stress_tensor[idx4(nid, s) + 2] *= particleWeight[s];
            internal_stress_tensor[idx4(nid, s) + 3] *= particleWeight[s];
            if (volSelect == 1)
            {
                internal_ndensity[idx1(nid, s)] *= volinv[nid];
                internal_J_theta_PIC[idx1(nid, s)] *= volinv[nid];
                internal_stress_tensor[idx4(nid, s) + 0] *= volinv[nid];
                internal_stress_tensor[idx4(nid, s) + 1] *= volinv[nid];
                internal_stress_tensor[idx4(nid, s) + 2] *= volinv[nid];
                internal_stress_tensor[idx4(nid, s) + 3] *= volinv[nid];
            }
        }
    }

    // build M, scalar p_e, and \pi = pI - S
    for (int iy = 0; iy < ny; iy++)
    {
        for (int ix = 0; ix < nx; ix++)
        {
            const int nid = igpt(ix, iy, iz);

            for (int s = 0; s < nspec; s++)
            {
                const double nlocal = internal_ndensity[idx1(nid, s)]; // this internal density already get weighted by the volume factor
                if (nlocal <= 0.0)
                {
                    continue;
                }
                const double mass = spec[s][specmass];

                uz = histVelField[idx6(nid, s) + 0];
                ur = histVelField[idx6(nid, s) + 1];
                uth = histVelField[idx6(nid, s) + 2];

                // Convective momentum flux M = m*n*u*u
                internal_momentum_azimuthal[idx4(nid, s) + 0] = mass * nlocal * ur * ur;   // M_rr
                internal_momentum_azimuthal[idx4(nid, s) + 1] = mass * nlocal * ur * uz;   // M_rz
                internal_momentum_azimuthal[idx4(nid, s) + 2] = mass * nlocal * uth * uth; // M_theta,theta
                internal_momentum_azimuthal[idx4(nid, s) + 3] = mass * nlocal * uz * uz;   // M_zz

                // Scalar pressure p = (mass/3) * (S_rr + S_zz + S_theta,theta)
                const double Srr = internal_stress_tensor[idx4(nid, s) + 0];
                const double Stt = internal_stress_tensor[idx4(nid, s) + 2];
                const double Szz = internal_stress_tensor[idx4(nid, s) + 3];
                const double pressure = (mass / 3.0) * (Srr + Szz + Stt);
                internal_p[idx1(nid, s)] = pressure;

                // \pi = pI - mass*S   (off-diagonal: -mass*S_rz)
                const double Srz = internal_stress_tensor[idx4(nid, s) + 1];
                internal_pi_tensor[idx4(nid, s) + 0] = pressure - mass * Srr; // \pi_rr
                internal_pi_tensor[idx4(nid, s) + 1] = -mass * Srz;           // \pi_rz
                internal_pi_tensor[idx4(nid, s) + 2] = pressure - mass * Stt; // \pi_theta,theta
                internal_pi_tensor[idx4(nid, s) + 3] = pressure - mass * Szz; // \pi_zz
            }
        }
    }

    for (int ix = 0; ix < nx; ++ix)
    {
        const int iy = 0;
        const int nid = igpt(ix, iy, iz);
        for (int s = 0; s < nspec; ++s)
        {
            // // isotropy & odd terms
            // const double Mavg = 0.5 * (internal_momentum_azimuthal[idx4(nid, s) + 0] + internal_momentum_azimuthal[idx4(nid, s) + 2]);
            internal_momentum_azimuthal[idx4(nid, s) + 0] = 0.0;
            internal_momentum_azimuthal[idx4(nid, s) + 2] = 0.0;
            internal_momentum_azimuthal[idx4(nid, s) + 1] = 0.0; // Mrz

            // const double piavg = 0.5 * (internal_pi_tensor[idx4(nid, s) + 0] + internal_pi_tensor[idx4(nid, s) + 2]);
            internal_pi_tensor[idx4(nid, s) + 0] = 0.0;
            internal_pi_tensor[idx4(nid, s) + 2] = 0.0;
            internal_pi_tensor[idx4(nid, s) + 1] = 0.0; // pirz

            // velocities at axis:
            internal_vel_azimuthal[idx1(nid, s)] = 0.0; // <v_theta>(axis) = 0
            // PIC Jtheta at axis:
            internal_J_theta_PIC[idx1(nid, s)] = 0.0;
        }
    }

    // lambda inline function, finite differences centered interior and one-sided at edges
    auto d_dr = [&](const double *A, int s, int ix, int iy) -> double
    {
        const int nid = igpt(ix, iy, iz);
        if (iy == 0)
        {
            // const int i1 = igpt(ix, iy + 1, iz), i2 = igpt(ix, iy + 2, iz);
            // return (-3.0 * A[idx1(nid, s)] + 4.0 * A[idx1(i1, s)] - A[idx1(i2, s)]) / (2.0 * dr);
            //  keep it symmetric for r=0
            return 0.0;
        }
        if (iy == ny - 1)
        {
            const int i1 = igpt(ix, iy - 1, iz), i2 = igpt(ix, iy - 2, iz);
            return (3.0 * A[idx1(nid, s)] - 4.0 * A[idx1(i1, s)] + A[idx1(i2, s)]) / (2.0 * dr);
        }
        return (A[idx1(igpt(ix, iy + 1, iz), s)] - A[idx1(igpt(ix, iy - 1, iz), s)]) / (2.0 * dr);
    };
    auto d_dz = [&](const double *A, int s, int ix, int iy) -> double
    {
        const int nid = igpt(ix, iy, iz);
        if (ix == 0)
        {
            const int i1 = igpt(ix + 1, iy, iz), i2 = igpt(ix + 2, iy, iz);
            return (-3.0 * A[idx1(nid, s)] + 4.0 * A[idx1(i1, s)] - A[idx1(i2, s)]) / (2.0 * dz);
        }
        if (ix == nx - 1)
        {
            const int i1 = igpt(ix - 1, iy, iz), i2 = igpt(ix - 2, iy, iz);
            return (3.0 * A[idx1(nid, s)] - 4.0 * A[idx1(i1, s)] + A[idx1(i2, s)]) / (2.0 * dz);
        }
        return (A[idx1(igpt(ix + 1, iy, iz), s)] - A[idx1(igpt(ix - 1, iy, iz), s)]) / (2 * dz);
    };
    auto d4_dr = [&](const double *A4, int s, int ix, int iy, int c) -> double
    {
        const int nid = igpt(ix, iy, iz);
        if (iy == 0)
        {
            const int i1 = igpt(ix, iy + 1, iz), i2 = igpt(ix, iy + 2, iz);
            return (-3.0 * A4[idx4(nid, s) + c] + 4.0 * A4[idx4(i1, s) + c] - A4[idx4(i2, s) + c]) / (2.0 * dr);
        }
        if (iy == ny - 1)
        {
            const int i1 = igpt(ix, iy - 1, iz), i2 = igpt(ix, iy - 2, iz);
            return (3.0 * A4[idx4(nid, s) + c] - 4.0 * A4[idx4(i1, s) + c] + A4[idx4(i2, s) + c]) / (2.0 * dr);
        }
        return (A4[idx4(igpt(ix, iy + 1, iz), s) + c] - A4[idx4(igpt(ix, iy - 1, iz), s) + c]) / (2 * dr);
    };
    auto d4_dz = [&](const double *A4, int s, int ix, int iy, int c) -> double
    {
        const int nid = igpt(ix, iy, iz);
        if (ix == 0)
        {
            const int i1 = igpt(ix + 1, iy, iz), i2 = igpt(ix + 2, iy, iz);
            return (-3.0 * A4[idx4(nid, s) + c] + 4.0 * A4[idx4(i1, s) + c] - A4[idx4(i2, s) + c]) / (2.0 * dz);
        }
        if (ix == nx - 1)
        {
            const int i1 = igpt(ix - 1, iy, iz), i2 = igpt(ix - 2, iy, iz);
            return (3.0 * A4[idx4(nid, s) + c] - 4.0 * A4[idx4(i1, s) + c] + A4[idx4(i2, s) + c]) / (2.0 * dz);
        }
        return (A4[idx4(igpt(ix + 1, iy, iz), s) + c] - A4[idx4(igpt(ix - 1, iy, iz), s) + c]) / (2 * dz);
    };
    auto safe_div = [&](double num, double den)
    {
        if (std::abs(den) < epsB)
        {
            return 0.0; // or std::copysign(epsB, den)
        }
        return num / den;
    };

    // components per Hu 2021 (Eqs. 20 and 21)
    for (int iy = 0; iy < ny; iy++)
    {
        // const double rc = domain[2] + (iy + 0.5) * dr;
        const double rnode = domain[2] + iy * dr;
        const double invr = (rnode > 0.0 ? 1.0 / rnode : 0.0);

        for (int ix = 0; ix < nx; ix++)
        {
            const int nid = igpt(ix, iy, iz);
            const double Bz = B[nid][0], Br = B[nid][1];
            const double Ez = E[nid][0], Er = E[nid][1];

            // decide which equation to use here
            const bool use20 = (eqSelect == 1) || (eqSelect == 0 && fabs(Bz) > 0.005);
            const bool use21 = (eqSelect == 2) || (eqSelect == 0 && fabs(Bz) <= 0.005);

            for (int s = 0; s < nspec; s++)
            {
                const double nlocal = internal_ndensity[idx1(nid, s)];
                if (nlocal <= 0.0)
                {
                    // keep zeros explicit
                    internal_ExB_comp[idx1(nid, s)] = 0.0;
                    internal_diamagnetic_comp[idx1(nid, s)] = 0.0;
                    internal_inertial_comp[idx1(nid, s)] = 0.0;
                    internal_stress_comp[idx1(nid, s)] = 0.0;
                    continue;
                }
                // const double eabs = fabs(spec[s][specchg]);
                const double charge_number = spec[s][specchg];
                double exb = 0.0, dia = 0.0, ine = 0.0, str = 0.0;

                // diamagnetic component
                const double dp_dr = d_dr(internal_p, s, ix, iy);
                const double dp_dz = d_dz(internal_p, s, ix, iy);

                // Momentum flux derivatives
                const double dMrr_dr = d4_dr(internal_momentum_azimuthal, s, ix, iy, 0);
                const double dMrz_dz = d4_dz(internal_momentum_azimuthal, s, ix, iy, 1);
                const double dMrz_dr = d4_dr(internal_momentum_azimuthal, s, ix, iy, 1);
                const double dMzz_dz = d4_dz(internal_momentum_azimuthal, s, ix, iy, 3);

                // \pi-derivatives
                const double dPirr_dr = d4_dr(internal_pi_tensor, s, ix, iy, 0);
                const double dPirz_dz = d4_dz(internal_pi_tensor, s, ix, iy, 1);
                const double dPirz_dr = d4_dr(internal_pi_tensor, s, ix, iy, 1);
                const double dPizz_dz = d4_dz(internal_pi_tensor, s, ix, iy, 3);

                const double Mrr = internal_momentum_azimuthal[idx4(nid, s) + 0];
                const double Mrz = internal_momentum_azimuthal[idx4(nid, s) + 1];
                const double Mtt = internal_momentum_azimuthal[idx4(nid, s) + 2];
                const double Pirr = internal_pi_tensor[idx4(nid, s) + 0];
                const double Pirz = internal_pi_tensor[idx4(nid, s) + 1];
                const double Pitt = internal_pi_tensor[idx4(nid, s) + 2];

                // double MrrMttOVERr, PirrPittOVERr, MrzOVERr, PirzOVERr;
                // if (iy == 0)
                // {
                //     // force axisymmetry
                //     MrrMttOVERr = 0.0;
                //     PirrPittOVERr = 0.0;
                //     MrzOVERr = 0.0;
                //     PirzOVERr = 0.0;

                //     /* MrrMttOVERr = d3_dr(internal_momentum_azimuthal, s, ix, iy, 0) - d3_dr(internal_momentum_azimuthal, s, ix, iy, 2); // \partial r(Mrr - Mtt)|0
                //     PirrPittOVERr = d3_dr(internal_pi_tensor, s, ix, iy, 0) - d3_dr(internal_pi_tensor, s, ix, iy, 2);                 // \partial r(πrr - πtt)|0
                //     MrzOVERr = d3_dr(internal_momentum_azimuthal, s, ix, iy, 1);                                                       // \partial r(Mrz)|0
                //     PirzOVERr = d3_dr(internal_pi_tensor, s, ix, iy, 1); */
                //     // \partial r(πrz)|0
                // }
                // else
                // {
                // r>0: ordinary nodal 1/r factors
                double MrrMttOVERr = (Mrr - Mtt) * invr;
                double PirrPittOVERr = (Pirr - Pitt) * invr;
                double MrzOVERr = Mrz * invr;
                double PirzOVERr = Pirz * invr;
                // }

                // Eq.(20): divide by Bz
                if (use20)
                {
                    exb = -charge_number * nlocal * safe_div(Er, Bz);
                    dia = safe_div(dp_dr, Bz);
                    ine = safe_div(dMrr_dr + dMrz_dz + MrrMttOVERr, Bz);
                    str = -safe_div(dPirr_dr + dPirz_dz + PirrPittOVERr, Bz);
                }
                // Eq.(21): divide by Br
                if (use21)
                {
                    // prefer switching entirely to Eq.(21) when chosen
                    exb = charge_number * nlocal * safe_div(Ez, Br);
                    dia = -safe_div(dp_dz, Br);
                    ine = -safe_div(dMzz_dz + dMrz_dr + MrzOVERr, Br);
                    str = safe_div(dPizz_dz + dPirz_dr + PirzOVERr, Br);
                }

                internal_ExB_comp[idx1(nid, s)] = exb;
                internal_diamagnetic_comp[idx1(nid, s)] = dia;
                internal_inertial_comp[idx1(nid, s)] = ine;
                internal_stress_comp[idx1(nid, s)] = str;
            }
        }
    }

    // -----------------------
    // pack results to vectors
    // -----------------------
    number_density.clear();
    vel_azimuthal.clear();
    momentum_azimuthal.clear();
    stress_azimuthal.clear();
    ExB_component.clear();
    diamagnetic_component.clear();
    inertial_component.clear();
    stressdiv_component.clear();
    Je_theta_total_fluid_fullmomentum.clear();
    Je_theta_PIC.clear();

    number_density.reserve(1 * nspec * nxyz);
    vel_azimuthal.reserve(1 * nspec * nxyz);
    ExB_component.reserve(1 * nspec * nxyz);
    diamagnetic_component.reserve(1 * nspec * nxyz);
    inertial_component.reserve(1 * nspec * nxyz);
    stressdiv_component.reserve(1 * nspec * nxyz);
    Je_theta_total_fluid_fullmomentum.reserve(1 * nspec * nxyz);
    Je_theta_PIC.reserve(1 * nspec * nxyz);
    momentum_azimuthal.reserve(4 * nspec * nxyz);
    stress_azimuthal.reserve(4 * nspec * nxyz);

    for (int nid = 0; nid < nxyz; nid++)
    {
        for (int s = 0; s < nspec; s++)
        {
            number_density.push_back(internal_ndensity[idx1(nid, s)]);
            vel_azimuthal.push_back(internal_vel_azimuthal[idx1(nid, s)]);
            ExB_component.push_back(internal_ExB_comp[idx1(nid, s)]);
            diamagnetic_component.push_back(internal_diamagnetic_comp[idx1(nid, s)]);
            inertial_component.push_back(internal_inertial_comp[idx1(nid, s)]);
            stressdiv_component.push_back(internal_stress_comp[idx1(nid, s)]);
            Je_theta_total_fluid_fullmomentum.push_back(internal_ExB_comp[idx1(nid, s)] + internal_diamagnetic_comp[idx1(nid, s)] + internal_inertial_comp[idx1(nid, s)] + internal_stress_comp[idx1(nid, s)]);
            Je_theta_PIC.push_back(internal_J_theta_PIC[idx1(nid, s)]);

            momentum_azimuthal.push_back(internal_momentum_azimuthal[idx4(nid, s) + 0]);
            momentum_azimuthal.push_back(internal_momentum_azimuthal[idx4(nid, s) + 1]);
            momentum_azimuthal.push_back(internal_momentum_azimuthal[idx4(nid, s) + 2]);
            momentum_azimuthal.push_back(internal_momentum_azimuthal[idx4(nid, s) + 3]);
            stress_azimuthal.push_back(internal_stress_tensor[idx4(nid, s) + 0]); // S_rr
            stress_azimuthal.push_back(internal_stress_tensor[idx4(nid, s) + 1]); // S_rz
            stress_azimuthal.push_back(internal_stress_tensor[idx4(nid, s) + 2]); // S_theta,theta
            stress_azimuthal.push_back(internal_stress_tensor[idx4(nid, s) + 3]); // S_zz
        }
    }

#undef igpt
#undef idx1
#undef idx3
#undef idx4
#undef idx6

    delete[] internal_ndensity;
    delete[] internal_vel_azimuthal;
    delete[] internal_momentum_azimuthal;
    delete[] internal_stress_tensor;
    delete[] internal_p;
    delete[] internal_pi_tensor;
    delete[] internal_ExB_comp;
    delete[] internal_diamagnetic_comp;
    delete[] internal_inertial_comp;
    delete[] internal_stress_comp;
    delete[] internal_J_theta_total;
    delete[] internal_J_theta_PIC;

    return;
}
