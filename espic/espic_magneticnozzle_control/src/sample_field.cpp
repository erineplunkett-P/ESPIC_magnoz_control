#include <iostream>
#include <fstream>
#include <sstream>
#include <cstdio>
#include <cstdlib>
#include <sys/stat.h>
#include <cstring>

#include "sample_field.h"
#include "pic_math.h"
#include "pic_io.h"

extern const double PI;
extern double *volinv;

// Constructor
SampleField::SampleField(const std::string &ifilename, int step, int num_spec, int ngpt_in[3])
{
    read_input(ifilename);

    nspec = num_spec;
    nsamples = 0;
    ngpt[0] = ngpt_in[0];
    ngpt[1] = ngpt_in[1];
    ngpt[2] = ngpt_in[2];
    nnodes = ngpt[0] * ngpt[1] * ngpt[2];
    spfield = new Sample[nspec * nnodes];
    if (!read_sample(step))
    {
        nsamples = 0;
        step_now = step;
        memset(spfield, 0, nspec * nnodes * sizeof(Sample));
    }
}

SampleField::~SampleField()
{
    delete[] spfield;
}

// public methods
bool SampleField::sample(int step, long nptot, double (*particles)[MAXNPARTPROPS])
{
    if (step <= step_start)
        return false; // sampling has NOT started

    if (step % nevery_sample == 0)
    {
        double boxlo[3] = {domain[0], domain[2], domain[4]};
        double boxhi[3] = {domain[1], domain[3], domain[5]};
        int ix, iy;
        double xpm, ypm;
        double ddx0, ddy0, ddx1, ddy1, w;
        double r2, rj2, rjp2, vinv; // axisymmetric
        long nid;

        step_now = step;
        ++nsamples;
        for (long ip = 0; ip < nptot; ++ip)
        {
            double r[3] = {particles[ip][xp], particles[ip][yp], particles[ip][zp]};
            double v[3] = {particles[ip][up], particles[ip][vp], particles[ip][wp]};
            int ispec = (int)(particles[ip][ptype] + 1e-8);

            for (int a = 0; a < 3; ++a)
            {
                if (r[a] >= boxhi[a])
                    r[a] = boxhi[a] - smallValue;
            }
            xpm = (r[0] - boxlo[0]) * dxinv[0];
            ypm = (r[1] - boxlo[1]) * dxinv[1];
            ix = (int)xpm;
            iy = (int)ypm;

            r2 = SQR(r[1]);
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

            if (ngpt[2] == 1)
            {
                nid = iy * ngpt[0] + ix;
                w = ddx1 * ddy1 * vinv;
                compute_moment_sum(&(spfield[ispec * nnodes + nid]), v, w);

                nid = iy * ngpt[0] + ix + 1;
                w = ddx0 * ddy1 * vinv;
                compute_moment_sum(&(spfield[ispec * nnodes + nid]), v, w);

                nid = (iy + 1) * ngpt[0] + ix;
                w = ddx1 * ddy0 * vinv;
                compute_moment_sum(&(spfield[ispec * nnodes + nid]), v, w);

                nid = (iy + 1) * ngpt[0] + ix + 1;
                w = ddx0 * ddy0 * vinv;
                compute_moment_sum(&(spfield[ispec * nnodes + nid]), v, w);
            }
            else
            {
                esPIC_exit("Sample field for 3D domain has NOT implemented");
            }
        }
    } // end if (step % nevery_sample == 0)

    if (step % nevery_write == 0)
        write_sample();

    return (step % nevery_postproc == 0);
}

void SampleField::postproc_field(const double dt, const bool HYBRID,
                                 const Vec &pin,
                                 const MapF &fixedBV,
                                 SpMat &Mat_A,
                                 PtrPoissonSolver poisson_solver,
                                 double (*B)[3],
                                 int eqSelect)
{
    if (ngpt[2] != 1)
    {
        esPIC_exit("Postprocess sample field for 3D domain has NOT implemented");
    }

    std::cout << "postprocessing field with " << nsamples << " samples";

    double **chgDenst = new double *[nnodes];
    for (long nid = 0; nid < nnodes; nid++)
    {
        chgDenst[nid] = new double[nspec + 1];
        for (int ispec = 0; ispec <= nspec; ++ispec)
            chgDenst[nid][ispec] = 0.0;
    }
    double (*E)[3] = new double[nnodes][3];
    memset(&(E[0][0]), 0, nnodes * 3 * sizeof(double));

    // only for 2D
    Sample *ptrsp;
    for (long nid = 0; nid < nnodes; ++nid)
    {
        for (int ispec = 0; ispec < nspec; ++ispec)
        {
            ptrsp = &spfield[ispec * nnodes + nid];
            chgDenst[nid][ispec] = (ptrsp->NSUM) * partWeight[ispec] * spec[ispec][specchg] * volinv[nid] / nsamples;
            chgDenst[nid][nspec] += chgDenst[nid][ispec];
        }
    }

    // output charge density field
    std::ofstream outfile;
    outfile.open("hist_ave_chgDenst.dat", std::ofstream::out | std::ofstream::app);
    outfile << step_now * dt << " "
            << ngpt[0] << " " << ngpt[1] << " " << ngpt[2] << " " << nspec + 1;
    for (long nid = 0; nid < nnodes; ++nid)
        for (int ispec = 0; ispec <= nspec; ++ispec)
            outfile << " " << chgDenst[nid][ispec];
    outfile << "\n";
    outfile.close();
    outfile.clear();

    Vec p(pin), b(nnodes);

    if (HYBRID)
    {
        solve_field_nonlinear(ngpt, nspec, Mat_A, fixedBV, chgDenst, p, b, E, E, poisson_solver);
    }
    else
    {
        solve_field_linear(ngpt, nspec, Mat_A, fixedBV, chgDenst, p, b, E, E, poisson_solver);
    }

    // output potential field
    outfile.open("hist_ave_phi.dat", std::ofstream::out | std::ofstream::app);
    outfile << step_now * dt << " "
            << ngpt[0] << " " << ngpt[1] << " " << ngpt[2];
    for (long nid = 0; nid < nnodes; ++nid)
        outfile << " " << p[nid];
    outfile << "\n";
    outfile.close();
    outfile.clear();

    // output velField
    outfile.open("hist_ave_velField.dat", std::ofstream::out | std::ofstream::app);
    outfile << step_now * dt << " "
            << ngpt[0] << " " << ngpt[1] << " " << ngpt[2] << " " << nspec;
    for (long nid = 0; nid < nnodes; ++nid)
    {
        for (int ispec = 0; ispec < nspec; ++ispec)
        {
            double v[3] = {0.0, 0.0, 0.0};
            double v2[3] = {0.0, 0.0, 0.0};

            ptrsp = &spfield[ispec * nnodes + nid];
            if (ptrsp->NSUM > 0)
            {
                for (int a = 0; a < 3; ++a)
                {
                    v[a] = ptrsp->VSUM[a] / ptrsp->NSUM;
                    v2[a] = ptrsp->V2SUM[a] / ptrsp->NSUM;
                }
            }

            for (int a = 0; a < 3; ++a)
            {
                outfile << " " << v[a];
            }
            for (int a = 0; a < 3; ++a)
            {
                outfile << " " << v2[a];
            }
        }
    }
    outfile << "\n";
    outfile.close();
    outfile.clear();

    // new hist
    outfile.open("hist_ave_new_velField.dat", std::ofstream::out | std::ofstream::app);
    outfile << step_now * dt << " "
            << ngpt[0] << " " << ngpt[1] << " " << ngpt[2] << " " << nspec;

    for (long nid = 0; nid < nnodes; ++nid)
    {
        for (int ispec = 0; ispec < nspec; ++ispec)
        {
            double v3m[5] = {0.0, 0.0, 0.0, 0.0, 0.0};

            ptrsp = &spfield[ispec * nnodes + nid];
            if (ptrsp->NSUM > 0)
            {
                for (int a = 0; a < 5; ++a)
                    v3m[a] = ptrsp->V3SUM[a] / ptrsp->NSUM;
            }

            for (int a = 0; a < 5; ++a)
            {
                outfile << " " << v3m[a];
            }
        }
    }
    outfile << "\n";
    outfile.close();
    outfile.clear();

    // new new hist
    outfile.open("hist_ave_new_new_velField.dat", std::ofstream::out | std::ofstream::app);
    outfile << step_now * dt << " "
            << ngpt[0] << " " << ngpt[1] << " " << ngpt[2] << " " << nspec;

    for (long nid = 0; nid < nnodes; ++nid)
    {
        for (int ispec = 0; ispec < nspec; ++ispec)
        {
            double v4m[3] = {0.0, 0.0, 0.0};

            ptrsp = &spfield[ispec * nnodes + nid];
            if (ptrsp->NSUM > 0)
            {
                for (int a = 0; a < 3; ++a)
                    v4m[a] = ptrsp->V4SUM[a] / ptrsp->NSUM;
            }

            for (int a = 0; a < 3; ++a)
            {
                outfile << " " << v4m[a];
            }
        }
    }
    outfile << "\n";
    outfile.close();
    outfile.clear();

    char c[3] = {'x', 'y', 'z'};
    long nid;

    // vx, vy field for each species
    for (int a = 0; a < 2; ++a)
    {
        for (int ispec = 0; ispec < nspec; ++ispec)
        {
            nid = 0;
            std::ostringstream fname;
            fname << "v" << c[a] << "_s" << ispec << ".dat";
            outfile.open(fname.str().c_str(), std::ofstream::out);

            for (int iy = 0; iy < ngpt[1]; ++iy)
            {
                for (int ix = 0; ix < ngpt[0]; ++ix)
                {
                    double v = 0.0;

                    ptrsp = &spfield[ispec * nnodes + nid];
                    if (ptrsp->NSUM > 0)
                        v = ptrsp->VSUM[a] / ptrsp->NSUM;

                    outfile << v << " ";
                    ++nid;
                }
                outfile << "\n";
            }
            outfile.close();
            outfile.clear();
        }
    }

    // Tx, Ty field for each species
    for (int a = 0; a < 2; ++a)
    {
        for (int ispec = 0; ispec < nspec; ++ispec)
        {
            nid = 0;
            std::ostringstream fname;
            fname << "T" << c[a] << "_s" << ispec << ".dat";
            outfile.open(fname.str().c_str(), std::ofstream::out);

            for (int iy = 0; iy < ngpt[1]; ++iy)
            {
                for (int ix = 0; ix < ngpt[0]; ++ix)
                {
                    double v = 0.0, v2 = 0.0;
                    double m = spec[ispec][specmass];

                    ptrsp = &spfield[ispec * nnodes + nid];
                    if (ptrsp->NSUM > 0)
                    {
                        v = ptrsp->VSUM[a] / ptrsp->NSUM;
                        v2 = ptrsp->V2SUM[a] / ptrsp->NSUM;
                    }

                    outfile << m * (v2 - v * v) << " ";
                    ++nid;
                }
                outfile << "\n";
            }
            outfile.close();
            outfile.clear();
        }
    }

    if (B != nullptr && nsamples > 0)
    {
        std::vector<double> number_density;
        std::vector<double> vel_azimuthal;
        std::vector<double> momentum_azimuthal;
        std::vector<double> stress_azimuthal;
        std::vector<double> ExB_component;
        std::vector<double> diamagnetic_component;
        std::vector<double> inertial_component;
        std::vector<double> stressdiv_component;
        std::vector<double> Je_theta_total_fluid_fullmomentum;
        std::vector<double> Je_theta_PIC;

        compute_time_averaged_azimuthal(E,
                                        B,
                                        eqSelect,
                                        number_density,
                                        vel_azimuthal,
                                        momentum_azimuthal,
                                        stress_azimuthal,
                                        ExB_component,
                                        diamagnetic_component,
                                        inertial_component,
                                        stressdiv_component,
                                        Je_theta_total_fluid_fullmomentum,
                                        Je_theta_PIC);

        write_time_averaged_azimuthal(dt,
                                      number_density,
                                      vel_azimuthal,
                                      momentum_azimuthal,
                                      stress_azimuthal,
                                      ExB_component,
                                      diamagnetic_component,
                                      inertial_component,
                                      stressdiv_component,
                                      Je_theta_total_fluid_fullmomentum,
                                      Je_theta_PIC);
    }

    std::cout << "\n";

    delete[] E;
    for (long nid = 0; nid < nnodes; nid++)
    {
        delete[] chgDenst[nid];
        chgDenst[nid] = NULL;
    }
    delete[] chgDenst;
}

// private methods
void SampleField::read_input(const std::string &ifilename)
{
    std::ifstream infile(ifilename.c_str());
    if (!infile)
    {
        std::cout << "Do NOT sample field." << std::endl;
        return;
    }

    std::cout << "Reading [" << ifilename << "]..." << std::endl;
    // first line
    std::string line;
    std::getline(infile, line);
    std::istringstream iss(line);
    iss >> step_start;

    std::getline(infile, line);
    iss.str(line);
    iss >> nevery_sample;

    std::getline(infile, line);
    iss.str(line);
    iss >> nevery_write;

    std::getline(infile, line);
    iss.str(line);
    iss >> nevery_postproc;
}

bool SampleField::read_sample(int step)
{
    if (!is_file_exist("esPIC.sample"))
    {
        return false;
    }

    // size compatibility check
    struct stat st;
    if (stat("esPIC.sample", &st) != 0)
    {
        return false;
    }

    const long expected =
        (long)sizeof(int) +
        (long)sizeof(uint32_t) +
        (long)sizeof(Sample) * (long)nspec * (long)nnodes;

    if (st.st_size != expected)
    {
        // Old Sample layout or corrupted file -> reset gracefully
        std::cout << "Sample file layout mismatch. Resetting sampled field...\n";
        return false;
    }

    FILE *ifp = fopen("esPIC.sample", "rb");
    if (ifp == NULL)
    {
        std::ostringstream oss;
        oss << "Failed to open sample file [esPIC.sample]";
        esPIC_exit(oss.str());
    }

    long size_read = 0;
    size_read = fread(&step_now, sizeof(int), 1, ifp);
    size_read = fread(&nsamples, sizeof(uint32_t), 1, ifp);

    size_read = fread(spfield, sizeof(Sample), nspec * nnodes, ifp);
    fclose(ifp);
    if (size_read != nspec * nnodes)
    {
        std::ostringstream oss;
        oss << "Incompatible or corrupted sample file [esPIC.sample]";
        esPIC_exit(oss.str());
    }

    if (step_start >= step)
    {
        printf("Time when sampling starts changed: Resetting sampled field ... \n");
        nsamples = 0;
        step_now = step;
        memset(spfield, 0, nspec * nnodes * sizeof(Sample));
        return true;
    }

    if (step_now > step)
    {
        printf("--> Waring: Time when last sample made is later than the current time. <--\n");
        step_now = step;
    }

    return true;
}

void SampleField::write_sample()
{
    FILE *ofp = fopen("esPIC.sample", "wb");

    if (ofp == NULL)
    {
        perror("Error opening esPIC.sample");
        return;
    }

    fwrite(&step_now, sizeof(int), 1, ofp);
    fwrite(&nsamples, sizeof(uint32_t), 1, ofp);
    fwrite(spfield, sizeof(Sample), nspec * nnodes, ofp);

    fclose(ofp);
}

void SampleField::compute_moment_sum(Sample *ptrsp, double *v, double wgt)
{
    ptrsp->NSUM += wgt;
    for (int a = 0; a < 3; ++a)
    {
        ptrsp->VSUM[a] += v[a] * wgt;
        ptrsp->V2SUM[a] += v[a] * v[a] * wgt;
    }
    const double v0 = v[0];
    const double v1 = v[1];

    // ptrsp->VZVRSUM += (v0 * v1) * wgt;

    // Matches evalNewHistVelField ordering
    ptrsp->V3SUM[0] += (v0 * v0 * v0) * wgt;
    ptrsp->V3SUM[1] += (v1 * v1 * v1) * wgt;
    ptrsp->V3SUM[2] += (v0 * v0 * v1) * wgt;
    ptrsp->V3SUM[3] += (v0 * v1 * v1) * wgt;
    ptrsp->V3SUM[4] += (v0 * v1) * wgt;

    // Matches evalNewNewHistVelField ordering
    ptrsp->V4SUM[0] += (v0 * v0 * v0 * v0) * wgt;
    ptrsp->V4SUM[1] += (v1 * v1 * v1 * v1) * wgt;
    ptrsp->V4SUM[2] += (v0 * v0 * v1 * v1) * wgt;
    return;
}

void SampleField::compute_time_averaged_azimuthal(double (*E)[3],
                                                  double (*B)[3],
                                                  int eqSelect,
                                                  std::vector<double> &number_density,
                                                  std::vector<double> &vel_azimuthal,
                                                  std::vector<double> &momentum_azimuthal,
                                                  std::vector<double> &stress_azimuthal,
                                                  std::vector<double> &ExB_component,
                                                  std::vector<double> &diamagnetic_component,
                                                  std::vector<double> &inertial_component,
                                                  std::vector<double> &stressdiv_component,
                                                  std::vector<double> &Je_theta_total_fluid_fullmomentum,
                                                  std::vector<double> &Je_theta_PIC)
{
    const int nx = ngpt[0], ny = ngpt[1], nz = ngpt[2];
    const int nxy = nx * ny;
    const int nxyz = nxy * nz;

    if (nz != 1)
    {
        esPIC_exit("SampleField::compute_time_averaged_azimuthal: nz must be 1.");
    }

    const double dz = dx[0], dr = dx[1];
    const double epsB = 1e-12;

#define igpt(i, j, k) ((k) * (nxy) + (j) * (nx) + (i))
#define idx1(nid, is) ((nid) * nspec + (is))
#define idx4(nid, is) (((nid) * nspec + (is)) * 4)
#define idx6(nid, is) (((nid) * nspec + (is)) * 6)

    // Internal arrays
    double *nd = new double[nspec * nxyz];
    double *uth = new double[nspec * nxyz];
    double *M = new double[4 * nspec * nxyz];
    double *S = new double[4 * nspec * nxyz];
    double *p = new double[nspec * nxyz];
    double *piT = new double[4 * nspec * nxyz];
    double *ExB = new double[nspec * nxyz];
    double *Dia = new double[nspec * nxyz];
    double *Ine = new double[nspec * nxyz];
    double *Str = new double[nspec * nxyz];
    double *Jpic = new double[nspec * nxyz];
    double *histVelAve = new double[6 * nspec * nxyz];

    memset(histVelAve, 0, 6 * nspec * nxyz * sizeof(double));
    memset(nd, 0, nspec * nxyz * sizeof(double));
    memset(uth, 0, nspec * nxyz * sizeof(double));
    memset(M, 0, 4 * nspec * nxyz * sizeof(double));
    memset(S, 0, 4 * nspec * nxyz * sizeof(double));
    memset(p, 0, nspec * nxyz * sizeof(double));
    memset(piT, 0, 4 * nspec * nxyz * sizeof(double));
    memset(ExB, 0, nspec * nxyz * sizeof(double));
    memset(Dia, 0, nspec * nxyz * sizeof(double));
    memset(Ine, 0, nspec * nxyz * sizeof(double));
    memset(Str, 0, nspec * nxyz * sizeof(double));
    memset(Jpic, 0, nspec * nxyz * sizeof(double));

    if (nsamples == 0)
    {
        // nothing to do
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

        delete[] nd;
        delete[] uth;
        delete[] M;
        delete[] S;
        delete[] p;
        delete[] piT;
        delete[] ExB;
        delete[] Dia;
        delete[] Ine;
        delete[] Str;
        delete[] Jpic;
        delete[] histVelAve;
        return;
    }

    // Build time-averaged moments from spfield
    for (long nid = 0; nid < nnodes; ++nid)
    {
        for (int s = 0; s < nspec; ++s)
        {
            Sample *ptrsp = &spfield[s * nnodes + nid];
            if (ptrsp->NSUM <= 0.0)
            {
                continue;
            }

            // Time-averaged number density
            const double nlocal = ptrsp->NSUM * partWeight[s] * volinv[nid] / (double)nsamples;
            nd[idx1(nid, s)] = nlocal;

            // Mean velocities
            const double uz_local = ptrsp->VSUM[0] / ptrsp->NSUM;
            const double ur_local = ptrsp->VSUM[1] / ptrsp->NSUM;
            const double uth_local = ptrsp->VSUM[2] / ptrsp->NSUM;
            uth[idx1(nid, s)] = uth_local;

            // Fill histVelAve
            histVelAve[idx6(nid, s) + 0] = uz_local;
            histVelAve[idx6(nid, s) + 1] = ur_local;
            histVelAve[idx6(nid, s) + 2] = uth_local;

            // Raw second moments
            const double vz2 = ptrsp->V2SUM[0] / ptrsp->NSUM;
            const double vr2 = ptrsp->V2SUM[1] / ptrsp->NSUM;
            const double vt2 = ptrsp->V2SUM[2] / ptrsp->NSUM;
            const double vzvr = ptrsp->V3SUM[4] / ptrsp->NSUM;
            // Thermal Stress tensor components
            // S_ij = n(<v_i v_j> - u_i u_j)
            const double Srr = nlocal * (vr2 - ur_local * ur_local);
            const double Srz = nlocal * (vzvr - ur_local * uz_local);
            const double Stt = nlocal * (vt2 - uth_local * uth_local);
            const double Szz = nlocal * (vz2 - uz_local * uz_local);
            S[idx4(nid, s) + 0] = Srr;
            S[idx4(nid, s) + 1] = Srz;
            S[idx4(nid, s) + 2] = Stt;
            S[idx4(nid, s) + 3] = Szz;

            // Convective momentum flux M = m*n*u*u
            const double mass = spec[s][specmass];
            M[idx4(nid, s) + 0] = mass * nlocal * ur_local * ur_local;
            M[idx4(nid, s) + 1] = mass * nlocal * ur_local * uz_local;
            M[idx4(nid, s) + 2] = mass * nlocal * uth_local * uth_local;
            M[idx4(nid, s) + 3] = mass * nlocal * uz_local * uz_local;

            // Scalar pressure p = (m/3)(Srr + Szz + Stt)
            const double pressure = (mass / 3.0) * (Srr + Szz + Stt);
            p[idx1(nid, s)] = pressure;

            // pi = pI - m*S
            piT[idx4(nid, s) + 0] = pressure - mass * Srr;
            piT[idx4(nid, s) + 1] = -mass * Srz;
            piT[idx4(nid, s) + 2] = pressure - mass * Stt;
            piT[idx4(nid, s) + 3] = pressure - mass * Szz;

            // Jtheta_PIC from averaged moments
            const double charge_number = spec[s][specchg];
            Jpic[idx1(nid, s)] = charge_number * nlocal * uth_local;
        }
    }

    // Axis enforcement (iy=0)
    for (int ix = 0; ix < nx; ++ix)
    {
        const int iy = 0;
        const int nid = igpt(ix, iy, 0);

        for (int s = 0; s < nspec; ++s)
        {
            M[idx4(nid, s) + 0] = 0.0;
            M[idx4(nid, s) + 1] = 0.0;
            M[idx4(nid, s) + 2] = 0.0;

            piT[idx4(nid, s) + 0] = 0.0;
            piT[idx4(nid, s) + 1] = 0.0;
            piT[idx4(nid, s) + 2] = 0.0;

            uth[idx1(nid, s)] = 0.0;
            histVelAve[idx6(nid, s) + 2] = 0.0;
            Jpic[idx1(nid, s)] = 0.0;
        }
    }

    // Derivative lambda functions
    auto d_dr = [&](const double *A, int s, int ix, int iy) -> double
    {
        const int nid = igpt(ix, iy, 0);
        if (iy == 0)
        {
            return 0.0;
        }
        if (iy == ny - 1)
        {
            const int i1 = igpt(ix, iy - 1, 0), i2 = igpt(ix, iy - 2, 0);
            return (3.0 * A[idx1(nid, s)] - 4.0 * A[idx1(i1, s)] + A[idx1(i2, s)]) / (2.0 * dr);
        }
        return (A[idx1(igpt(ix, iy + 1, 0), s)] - A[idx1(igpt(ix, iy - 1, 0), s)]) / (2.0 * dr);
    };
    auto d_dz = [&](const double *A, int s, int ix, int iy) -> double
    {
        const int nid = igpt(ix, iy, 0);
        if (ix == 0)
        {
            const int i1 = igpt(ix + 1, iy, 0), i2 = igpt(ix + 2, iy, 0);
            return (-3.0 * A[idx1(nid, s)] + 4.0 * A[idx1(i1, s)] - A[idx1(i2, s)]) / (2.0 * dz);
        }
        if (ix == nx - 1)
        {
            const int i1 = igpt(ix - 1, iy, 0), i2 = igpt(ix - 2, iy, 0);
            return (3.0 * A[idx1(nid, s)] - 4.0 * A[idx1(i1, s)] + A[idx1(i2, s)]) / (2.0 * dz);
        }
        return (A[idx1(igpt(ix + 1, iy, 0), s)] - A[idx1(igpt(ix - 1, iy, 0), s)]) / (2.0 * dz);
    };
    auto d4_dr = [&](const double *A4, int s, int ix, int iy, int c) -> double
    {
        const int nid = igpt(ix, iy, 0);
        if (iy == 0)
        {
            const int i1 = igpt(ix, iy + 1, 0), i2 = igpt(ix, iy + 2, 0);
            return (-3.0 * A4[idx4(nid, s) + c] + 4.0 * A4[idx4(i1, s) + c] - A4[idx4(i2, s) + c]) / (2.0 * dr);
        }
        if (iy == ny - 1)
        {
            const int i1 = igpt(ix, iy - 1, 0), i2 = igpt(ix, iy - 2, 0);
            return (3.0 * A4[idx4(nid, s) + c] - 4.0 * A4[idx4(i1, s) + c] + A4[idx4(i2, s) + c]) / (2.0 * dr);
        }
        return (A4[idx4(igpt(ix, iy + 1, 0), s) + c] - A4[idx4(igpt(ix, iy - 1, 0), s) + c]) / (2.0 * dr);
    };
    auto d4_dz = [&](const double *A4, int s, int ix, int iy, int c) -> double
    {
        const int nid = igpt(ix, iy, 0);
        if (ix == 0)
        {
            const int i1 = igpt(ix + 1, iy, 0), i2 = igpt(ix + 2, iy, 0);
            return (-3.0 * A4[idx4(nid, s) + c] + 4.0 * A4[idx4(i1, s) + c] - A4[idx4(i2, s) + c]) / (2.0 * dz);
        }
        if (ix == nx - 1)
        {
            const int i1 = igpt(ix - 1, iy, 0), i2 = igpt(ix - 2, iy, 0);
            return (3.0 * A4[idx4(nid, s) + c] - 4.0 * A4[idx4(i1, s) + c] + A4[idx4(i2, s) + c]) / (2.0 * dz);
        }
        return (A4[idx4(igpt(ix + 1, iy, 0), s) + c] - A4[idx4(igpt(ix - 1, iy, 0), s) + c]) / (2.0 * dz);
    };
    auto safe_div = [&](double num, double den) -> double
    {
        if (std::abs(den) < epsB)
        {
            return 0.0;
        }
        return num / den;
    };

    // Compute components
    for (int iy = 0; iy < ny; ++iy)
    {
        const double rnode = domain[2] + iy * dr;
        const double invr = (rnode > 0.0 ? 1.0 / rnode : 0.0);

        for (int ix = 0; ix < nx; ++ix)
        {
            const int nid = igpt(ix, iy, 0);

            const double Bz = B[nid][0], Br = B[nid][1];
            const double Ez = E[nid][0], Er = E[nid][1];

            const bool use20 = (eqSelect == 1) || (eqSelect == 0 && fabs(Bz) > 0.005);
            const bool use21 = (eqSelect == 2) || (eqSelect == 0 && fabs(Bz) <= 0.005);
            for (int s = 0; s < nspec; ++s)
            {
                const double nlocal = nd[idx1(nid, s)];
                if (nlocal <= 0.0)
                {
                    continue;
                }
                const double charge_number = spec[s][specchg];
                const double dp_dr = d_dr(p, s, ix, iy);
                const double dp_dz = d_dz(p, s, ix, iy);
                const double dMrr_dr = d4_dr(M, s, ix, iy, 0);
                const double dMrz_dz = d4_dz(M, s, ix, iy, 1);
                const double dMrz_dr = d4_dr(M, s, ix, iy, 1);
                const double dMzz_dz = d4_dz(M, s, ix, iy, 3);
                const double dPirr_dr = d4_dr(piT, s, ix, iy, 0);
                const double dPirz_dz = d4_dz(piT, s, ix, iy, 1);
                const double dPirz_dr = d4_dr(piT, s, ix, iy, 1);
                const double dPizz_dz = d4_dz(piT, s, ix, iy, 3);
                const double Mrr = M[idx4(nid, s) + 0];
                const double Mrz = M[idx4(nid, s) + 1];
                const double Mtt = M[idx4(nid, s) + 2];
                const double Pirr = piT[idx4(nid, s) + 0];
                const double Pirz = piT[idx4(nid, s) + 1];
                const double Pitt = piT[idx4(nid, s) + 2];
                const double MrrMttOVERr = (Mrr - Mtt) * invr;
                const double PirrPittOVERr = (Pirr - Pitt) * invr;
                const double MrzOVERr = Mrz * invr;
                const double PirzOVERr = Pirz * invr;
                double exb = 0.0, dia = 0.0, ine = 0.0, strc = 0.0;
                if (use20)
                {
                    exb = -charge_number * nlocal * safe_div(Er, Bz);
                    dia = safe_div(dp_dr, Bz);
                    ine = safe_div(dMrr_dr + dMrz_dz + MrrMttOVERr, Bz);
                    strc = -safe_div(dPirr_dr + dPirz_dz + PirrPittOVERr, Bz);
                }
                if (use21)
                {
                    exb = charge_number * nlocal * safe_div(Ez, Br);
                    dia = -safe_div(dp_dz, Br);
                    ine = -safe_div(dMzz_dz + dMrz_dr + MrzOVERr, Br);
                    strc = safe_div(dPizz_dz + dPirz_dr + PirzOVERr, Br);
                }
                ExB[idx1(nid, s)] = exb;
                Dia[idx1(nid, s)] = dia;
                Ine[idx1(nid, s)] = ine;
                Str[idx1(nid, s)] = strc;
            }
        }
    }

    // PACK RESULTS
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

    number_density.reserve(nspec * nxyz);
    vel_azimuthal.reserve(nspec * nxyz);
    ExB_component.reserve(nspec * nxyz);
    diamagnetic_component.reserve(nspec * nxyz);
    inertial_component.reserve(nspec * nxyz);
    stressdiv_component.reserve(nspec * nxyz);
    Je_theta_total_fluid_fullmomentum.reserve(nspec * nxyz);
    Je_theta_PIC.reserve(nspec * nxyz);
    momentum_azimuthal.reserve(4 * nspec * nxyz);
    stress_azimuthal.reserve(4 * nspec * nxyz);

    for (int nid = 0; nid < nxyz; ++nid)
    {
        for (int s = 0; s < nspec; ++s)
        {
            number_density.push_back(nd[idx1(nid, s)]);
            vel_azimuthal.push_back(uth[idx1(nid, s)]);
            ExB_component.push_back(ExB[idx1(nid, s)]);
            diamagnetic_component.push_back(Dia[idx1(nid, s)]);
            inertial_component.push_back(Ine[idx1(nid, s)]);
            stressdiv_component.push_back(Str[idx1(nid, s)]);
            Je_theta_total_fluid_fullmomentum.push_back(ExB[idx1(nid, s)] + Dia[idx1(nid, s)] + Ine[idx1(nid, s)] + Str[idx1(nid, s)]);
            Je_theta_PIC.push_back(Jpic[idx1(nid, s)]);
            momentum_azimuthal.push_back(M[idx4(nid, s) + 0]);
            momentum_azimuthal.push_back(M[idx4(nid, s) + 1]);
            momentum_azimuthal.push_back(M[idx4(nid, s) + 2]);
            momentum_azimuthal.push_back(M[idx4(nid, s) + 3]);
            stress_azimuthal.push_back(S[idx4(nid, s) + 0]);
            stress_azimuthal.push_back(S[idx4(nid, s) + 1]);
            stress_azimuthal.push_back(S[idx4(nid, s) + 2]);
            stress_azimuthal.push_back(S[idx4(nid, s) + 3]);
        }
    }

#undef igpt
#undef idx1
#undef idx4
#undef idx6

    delete[] nd;
    delete[] uth;
    delete[] M;
    delete[] S;
    delete[] p;
    delete[] piT;
    delete[] ExB;
    delete[] Dia;
    delete[] Ine;
    delete[] Str;
    delete[] Jpic;
    delete[] histVelAve;

    return;
}

void SampleField::write_time_averaged_azimuthal(
    const double dt,
    const std::vector<double> &number_density,
    const std::vector<double> &vel_azimuthal,
    const std::vector<double> &momentum_azimuthal,
    const std::vector<double> &stress_azimuthal,
    const std::vector<double> &ExB_component,
    const std::vector<double> &diamagnetic_component,
    const std::vector<double> &inertial_component,
    const std::vector<double> &stressdiv_component,
    const std::vector<double> &Je_theta_total_fluid_fullmomentum,
    const std::vector<double> &Je_theta_PIC)
{
    const std::string prefix = "hist_ave_azimuthalCurrent_";
    auto write_vec = [&](const std::string &fname,
                         const std::vector<double> &vec)
    {
        std::ofstream out(fname.c_str(), std::ofstream::out | std::ofstream::app);

        out << step_now * dt << " "
            << ngpt[0] << " " << ngpt[1] << " " << ngpt[2] << " " << nspec;

        for (double x : vec)
        {
            out << " " << x;
        }

        out << "\n";
        out.close();
        out.clear();
    };
    // Scalars (size = nspec * nxyz)
    write_vec(prefix + "number_density.dat", number_density);
    write_vec(prefix + "vel_theta.dat", vel_azimuthal);

    // 4-comp tensors
    write_vec(prefix + "momentum_tensor.dat", momentum_azimuthal);
    write_vec(prefix + "stress_tensor.dat", stress_azimuthal);

    // Current decomposition scalars
    write_vec(prefix + "ExB_component.dat", ExB_component);
    write_vec(prefix + "diamagnetic_component.dat", diamagnetic_component);
    write_vec(prefix + "inertial_component.dat", inertial_component);
    write_vec(prefix + "stressdiv_component.dat", stressdiv_component);

    // Totals
    write_vec(prefix + "Je_theta_total_fluid.dat", Je_theta_total_fluid_fullmomentum);
    write_vec(prefix + "Je_theta_PIC.dat", Je_theta_PIC);
    return;
}