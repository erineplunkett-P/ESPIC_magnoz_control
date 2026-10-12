#include "Eigen/Sparse"
typedef Eigen::VectorXd Vec;

#include "beams_in_domain.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include <math.h>

#include "global.h"
#include "pic_io.h"
#include "pic_math.h"
#include "pic_part.h"

Beams::Beams(const std::string &ifilename) : beam(NULL), nbeams(0)
{
    std::ifstream infile(ifilename.c_str());
    if (!infile)
    {
        std::cout << "No beam initialized." << std::endl;
        return;
    }

    const double PI = 2.0 * asin(1.0);
    std::string line;

    // first line -- # of beams
    std::getline(infile, line);
    std::istringstream iss(line);
    iss >> nbeams;

    beam = new OneBeam[nbeams];
    for (int ispec = 0; ispec < MAXNSPECTYPES; ispec++)
        fnparts[ispec] = 0.0;

    size_type ib;
    for (ib = 0; ib < nbeams; ib++)
    {
        // read coordinate of beam injection plane center
        if (!std::getline(infile, line))
            break;
        iss.str(line);
        iss >> beam[ib].xc[0] >> beam[ib].xc[1] >> beam[ib].xc[2];

        // read dimension and radius
        if (!std::getline(infile, line))
            break;
        iss.str(line);
        iss >> beam[ib].dim >> beam[ib].r;

        if (beam[ib].dim == 2)
        {
            beam[ib].area = PI * (SQR(beam[ib].xc[1] + beam[ib].r) - SQR(beam[ib].xc[1] - beam[ib].r));
        }
        else if (beam[ib].dim == 3)
            beam[ib].area = PI * beam[ib].r * beam[ib].r;
        else
        {
            esPIC_exit("Beam dimension must be 2 or 3");
        }

        // read normal vector of injection plane
        if (!std::getline(infile, line))
            break;
        iss.str(line);
        iss >> beam[ib].nvec[0] >> beam[ib].nvec[1] >> beam[ib].nvec[2];

        // read one tangential vector on injection plane
        if (!std::getline(infile, line))
            break;
        iss.str(line);
        iss >> beam[ib].tvec[0][0] >> beam[ib].tvec[0][1] >> beam[ib].tvec[0][2];

        // read beam plasma properties
        if (!std::getline(infile, line))
            break;
        iss.str(line);
        iss >> beam[ib].specid >> beam[ib].ndens >> beam[ib].temp >> beam[ib].vd[0] >> beam[ib].vd[1] >> beam[ib].vd[2];
        beam[ib].specid--;
        beam[ib].wgt = partWeight[beam[ib].specid];
    }

    if (ib != nbeams)
        esPIC_exit("reading beam input file failed");

    inject_init();
}

void Beams::inject_init()
{
    double beta, sn, h;
    double mass, T;
    double nvec_m, tvec_m;

    // pre-compute coefficients used for injecting particles
    for (size_type ib = 0; ib < nbeams; ib++)
    {
        mass = spec[beam[ib].specid][specmass];
        T = beam[ib].temp;
        beta = sqrt(0.5 * mass / T); // most probable speed

        // 2D simulation only in x-y plane
        if (beam[ib].dim == 2)
            beam[ib].nvec[2] = 0.0;

        // normalize n_vector of inj plane
        nvec_m = 0.0;
        for (int a = 0; a < 3; a++)
            nvec_m += (beam[ib].nvec[a] * beam[ib].nvec[a]);
        if (nvec_m == 0.0)
        {
            esPIC_exit("Beam injection plane normal vector cannot be 0");
        }
        nvec_m = sqrt(nvec_m);
        for (int a = 0; a < 3; a++)
            beam[ib].nvec[a] /= nvec_m;

        // velocity projection to normal vector of inj plane
        sn = 0.0;
        for (int a = 0; a < 3; a++)
            sn += (beam[ib].vd[a] * beam[ib].nvec[a]);
        sn *= beta; // speed ratio

        beam[ib].sn = sn;
        h = sqrt(sn * sn + 2.0);
        beam[ib].p_norm = 2.0 / (sn + h) * exp(0.5 + 0.5 * sn * (sn - h));

        // compute tangential vector of injection plane
        if (beam[ib].dim == 2)
        { // compute and update t0 for 2D anyway
            // rotate nvec by 90 deg counter-clockwise
            beam[ib].tvec[0][0] = -beam[ib].nvec[1];
            beam[ib].tvec[0][1] = beam[ib].nvec[0];
            beam[ib].tvec[0][2] = 0.0;

            // t1_vec = nvec cross t0_vec
            // 2D case: t1_vec must be (0, 0, 1) anyway
            beam[ib].tvec[1][0] = beam[ib].nvec[1] * beam[ib].tvec[0][2] - beam[ib].nvec[2] * beam[ib].tvec[0][1];
            beam[ib].tvec[1][1] = beam[ib].nvec[2] * beam[ib].tvec[0][0] - beam[ib].nvec[0] * beam[ib].tvec[0][2];
            beam[ib].tvec[1][2] = beam[ib].nvec[0] * beam[ib].tvec[0][1] - beam[ib].nvec[1] * beam[ib].tvec[0][0];
        }
        else if (beam[ib].dim == 3)
        { // compute t2
            // normalize 1st tangential vector
            tvec_m = 0.0;
            for (int a = 0; a < 3; a++)
                tvec_m += (beam[ib].tvec[0][a] * beam[ib].tvec[0][a]);
            if (tvec_m == 0.0)
            {
                esPIC_exit("Beam injection plane tangential vector cannot be 0");
            }
            tvec_m = sqrt(tvec_m);
            for (int a = 0; a < 3; a++)
                beam[ib].tvec[0][a] /= tvec_m;

            // 2nd tangential vector normal to both n and t0
            // t1_vec = nvec cross t0_vec
            beam[ib].tvec[1][0] = beam[ib].nvec[1] * beam[ib].tvec[0][2] - beam[ib].nvec[2] * beam[ib].tvec[0][1];
            beam[ib].tvec[1][1] = beam[ib].nvec[2] * beam[ib].tvec[0][0] - beam[ib].nvec[0] * beam[ib].tvec[0][2];
            beam[ib].tvec[1][2] = beam[ib].nvec[0] * beam[ib].tvec[0][1] - beam[ib].nvec[1] * beam[ib].tvec[0][0];
        }
    }
}

Beams::Beams(const Beams &source) : beam(NULL), nbeams(source.nbeams)
{
    beam = new OneBeam[nbeams];

    for (int ispec = 0; ispec < MAXNSPECTYPES; ispec++)
        fnparts[ispec] = source.fnparts[ispec];

    for (size_type ib = 0; ib < nbeams; ib++)
    {
        OneBeam *src_beam = &source.beam[ib];

        for (int a = 0; a < 3; a++)
            beam[ib].xc[a] = src_beam->xc[a];
        beam[ib].dim = src_beam->dim;
        beam[ib].r = src_beam->r;
        beam[ib].area = src_beam->area;

        for (int a = 0; a < 3; a++)
        {
            beam[ib].nvec[a] = src_beam->nvec[a];
            beam[ib].tvec[0][a] = src_beam->tvec[0][a];
            beam[ib].tvec[1][a] = src_beam->tvec[1][a];
        }

        beam[ib].specid = src_beam->specid;
        beam[ib].ndens = src_beam->ndens;
        beam[ib].temp = src_beam->temp;

        for (int a = 0; a < 3; a++)
            beam[ib].vd[a] = src_beam->vd[a];

        beam[ib].wgt = src_beam->wgt;
        beam[ib].sn = src_beam->sn;
        beam[ib].p_norm = src_beam->p_norm;
    }
}

Beams::~Beams()
{
    if (nbeams > 0)
    {
        nbeams = 0;
        delete[] beam;
        beam = NULL;
    }
}

// void Beams::inject_particles(double (*part)[MAXNPARTPROPS], long *nparts)
void Beams::inject_particles(double dt, int istep, long *nparts, const Vec &pot)
{
    const double PI = 2.0 * asin(1.0);
    const double PI2 = 4.0 * asin(1.0);
    double part_tmp[MAXNPARTPROPS];
    long nparts_tot = nparts[MAXNSPECTYPES];
    long nparts_gen;
    OneBeam *one;
    int ispec;
    double vth(0.), Sproj(0.), nflux(0.);
    double f, Snew;
    double rf1, rf2;
    double uprim[3], v[3], xd, xs, x[3] = {0.0, 0.0, 0.0};
    double rlosq, rhisq, drsq; // axisymmetric

    for (int ib = 0; ib < nbeams; ib++)
    {
        one = &beam[ib];

        rlosq = SQR(one->xc[1] - one->r); // axisymmetric
        rhisq = SQR(one->xc[1] + one->r);
        drsq = rhisq - rlosq;

        ispec = one->specid;
        vth = sqrt(2.0 * one->temp / spec[ispec][specmass]);
        Sproj = one->sn;
        if (Sproj >= 0.0)
        {
            nflux = 0.5 * dt * vth * one->ndens * one->area * (exp(-SQR(Sproj)) / sqrt(PI) + Sproj * (1.0 + erf(Sproj))) / partWeight[ispec];
        }
        else
        {
            esPIC_exit("Sproj < 0");
        }

        nparts_gen = (long)(nflux + fnparts[ispec]);
        fnparts[ispec] += (nflux - nparts_gen);

        if (nparts_gen <= 0)
            continue;
        nparts[ispec] += nparts_gen;
        nparts[MAXNSPECTYPES] += nparts_gen;

        while (nparts[MAXNSPECTYPES] > MAXNPARTS)
            adjustPartMem();

        int nx = (int)((domain[1] - domain[0]) * dxinv[0] + 0.00001) + 1;

        for (long ipart = 0; ipart < nparts_gen; ipart++)
        {
            /* Acceptance-rejection method to sampling drifting injection distribution */
            while (1)
            {
                Snew = -3.0 + 6.0 * ranf();
                if ((Snew + Sproj) <= 0.0)
                    continue;
                f = one->p_norm * (Snew + Sproj) * exp(-SQR(Snew));
                if (f > ranf())
                    break;
            }
            uprim[0] = Snew * vth;
            xd = (Snew + Sproj) * vth * dt * ranf();

            /* Box-Muller algorithm to genrate particle V from Maxwellian in tangetial direction*/
            rf1 = vth * sqrt(-log(ranf() + smallValue));
            rf2 = PI2 * ranf();
            uprim[1] = rf1 * cos(rf2);
            uprim[2] = rf1 * sin(rf2);

            // project thermal velocity on beam plane coordinate to Cartesian
            v[0] = uprim[0] * one->nvec[0] + uprim[1] * one->tvec[0][0] + uprim[2] * one->tvec[1][0]; // vx
            v[1] = uprim[0] * one->nvec[1] + uprim[1] * one->tvec[0][1] + uprim[2] * one->tvec[1][1]; // vy
            v[2] = uprim[0] * one->nvec[2] + uprim[1] * one->tvec[0][2] + uprim[2] * one->tvec[1][2]; // vz

            // project particle position on beam plane coordinate to Cartesian
            if (one->dim == 2)
            {
                //        xs = 2.0*one->r*(ranf() - 0.5);
                // axisymmetric
                xs = sqrt(rlosq + drsq * ranf());
                x[0] = xd * one->nvec[0] + xs * one->tvec[0][0] + xs * one->tvec[1][0];
                x[1] = xd * one->nvec[1] + xs * one->tvec[0][1] + xs * one->tvec[1][1];
                x[0] += one->xc[0];
            }

            for (int a = 0; a < 3; a++)
            {
                // axisymmetric
                //        part_tmp[xp+a] = x[a] + one->xc[a];
                part_tmp[xp + a] = x[a];
                part_tmp[up + a] = v[a] + one->vd[a];
            }
            part_tmp[ptype] = ispec;
            part_tmp[pid] = istep;

            // initial info of a particle
            part_tmp[xp0] = part_tmp[xp];
            part_tmp[yp0] = part_tmp[yp];
            part_tmp[up0] = part_tmp[up];
            part_tmp[vp0] = part_tmp[vp];

            // get potenial of nearest nodes
            double xd[2] = {part_tmp[xp0] - domain[0], part_tmp[yp0] - domain[2]};
            int ix = (int)(xd[0] * dxinv[0] + 0.5), iy = (int)(xd[1] * dxinv[1] + 0.5);
            part_tmp[Epot] = spec[ispec][specchg] * pot(iy * nx + ix);

            for (int a = 0; a < partPropN; a++)
                part[nparts_tot][a] = part_tmp[a];

            ++nparts_tot;
        } // end for (long ipart = 0; ipart < nparts_gen; ipart++)
    } // end for (int ib = 0; ib < nbeams; ib++)
}
