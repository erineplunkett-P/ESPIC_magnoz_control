#include <iostream>
#include <vector>
#include <string>
#include <sstream>
#include <fstream>
#include <cstdlib>
#include <cmath>
#include <memory.h>

#include "global.h"
#include "pic_io.h"
#include "pic_part.h"
#include "pic_math.h"

#include "Eigen/Sparse"
typedef Eigen::VectorXd Vec;

using namespace std;

const double PI = 2.0 * asin(1.0);
const double PI2 = 4.0 * asin(1.0);

vector<int> infl_outfl_bc[6];
vector<int> infl_outfl_spectype;
vector<double> infl_outfl_weight;
vector<double> infl_outfl_ndens;
vector<double> infl_outfl_T;
vector<double> infl_outfl_V[3];
double infl_outfl_buf[6][6];

int nflows = 0;

double ranf()
{
    double rf;

    rf = ((double)rand() / (double)RAND_MAX);

    return rf;
}

int initFlowInOut(double *partWeight)
{
    int iflow = 0;
    int bc_flag[6];
    int spec_type;
    double ndens, T, V[3];

    ifstream infile("flowInOut.in");
    if (!infile)
        return -1; // return flag (-1 - no in_out_flow injection)
    else
        cout << "Reading [flowInOut.in]..." << endl;

    string line;
    while (getline(infile, line))
    {
        istringstream iss(line);

        iss >> bc_flag[0] >> bc_flag[1] >> bc_flag[2] >> bc_flag[3] >> bc_flag[4] >> bc_flag[5];
        for (int a = 0; a < 6; a++)
            infl_outfl_bc[a].push_back(bc_flag[a]);

        getline(infile, line);
        iss.str(line);
        iss >> spec_type >> ndens >> T >> V[0] >> V[1] >> V[2];
        --spec_type; // in [flowInOut.in], index starts from 1

        infl_outfl_spectype.push_back(spec_type);
        infl_outfl_ndens.push_back(ndens);
        infl_outfl_T.push_back(T);
        infl_outfl_weight.push_back(partWeight[spec_type]);
        for (int a = 0; a < 3; a++)
            infl_outfl_V[a].push_back(V[a]);

        ++iflow;
    }

    if (iflow >= MAXNFLOWS)
    {
        esPIC_exit("# of flow condition exceeds maximum");
    }
    infile.close();
    infile.clear();

    return iflow;
}

// double bcTransM[6][3][3];
double Sn[MAXNFLOWS][6];
double Pnorm[MAXNFLOWS][6];
std::vector<double> fnp_gen[MAXNFLOWS][6]; // (fractional) number of particles needs to be generated
std::vector<double> fnp_res[MAXNFLOWS][6]; // (fractional) residual number of particles

void setInjectPrecomp(double dt, int ngpt[3])
{
    /* compute area at each bc for 3D injection */
    //  for(int a = 0; a < 6; a++) {
    //    bcArea[a] = (domain[(2*(a/2)+3)%6] - domain[(2*(a/2)+2)%6])
    //              * (domain[(2*(a/2)+5)%6] - domain[(2*(a/2)+4)%6]);
    //  }

    for (int iflow = 0; iflow < nflows; iflow++)
    {
        int itype = infl_outfl_spectype.at(iflow);
        double beta = sqrt(0.5 * spec[itype][specmass] / infl_outfl_T.at(iflow));
        for (int ibound = 0; ibound < 6; ibound++)
        {
            if (infl_outfl_bc[ibound].at(iflow) != 1)
                continue;

            double vd[3];
            for (int a = 0; a < 3; a++)
                vd[a] = infl_outfl_V[a].at(iflow);

            /* normal component of drifting vel to the bound */
            /* normalized by thermal speed (2kT/m)^(1/2), different from the buf region method */
            double s = beta * vd[ibound / 2] * (1 - 2 * (ibound % 2));
            Sn[iflow][ibound] = s;

            double h = sqrt(s * s + 2.0);
            double norm = 2.0 / (s + h) * exp(0.5 + 0.5 * s * (s - h));
            Pnorm[iflow][ibound] = norm;

            double vth = sqrt(2.0 * infl_outfl_T.at(iflow) / spec[itype][specmass]);
            double Sproj = Sn[iflow][ibound], nflux(0.0);

            int iside = ibound / 2;
            int isside = (iside + 1) % 2;
            int ncells = ngpt[isside] - 1;
            double A = dx[isside]; // 2d only, if dir is x, then A = dx[y-direction]
                                   //          if dir is y, then A = dx[x-direction]

            if (Sproj >= 0.0)
            {
                nflux = 0.5 * vth * infl_outfl_ndens.at(iflow) * A * dt * (exp(-SQR(Sproj)) / sqrt(PI) + Sproj * (1.0 + erf(Sproj))) / partWeight[itype];
            }
            else
            {
                nflux = 0.5 * vth * infl_outfl_ndens.at(iflow) * A * dt * (exp(-SQR(Sproj)) / sqrt(PI) + Sproj * erfc(-Sproj)) / partWeight[itype];
            }

            fnp_gen[iflow][ibound].resize(ncells, nflux);
            fnp_res[iflow][ibound].resize(ncells, 0.0);

            //      for (int icell = 0; icell < ncells; icell++) {
            //        printf("fnp_gen[%d][%d][%d] = %g\n", iflow, ibound, icell, fnp_gen[iflow][ibound][icell]);
            //      }
        } // end for (int ibound = 0; ibound < 6; ibound++)
    } // end for (int iflow = 0; iflow < nflows; iflow++)
}

void injectFlowInOut(int nspec, long *nparts, int ngpt[3], double *partWeight, double dt, int istep, const Vec &pot)
{
    static int first = 1;

    if (first)
    {
        nflows = initFlowInOut(partWeight);
        if (nflows > 0)
        {
            setInjectPrecomp(dt, ngpt);
            //      srand((int)time(NULL));
        }

        cout << "--> 2D particle injection launched <--" << endl;
        first = 0;
    }

    double part_tmp[MAXNPARTPROPS];
    long npartsTot = nparts[MAXNSPECTYPES];
    for (int iflow = 0; iflow < nflows; iflow++)
    {
        for (int ibound = 0; ibound < 6; ibound++)
        {
            if (infl_outfl_bc[ibound].at(iflow) != 1)
                continue;

            int itype = infl_outfl_spectype.at(iflow);
            if (itype >= nspec)
            {
                ostringstream oss;
                oss << "Inflow(" << iflow + 1 << "): (itype = " << itype + 1 << ") beyonds (nspec = " << nspec << ") read from input";
                esPIC_exit(oss.str());
            };

            double vth = sqrt(2.0 * infl_outfl_T.at(iflow) / spec[itype][specmass]);
            double Sproj = Sn[iflow][ibound];
            int iside = ibound / 2;
            int isside = (iside + 1) % 2;
            int ncells = ngpt[isside] - 1;

            for (int icell = 0; icell < ncells; ++icell)
            {
                int np_gen = (int)(fnp_gen[iflow][ibound][icell] + fnp_res[iflow][ibound][icell] + ranf());
                fnp_res[iflow][ibound][icell] += (fnp_gen[iflow][ibound][icell] - np_gen);

                if (np_gen < 0)
                    continue;

                nparts[itype] += np_gen;
                nparts[MAXNSPECTYPES] += np_gen;
                while (nparts[MAXNSPECTYPES] > MAXNPARTS)
                    adjustPartMem();

                double f, Snew;
                double rf1, rf2;
                double uprim[3];
                int nx = (int)((domain[1] - domain[0]) * dxinv[0] + 0.00001) + 1;

                for (int ipart = 0; ipart < np_gen; ipart++)
                {
                    /* Acceptance-rejection method to sampling drifting injection distribution */
                    while (1)
                    {
                        Snew = -3.0 + 6.0 * ranf();
                        if ((Snew + Sproj) <= 0.0)
                            continue;
                        f = Pnorm[iflow][ibound] * (Snew + Sproj) * exp(-SQR(Snew));
                        if (f > ranf())
                            break;
                    }
                    /*
                              do {
                                Snew = -3.0 + 6.0*ranf();
                                f = Pnorm[iflow][ibound]*(Snew + Sproj)*exp(-SQR(Snew));
                              } while(f < ranf());
                    */
                    uprim[0] = Snew * vth;

                    /* Box-Muller algorithm to genrate particle V from Maxwellian in tangetial direction*/
                    rf1 = vth * sqrt(-log(ranf() + smallValue));
                    rf2 = PI2 * ranf();
                    uprim[1] = rf1 * cos(rf2);
                    uprim[2] = rf1 * sin(rf2);

                    for (int a = 0; a < 3; a++)
                        part_tmp[up + a] = infl_outfl_V[a].at(iflow) + (1 - 2 * (ibound % 2)) * uprim[(a - (ibound / 2) + 3) % 3];

                    part_tmp[xp + iside] = domain[ibound] + dt * part_tmp[up + iside] * ranf();
                    part_tmp[xp + isside] = (icell + ranf()) * dx[isside];
                    part_tmp[xp + 2] = 0.0; // 2D injection, particle_z = 0

                    part_tmp[ptype] = itype;
                    part_tmp[pid] = istep; // group ID for each particle

                    // initial info of a particle
                    part_tmp[xp0] = part_tmp[xp];
                    part_tmp[yp0] = part_tmp[yp];
                    part_tmp[up0] = part_tmp[up];
                    part_tmp[vp0] = part_tmp[vp];

                    // get potenial of nearest nodes
                    double xd[2] = {part_tmp[xp0] - domain[0], part_tmp[yp0] - domain[2]};
                    int ix = (int)(xd[0] * dxinv[0] + 0.5), iy = (int)(xd[1] * dxinv[1] + 0.5);
                    part_tmp[Epot] = spec[itype][specchg] * pot(iy * nx + ix);

                    for (int a = 0; a < partPropN; a++)
                        part[npartsTot][a] = part_tmp[a];

                    ++npartsTot;
                } // end for(int ipart = 0; ipart < np_gen; ipart++)
            } // end for (int icell = 0; icell < ncells; ++icell)
        } // end for(int ibount = 0; ibound < 6; ibound++)
    } // end for(int iflow = 0; iflow < nflows; iflow++)
}
