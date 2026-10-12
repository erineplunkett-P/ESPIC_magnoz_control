#include <cstdio>
#include <cstdlib>
#include <cstddef>
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <time.h>
#include <math.h>

#include "part_spec.h"

#define SQR(a) ((a) * (a))

using namespace std;

const double PI2 = 4.0 * asin(1.0);
const double smallValue = 1e-20;

double dx[3] = {1., 1., 1.};
double dxinv[3];
double domain[6] = {0.0, 1.0, 0.0, 1.0, 0.0, 1.0};

double spec[MAXNSPECTYPES][MAXNSPECPROPS];

double ranf()
{
    double rf;

    rf = ((double)rand() / (double)RAND_MAX);

    return rf;
}

void espic_exit(const string &erro_info)
{
    cerr << erro_info << "...\nExiting..." << endl;
    exit(0);
}

long initSpec()
{
    long nspec = 0;
    /* read in [spec.in] to properties of different species */
    /* hard coded here - electron: type 0; ion: type 1 */

    ifstream infile("specProp.in");
    if (!infile)
        espic_exit("Input file [specProp.in] missing");
    else
        cout << "Reading [specProp.in]..." << endl;

    string line;
    getline(infile, line);
    while (getline(infile, line))
    {
        istringstream iss(line);
        iss >> spec[nspec][specmass] >> spec[nspec][specchg] >> spec[nspec][specFixed];
        spec[nspec][specQM] = spec[nspec][specchg] / spec[nspec][specmass];
        ++nspec;
    }

    return nspec;
}

int main(int argc, char *argv[])
{
    long nspec(0);
    nspec = initSpec(); // read species

    ifstream infile("ambientCtrl.in");
    ofstream outfile("ambient.part", ios::binary);
    ofstream outfile1("weight.part");
    outfile << "# x y z u v w itype\n";
    outfile1 << "# part_weight\n";

    if (!infile)
        return 0; // return generate ambient particles (0 - no mesh generated)
    else
        cout << "Generating ambient partciles according to [ambientCtrl.in]..." << endl;

    int seed;
    {
        string line;
        getline(infile, line);

        istringstream iss(line);
        iss >> seed;
    }
    if (seed < 0)
        srand((int)time(NULL));
    else if (seed > 0)
        srand(seed);

    for (long ispec = 0; ispec < nspec; ispec++)
    {
        string line;
        getline(infile, line);
        istringstream iss(line);
        int nc[3], ngpt[3];
        double domain[6], dx[3];

        iss >> domain[0] >> domain[1] >> domain[2] >> domain[3] >> domain[4] >> domain[5];

        getline(infile, line);
        iss.str(line);
        iss >> dx[0] >> dx[1] >> dx[2];

        for (int a = 0; a < 3; a++)
        {
            nc[a] = (long)((domain[2 * a + 1] - domain[2 * a]) / dx[a] + 1e-5);
            ngpt[a] = nc[a] + 1;
        }

        int npartCell;
        double T, n, vth, U, V, W;

        getline(infile, line);
        iss.str(line);
        iss >> npartCell >> n >> T >> U >> V >> W;
        vth = sqrt(T / spec[ispec][specmass]);

        double dz;
        if (ngpt[2] == 1)
            dz = 1.0;
        else
            dz = dx[2];
        cout << "weight[ispec=" << ispec << "] = " << (double)(n * dx[0] * dx[1] * dz) / npartCell << endl;
        outfile1 << (double)(n * dx[0] * dx[1] * dz) / npartCell << endl;

        if (ngpt[2] == 1)
        {
            if (spec[ispec][specFixed] > 0.0)
            {
                int nparts_unif = (int)sqrt(npartCell);
                int nparts_rest = npartCell - SQR(nparts_unif);
                double pdx = dx[0] / nparts_unif, pdy = dx[1] / nparts_unif;
                double pdx_half = 0.5 * pdx, pdy_half = 0.5 * pdy;

                //        printf("(bool)spec[ispec=%ld][specFixed=%d] = %d, nparts_rest = %d\n",
                //               ispec, specFixed, (bool)spec[ispec][specFixed], nparts_rest);
                printf("spec[ispec=%ld][specFixed=%d] = %g, nparts_rest = %d\n",
                       ispec, specFixed, spec[ispec][specFixed], nparts_rest);

                for (int j = 0; j < nc[1]; j++)
                    for (int i = 0; i < nc[0]; i++)
                    {
                        //            if(i != 0) continue;    // only load particles for first col

                        double xc0, yc0;
                        double x, y, z, u, v, w;
                        double rf1, rf2;

                        z = 0.0;
                        xc0 = i * dx[0] + pdx_half + domain[0];
                        yc0 = j * dx[1] + pdy_half + domain[2];

                        /* uniformly distribute npart_unif^2 particles */
                        for (int ip_y = 0; ip_y != nparts_unif; ip_y++)
                            for (int ip_x = 0; ip_x != nparts_unif; ip_x++)
                            {
                                x = xc0 + ip_x * pdx;
                                y = yc0 + ip_y * pdy;
                                rf1 = vth * sqrt(-2.0 * log(ranf() + smallValue));
                                rf2 = PI2 * ranf();
                                u = rf1 * cos(rf2) + U;
                                v = rf1 * sin(rf2) + V;
                                rf1 = vth * sqrt(-2.0 * log(ranf() + smallValue));
                                rf2 = PI2 * ranf();
                                w = rf1 * cos(rf2) + W;
                                outfile << x << " " << y << " " << z << " "
                                        << u << " " << v << " " << w << " " << ispec << "\n";
                            }
                        /* randomly distribute rest of the particles */
                        for (int ipart = 0; ipart != nparts_rest; ipart++)
                        {
                            x = (xc0 - pdx_half) + dx[0] * ranf();
                            y = (yc0 - pdy_half) + dx[1] * ranf();
                            rf1 = vth * sqrt(-2.0 * log(ranf() + smallValue));
                            rf2 = PI2 * ranf();
                            u = rf1 * cos(rf2) + U;
                            v = rf1 * sin(rf2) + V;
                            rf1 = vth * sqrt(-2.0 * log(ranf() + smallValue));
                            rf2 = PI2 * ranf();
                            w = rf1 * cos(rf2) + W;
                            outfile << x << " " << y << " " << z << " "
                                    << u << " " << v << " " << w << " " << ispec << "\n";
                        } // end for(int ipart = 0; ipart < nparts_rest; ipart++)
                    } // end for(int i = 0; i < nc[0]; i++)
            } // end if(spec[ispec][specFixed] > 0.0)
            else
            {
                for (int j = 0; j < nc[1]; j++)
                    for (int i = 0; i < nc[0]; i++)
                    {
                        //            if(i != 0) continue;    // only load particles for first col

                        double xc0, yc0;
                        double x, y, z, u, v, w;
                        double rf1, rf2;

                        z = 0.0;
                        xc0 = i * dx[0] + domain[0];
                        yc0 = j * dx[1] + domain[2];
                        for (int ipart = 0; ipart < npartCell; ipart++)
                        {
                            x = xc0 + dx[0] * ranf();
                            y = yc0 + dx[1] * ranf();
                            rf1 = vth * sqrt(-2.0 * log(ranf() + smallValue));
                            rf2 = PI2 * ranf();
                            u = rf1 * cos(rf2) + U;
                            v = rf1 * sin(rf2) + V;
                            rf1 = vth * sqrt(-2.0 * log(ranf() + smallValue));
                            rf2 = PI2 * ranf();
                            w = rf1 * cos(rf2) + W;
                            outfile << x << " " << y << " " << z << " "
                                    << u << " " << v << " " << w << " " << ispec << "\n";
                        } // end for(int ipart = 0; ipart < npartCell; ipart++)
                    } // end for(int i = 0; i < nc[0]; i++)
            } // end if(spec[ispec][specFixed] > 0.0) else
        } // end if(ngpt[2] == 1)
        else
        {
            for (int k = 0; k < nc[2]; k++)
                for (int j = 0; j < nc[1]; j++)
                    for (int i = 0; i < nc[0]; i++)
                    {
                        double xc0, yc0, zc0;
                        double x, y, z, u, v, w;
                        double rf1, rf2;

                        xc0 = i * dx[0] + domain[0];
                        yc0 = j * dx[1] + domain[2];
                        zc0 = k * dx[2] + domain[4];
                        for (int ipart = 0; ipart < npartCell; ipart++)
                        {
                            x = xc0 + dx[0] * ranf();
                            y = yc0 + dx[1] * ranf();
                            z = zc0 + dx[2] * ranf();
                            rf1 = vth * sqrt(-2.0 * log(ranf() + smallValue));
                            rf2 = PI2 * ranf();
                            u = rf1 * cos(rf2) + U;
                            v = rf1 * sin(rf2) + V;
                            rf1 = vth * sqrt(-2.0 * log(ranf() + smallValue));
                            rf2 = PI2 * ranf();
                            w = rf1 * cos(rf2) + W;
                            outfile << x << " " << y << " " << z << " "
                                    << u << " " << v << " " << w << " " << ispec << "\n";
                        } // end for(int ipart = 0; ipart < npartCell; ipart++)
                    } // end for(long i = 0; i < nc[0]; i++)
        } // end if(ngpt[2] == 1) else
    } // end for(long ispec = 0; ispec < nspec; ispec++)
    outfile.close();
    outfile.clear();
    outfile1.close();
    outfile1.clear();
    infile.close();
    infile.clear();

    return 0;
}
