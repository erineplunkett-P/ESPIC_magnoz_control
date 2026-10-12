#include <iostream>
#include <fstream>
#include <sstream>
#include <cstdio>
#include <cstdlib>
#include <sys/stat.h>

#include "pic_io.h"
#include "vdf.h"

using std::cout;
using std::endl;
using std::ifstream;
using std::istringstream;
using std::ofstream;
using std::ostringstream;
using std::string;

// pulic methods
//
/* -----------------------------------------------------------------------*/
VelDist::VelDist(const string &ifilename,
                 int step,
                 int num_spec,
                 int ngpt[3])
    : step_now(step),
      nspec(num_spec),
      ncells(0),
      ncells_sample(0)
{
    dim_cell[0] = ngpt[0] - 1;
    dim_cell[1] = ngpt[1] - 1;
    if (ngpt[2] == 1)
        dim_cell[2] = 1; // 2d
    else
        dim_cell[2] = ngpt[2] - 1;
    ncells = dim_cell[0] * dim_cell[1] * dim_cell[2];

    init(ifilename);
    read_restart();
}

/* -----------------------------------------------------------------------*/
VelDist::~VelDist()
{
    delete[] cell;
    for (int icell = 0; icell < ncells_sample; ++icell)
    {
        delete[] ncount[icell];
        delete[] hist_array[icell];
    }
    delete[] hist_array;
    delete[] ncount;
}

void VelDist::sample(int step, long nptot, double (*particles)[MAXNPARTPROPS])
{
    // if (step_start < 0)
    // sample for unsteady state;
    // sample starts immediately, ignoring step_start
    if (step_start >= 0)
        if (step <= step_start)
            return; // sampling has NOT started

    if (step % nevery_sample == 0)
    {
        step_now = step;

        // sampling for unsteady case
        if (step_start < 0)
            clear_vdf();

        double boxlo[3] = {domain[0], domain[2], domain[4]};
        double vel, vb0, vb1, dvb;
        int nx = dim_cell[0], ny = dim_cell[1];
        int ix, iy, cellid, histid, ibin;

#define ID(i, j) ((j) * nx + (i))

        for (long ip = 0; ip < nptot; ++ip)
        {
            int ispec = (int)(particles[ip][ptype] + 1e-8);
            if (!spec_exist[ispec])
                continue;

            double r[3] = {particles[ip][xp], particles[ip][yp], particles[ip][zp]};
            ix = (int)((r[0] - boxlo[0]) * dxinv[0]);
            iy = (int)((r[1] - boxlo[1]) * dxinv[1]);

            ix = (ix >= nx ? nx - 1 : ix);
            iy = (iy >= ny ? ny - 1 : iy);
            cellid = cell[ID(ix, iy)];
            if (cellid < 0)
                continue;

            double v[3] = {particles[ip][up], particles[ip][vp], particles[ip][wp]};

            for (auto it_vid = ivel[ispec].cbegin(); it_vid != ivel[ispec].cend(); ++it_vid)
            {
                vel = v[*it_vid];
                histid = map_histid[ispec][*it_vid];
                dvb = dv[histid];
                vb0 = vmin[histid] - 0.5 * dvb;
                vb1 = vmax[histid] + 0.5 * dvb;
                ++ncount[cellid][histid];
                if (vel <= vb0)
                    ++hist_array[cellid][histid][0];
                else if (vel >= vb1)
                    ++hist_array[cellid][histid][nbin[histid] - 1];
                else
                {
                    ibin = (int)((vel - vb0) / dvb);
                    ++hist_array[cellid][histid][ibin];
                }
            }
        }

#undef ID
    }

    if (step % nevery_restart == 0)
        write_restart();

    if (step % nevery_save == 0)
        save_vdf();

    return;
}

// private methods
//
/* -----------------------------------------------------------------------*/
void VelDist::init(const string &ifilename)
{
    ifstream infile(ifilename.c_str());
    if (!infile)
    {
        cout << "Do NOT compute vdf." << endl;
        return;
    }
    cout << "Reading [" << ifilename << "]..." << endl;

    // first line
    string line;
    std::getline(infile, line);
    istringstream iss(line);
    iss >> step_start;

    std::getline(infile, line);
    iss.str(line);
    iss >> nevery_sample;

    std::getline(infile, line);
    iss.str(line);
    iss >> nevery_restart;

    std::getline(infile, line);
    iss.str(line);
    iss >> nevery_save;

    std::getline(infile, line);
    iss.str(line);

    int nblocks;
    iss >> nblocks;

    int(*clo)[3] = new int[nblocks][3];
    int(*chi)[3] = new int[nblocks][3];
    for (int i = 0; i < nblocks; ++i)
    {
        iss >> clo[i][0] >> clo[i][1] >> clo[i][2] >> chi[i][0] >> chi[i][1] >> chi[i][2];

        for (int k = 0; k < 3; ++k)
        {
            if (clo[i][k] < 0 || clo[i][k] >= dim_cell[k])
                esPIC_exit("Invalid parameters in vdf input");
            if (chi[i][k] < 0 || chi[i][k] >= dim_cell[k])
                esPIC_exit("Invalid parameters in vdf input");
        }
    }

    int nx = dim_cell[0], ny = dim_cell[1];

#define ID(i, j, k) (((k) * (ny) + (j)) * (nx) + (i))

    cell = new int[ncells];
    for (int icell = 0; icell < ncells; ++icell)
        cell[icell] = -1;

    for (int ib = 0; ib < nblocks; ++ib)
    {
        for (int k = clo[ib][2]; k <= chi[ib][2]; ++k)
            for (int j = clo[ib][1]; j <= chi[ib][1]; ++j)
            {
                for (int i = clo[ib][0]; i <= chi[ib][0]; ++i)
                {
                    cell[ID(i, j, k)] = ncells_sample++;
                }
            }
    }

#undef ID

    delete[] clo;
    delete[] chi;

    spec_exist.resize(nspec, false);
    for (int i = 0; i < MAXNSPECTYPES; ++i)
        for (int j = 0; j < 3; ++j)
            map_histid[i][j] = -1;

    int histid = 0;
    while (std::getline(infile, line))
    {
        iss.str(line);

        int ispec;
        char d;
        double v[3];

        iss >> ispec >> d >> v[0] >> v[1] >> v[2];
        --ispec;

        if (ispec > nspec)
            esPIC_exit("Invalid spec ID in vdf input");

        spec_exist[ispec] = true;

        vmin.push_back(v[0]);
        vmax.push_back(v[1]);
        dv.push_back(v[2]);

        if (d == 'x')
            ivel[ispec].push_back(0);
        else if (d == 'y')
            ivel[ispec].push_back(1);
        else if (d == 'z')
            ivel[ispec].push_back(2);
        else
            esPIC_exit("Invalid velocity component in vdf input");

        map_histid[ispec][ivel[ispec].back()] = histid++;
    }

    int nhist = histid;

    nbin.resize(nhist, 0);
    for (histid = 0; histid < nhist; ++histid)
        nbin[histid] = (int)((vmax[histid] - vmin[histid]) / dv[histid]) + 1;

    ncount = new int *[ncells_sample];
    hist_array = new vector<int> *[ncells_sample];
    for (int icell = 0; icell < ncells_sample; ++icell)
    {
        ncount[icell] = new int[nhist];
        hist_array[icell] = new vector<int>[nhist];

        for (histid = 0; histid < nhist; ++histid)
        {
            hist_array[icell][histid].resize(nbin[histid], 0);
            ncount[icell][histid] = 0;
        }
    }
}

/* -----------------------------------------------------------------------*/
void VelDist::write_restart()
{
    FILE *ofp;
    ofp = fopen("esPICVDF.restart", "wb");

    fwrite(&step_now, sizeof(int), 1, ofp);
    fwrite(&nspec, sizeof(int), 1, ofp);
    for (auto it = spec_exist.cbegin(); it != spec_exist.cend(); ++it)
    {
        bool flag = *it;
        fwrite(&flag, sizeof(bool), 1, ofp);
    }

    fwrite(&dim_cell[0], sizeof(int), 3, ofp);
    fwrite(&ncells_sample, sizeof(int), 1, ofp);
    fwrite(&cell[0], sizeof(int), ncells, ofp);

    int maxnspec = MAXNSPECTYPES;
    int nhist, ispec;
    fwrite(&maxnspec, sizeof(int), 1, ofp);
    for (ispec = 0; ispec < maxnspec; ++ispec)
    {
        nhist = (int)ivel[ispec].size();
        fwrite(&nhist, sizeof(int), 1, ofp);

        for (int i = 0; i < nhist; ++i)
        {
            int vid = ivel[ispec][i];
            fwrite(&vid, sizeof(int), 1, ofp);
        }
    }

    fwrite(&map_histid[0][0], sizeof(int), 3 * maxnspec, ofp);

    nhist = (int)vmin.size();
    fwrite(&nhist, sizeof(int), 1, ofp);

    for (int i = 0; i < nhist; ++i)
    {
        double v = vmin[i];
        fwrite(&v, sizeof(double), 1, ofp);
    }

    for (int i = 0; i < nhist; ++i)
    {
        double v = vmax[i];
        fwrite(&v, sizeof(double), 1, ofp);
    }

    for (int i = 0; i < nhist; ++i)
    {
        double v = dv[i];
        fwrite(&v, sizeof(double), 1, ofp);
    }

    for (int icell = 0; icell < ncells_sample; ++icell)
    {
        for (int histid = 0; histid < nhist; ++histid)
        {
            for (int ibin = 0; ibin < nbin[histid]; ++ibin)
            {
                int n = hist_array[icell][histid][ibin];
                fwrite(&n, sizeof(int), 1, ofp);
            }
        }
    }

    for (int icell = 0; icell < ncells_sample; ++icell)
        fwrite(&ncount[icell][0], sizeof(int), nhist, ofp);

    fclose(ofp);
    return;
}

/* -----------------------------------------------------------------------*/
void VelDist::read_restart()
{
    if (!is_file_exist("esPICVDF.restart"))
        return; // no vdf data saved

    FILE *ifp;
    ifp = fopen("esPICVDF.restart", "rb");
    if (ifp == NULL)
    {
        esPIC_exit("Failed to open vdf restart file[esPICVDF.restart]");
    }

    long size_read;

    int step_old;
    size_read = fread(&step_old, sizeof(int), 1, ifp);

    bool flag_reset = step_now < step_old;
    if (flag_reset)
    {
        cout << "Warning: reset data saved for vdf.\n";

        fclose(ifp);
        return;
    }

    int nspec_old;
    size_read = fread(&nspec_old, sizeof(int), 1, ifp);
    if (nspec_old != nspec)
    {
        esPIC_exit("Failed to read [esPICVDF.restart]: # of species does not match");
    }

    bool flag_spec[MAXNSPECTYPES];
    size_read = fread(flag_spec, sizeof(bool), nspec, ifp);

    if (size_read != nspec)
        esPIC_exit("Bad restart file [esPICVDF.restart]");

    int ispec;
    for (ispec = 0; ispec < nspec; ++ispec)
    {
        if (flag_spec[ispec] != spec_exist[ispec])
            esPIC_exit("Failed to read [esPICVDF.restart]: species to be processed do not match");
    }

    int dim_old[3];
    size_read = fread(&dim_old[0], sizeof(int), 3, ifp);
    if (size_read != 3)
        esPIC_exit("Bad restart file [esPICVDF.restart]");

    for (int i = 0; i < 3; ++i)
    {
        if (dim_old[i] != dim_cell[i])
            esPIC_exit("Failed to read [esPICVDF.restart]: # of cells does not match");
    }

    int ncells_sample_old;
    size_read = fread(&ncells_sample_old, sizeof(int), 1, ifp);

    if (ncells_sample_old != ncells_sample)
        esPIC_exit("Failed to read [esPICVDF.restart]: # of cells to be processed does not match");

    int *cell_old = new int[ncells];
    size_read = fread(&cell_old[0], sizeof(int), ncells, ifp);

    if (size_read != ncells)
        esPIC_exit("Bad restart file [esPICVDF.restart]");

    for (int icell = 0; icell < ncells; ++icell)
    {
        if (cell_old[icell] != cell[icell])
            esPIC_exit("Failed to read [esPICVDF.restart]: cell to be processed does not match");
    }
    delete[] cell_old;

    int maxnspec;
    int nhist;
    size_read = fread(&maxnspec, sizeof(int), 1, ifp);

    if (maxnspec != MAXNSPECTYPES)
        esPIC_exit("Failed to read [esPICVDF.restart]: max # of species does not match");

    for (ispec = 0; ispec < maxnspec; ++ispec)
    {
        size_read = fread(&nhist, sizeof(int), 1, ifp);

        if (nhist != (int)ivel[ispec].size())
            esPIC_exit("Failed to read [esPICVDF.restart]: # of vdf to be processed does not match");

        for (int i = 0; i < nhist; ++i)
        {
            int vid;
            size_read = fread(&vid, sizeof(int), 1, ifp);
            if (vid != ivel[ispec][i])
                esPIC_exit("Failed to read [esPICVDF.restart]: component of vdf to be processed does not match");
        }
    }

    int map_histid_old[MAXNSPECTYPES][3];
    size_read = fread(&map_histid_old[0][0], sizeof(int), 3 * maxnspec, ifp);

    if (size_read != 3 * maxnspec)
        esPIC_exit("Bad restart file [esPICVDF.restart]");

    for (ispec = 0; ispec < maxnspec; ++ispec)
        for (int i = 0; i < 3; ++i)
        {
            if (map_histid_old[ispec][i] != map_histid[ispec][i])
                esPIC_exit("Failed to read [esPICVDF.restart]: order of component of vdf to be processed does not match");
        }

    size_read = fread(&nhist, sizeof(int), 1, ifp);
    if (nhist != (int)vmin.size())
        esPIC_exit("Failed to read [esPICVDF.restart]: # of vdf to be processed does not match");

    double *v = new double[nhist];

    size_read = fread(&v[0], sizeof(double), nhist, ifp);
    if (size_read != nhist)
        esPIC_exit("Bad restart file [esPICVDF.restart]");

    for (int i = 0; i < nhist; ++i)
    {
        if (v[i] != vmin[i])
            esPIC_exit("Failed to read [esPICVDF.restart]: vmin does not match");
    }

    size_read = fread(&v[0], sizeof(double), nhist, ifp);
    if (size_read != nhist)
        esPIC_exit("Bad restart file [esPICVDF.restart]");

    for (int i = 0; i < nhist; ++i)
    {
        if (v[i] != vmax[i])
            esPIC_exit("Failed to read [esPICVDF.restart]: vmax does not match");
    }

    size_read = fread(&v[0], sizeof(double), nhist, ifp);
    if (size_read != nhist)
        esPIC_exit("Bad restart file [esPICVDF.restart]");

    for (int i = 0; i < nhist; ++i)
    {
        if (v[i] != dv[i])
            esPIC_exit("Failed to read [esPICVDF.restart]: vmax does not match");
    }

    delete[] v;

    int *nbin_old = new int[nhist];
    int maxnbin = 0;
    for (int histid = 0; histid < nhist; ++histid)
    {
        nbin_old[histid] = (int)((vmax[histid] - vmin[histid]) / dv[histid]) + 1;
        maxnbin = (maxnbin > nbin_old[histid] ? maxnbin : nbin_old[histid]);
    }
    int *hist = new int[maxnbin];

    for (int icell = 0; icell < ncells_sample; ++icell)
    {
        for (int histid = 0; histid < nhist; ++histid)
        {

            size_read = fread(&hist[0], sizeof(int), nbin_old[histid], ifp);
            if (size_read != nbin_old[histid])
                esPIC_exit("Failed to read [esPICVDF.restart]: failed to read data from hist_array");

            for (int ibin = 0; ibin < nbin_old[histid]; ++ibin)
                hist_array[icell][histid][ibin] = hist[ibin];
        }
    }

    for (int icell = 0; icell < ncells_sample; ++icell)
    {
        size_read = fread(&ncount[icell][0], sizeof(int), nhist, ifp);
        if (size_read != nhist)
            esPIC_exit("Bad restart file [esPICVDF.restart]");
    }

    delete[] hist;
    delete[] nbin_old;

    fclose(ifp);
    return;
}

void VelDist::save_vdf()
{
    int nx = dim_cell[0], ny = dim_cell[1];
    int icell, cellid, histid, ibin;

    for (int ispec = 0; ispec < nspec; ++ispec)
    {
        if (!spec_exist[ispec])
            continue;

        for (auto it_vid = ivel[ispec].cbegin(); it_vid != ivel[ispec].cend(); ++it_vid)
        {
            histid = map_histid[ispec][*it_vid];

            ostringstream oss;
            oss << "vbin_s" << ispec << "_v" << *it_vid << ".dat";

            ofstream outfile;
            outfile.open(oss.str().c_str(), ofstream ::out);

            for (ibin = 0; ibin + 1 < nbin[histid]; ++ibin)
                outfile << vmin[histid] + ibin * dv[histid] << "\n";
            outfile << vmax[histid];
            outfile.close();

            oss.str("");
            if (step_start >= 0) // steady-state sampling
                oss << "vdf_hist_s" << ispec << "_v" << *it_vid << ".dat";
            else
                oss << "vdf_hist_s" << ispec << "_v" << *it_vid << "_step" << step_now << ".dat";

            outfile.clear();
            outfile.open(oss.str().c_str(), ofstream::out);
            icell = 0;
            for (int iy = 0; iy < ny; ++iy)
            {
                for (int ix = 0; ix < nx; ++ix)
                {
                    cellid = cell[icell];
                    if (cellid >= 0)
                    {
                        outfile << ix << " " << iy << " " << ncount[cellid][histid];
                        for (ibin = 0; ibin < nbin[histid]; ++ibin)
                            outfile << " " << hist_array[cellid][histid][ibin];
                        outfile << "\n";
                    }
                    ++icell;
                }
            }
            outfile.close();
        }
    }

    return;
}

void VelDist::clear_vdf()
{
    auto nhist = nbin.size();
    decltype(nhist) histid = 0;

    for (int icell = 0; icell < ncells_sample; ++icell)
    {
        for (histid = 0; histid < nhist; ++histid)
        {
            hist_array[icell][histid].assign(nbin[histid], 0);
            ncount[icell][histid] = 0;
        }
    }

    return;
}
