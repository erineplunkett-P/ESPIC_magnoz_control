#include <iostream>
#include <fstream>
#include <sstream>
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <sys/stat.h>

#include "pic_io.h"
#include "edf.h"

using std::cout;
using std::endl;
using std::ifstream;
using std::istringstream;
using std::ofstream;
using std::ostringstream;
using std::string;

// pulic methods
EnergyDist::EnergyDist(const string &ifilename,
                       int step,
                       int num_spec,
                       int ngpt[3])
    : step_now(step),
      nspec(num_spec),
      ncells(0),
      ncells_sample(0),
      cell(nullptr),
      hist_array_full(nullptr),
      hist_array_axis(nullptr),
      ncount(nullptr)
{
    dim_cell[0] = ngpt[0] - 1;
    dim_cell[1] = ngpt[1] - 1;
    if (ngpt[2] == 1)
    {
        dim_cell[2] = 1;
    }
    else
    {
        dim_cell[2] = ngpt[2] - 1;
    }
    ncells = dim_cell[0] * dim_cell[1] * dim_cell[2];

    init(ifilename);
    read_restart();
}

EnergyDist::~EnergyDist()
{
    delete[] cell;

    for (int icell = 0; icell < ncells_sample; ++icell)
    {
        delete[] ncount[icell];
        delete[] hist_array_full[icell];
        delete[] hist_array_axis[icell];
    }

    delete[] hist_array_full;
    delete[] hist_array_axis;
    delete[] ncount;
}

void EnergyDist::sample(int step, long nptot, double (*particles)[MAXNPARTPROPS])
{
    // step_start < 0 for unsteady sampling: clear every sample
    if (step_start >= 0)
    {
        if (step <= step_start)
        {
            return;
        }
    }
    if (step % nevery_sample == 0)
    {
        step_now = step;
        if (step_start < 0)
        {
            clear_edf();
        }
        double boxlo[3] = {domain[0], domain[2], domain[4]};
        int nx = dim_cell[0], ny = dim_cell[1];
        int ix, iy, cellid, histid;
#define ID(i, j) ((j) * nx + (i))
        for (long ip = 0; ip < nptot; ++ip)
        {
            int ispec = (int)(particles[ip][ptype] + 1e-8);
            if (ispec < 0 || ispec >= nspec)
            {
                continue;
            }
            if (!spec_exist[ispec])
            {
                continue;
            }
            double r[3] = {particles[ip][xp], particles[ip][yp], particles[ip][zp]};
            ix = (int)((r[0] - boxlo[0]) * dxinv[0]);
            iy = (int)((r[1] - boxlo[1]) * dxinv[1]);
            ix = (ix >= nx ? nx - 1 : ix);
            iy = (iy >= ny ? ny - 1 : iy);
            ix = (ix < 0 ? 0 : ix);
            iy = (iy < 0 ? 0 : iy);
            cellid = cell[ID(ix, iy)];
            if (cellid < 0)
            {
                continue;
            }
            histid = map_histid[ispec];
            if (histid < 0)
            {
                continue;
            }
            double v[3] = {particles[ip][up], particles[ip][vp], particles[ip][wp]};
            const double m = spec[ispec][specmass];
            const int mode = emode[histid];
            const bool sample_full = (mode == ENERGY_FULL3 || mode == ENERGY_BOTH);
            const bool sample_axis = (mode == ENERGY_AXIS_ZR || mode == ENERGY_BOTH);
            const double v2_planar = v[0] * v[0] + v[1] * v[1];
            const double ek_axis = 0.5 * m * v2_planar;
            const double ek_full = (sample_full ? 0.5 * m * (v2_planar + v[2] * v[2]) : 0.0);
            const double deb = de[histid];
            const double eb0 = emin[histid] - 0.5 * deb;
            const double eb1 = emax[histid] + 0.5 * deb;
            ++ncount[cellid][histid];
            auto tally_energy = [&](double ek, vector<int> &bins) {
                if (ek <= eb0)
                {
                    ++bins[0];
                    return;
                }
                if (ek >= eb1)
                {
                    ++bins[nbin[histid] - 1];
                    return;
                }
                int ibin_local = (int)((ek - eb0) / deb);
                if (ibin_local < 0)
                {
                    ibin_local = 0;
                }
                if (ibin_local >= nbin[histid])
                {
                    ibin_local = nbin[histid] - 1;
                }
                ++bins[ibin_local];
            };
            if (sample_axis)
            {
                tally_energy(ek_axis, hist_array_axis[cellid][histid]);
            }
            if (sample_full)
            {
                tally_energy(ek_full, hist_array_full[cellid][histid]);
            }
        }
#undef ID
    }
    if (step % nevery_restart == 0)
    {
        write_restart();
    }
    if (step % nevery_save == 0)
    {
        save_edf_epf();
    }
    return;
}

// private methods
void EnergyDist::init(const string &ifilename)
{
    ifstream infile(ifilename.c_str());
    if (!infile)
    {
        cout << "Do NOT compute edf/epf." << endl;
        return;
    }

    cout << "Reading [" << ifilename << "]..." << endl;

    string line;
    istringstream iss;

    // step_start
    std::getline(infile, line);
    iss.str(line);
    iss >> step_start;
    // nevery_sample
    std::getline(infile, line);
    iss.clear();
    iss.str(line);
    iss >> nevery_sample;
    // nevery_restart
    std::getline(infile, line);
    iss.clear();
    iss.str(line);
    iss >> nevery_restart;
    // nevery_save
    std::getline(infile, line);
    iss.clear();
    iss.str(line);
    iss >> nevery_save;
    // blocks
    std::getline(infile, line);
    iss.clear();
    iss.str(line);
    int nblocks;
    iss >> nblocks;
    int (*clo)[3] = new int[nblocks][3];
    int (*chi)[3] = new int[nblocks][3];
    for (int i = 0; i < nblocks; ++i)
    {
        iss >> clo[i][0] >> clo[i][1] >> clo[i][2] >> chi[i][0] >> chi[i][1] >> chi[i][2];
        for (int k = 0; k < 3; ++k)
        {
            if (clo[i][k] < 0 || clo[i][k] >= dim_cell[k])
            {
                esPIC_exit("Invalid parameters in edf input");
            }
            if (chi[i][k] < 0 || chi[i][k] >= dim_cell[k])
            {
                esPIC_exit("Invalid parameters in edf input");
            }
        }
    }
    int nx = dim_cell[0], ny = dim_cell[1]; //, nz = dim_cell[2];
#define ID(i, j, k) (((k) * (ny) + (j)) * (nx) + (i))
    cell = new int[ncells];
    for (int icell = 0; icell < ncells; ++icell)
    {
        cell[icell] = -1;
    }
    for (int ib = 0; ib < nblocks; ++ib)
    {
        for (int k = clo[ib][2]; k <= chi[ib][2]; ++k)
        {
            for (int j = clo[ib][1]; j <= chi[ib][1]; ++j)
            {
                for (int i = clo[ib][0]; i <= chi[ib][0]; ++i)
                {
                    cell[ID(i, j, k)] = ncells_sample++;
                }
            }
        }
    }
#undef ID
    delete[] clo;
    delete[] chi;
    spec_exist.resize(nspec, false);
    for (int i = 0; i < MAXNSPECTYPES; ++i)
    {
        map_histid[i] = -1;
    }
    // ispec emin emax de
    int histid = 0;
    while (std::getline(infile, line))
    {
        if (line.size() == 0)
        {
            continue;
        }
        iss.clear();
        iss.str(line);
        int ispec;
        double e0, e1, de_in;
        std::string mode_str;
        iss >> ispec >> e0 >> e1 >> de_in;
        --ispec;
        if (ispec < 0 || ispec >= nspec)
        {
            esPIC_exit("Invalid spec ID in edf input");
        }
        if (de_in <= 0.0 || e1 <= e0)
        {
            esPIC_exit("Invalid energy bin parameters in edf input");
        }
        int mode = ENERGY_FULL3; // default
        if (iss >> mode_str)
        {
            for (auto &ch : mode_str)
                if (ch >= 'A' && ch <= 'Z')
                    ch = char(ch - 'A' + 'a');

            if (mode_str == "zr" || mode_str == "axis" || mode_str == "2d" ||
                mode_str == "no_theta" || mode_str == "a")
            {
                mode = ENERGY_AXIS_ZR;
            }
            else if (mode_str == "full" || mode_str == "3d" || mode_str == "xyz" ||
                     mode_str == "all" || mode_str == "f")
            {
                mode = ENERGY_FULL3;
            }
            else if (mode_str == "both" || mode_str == "dual" || mode_str == "fa")
            {
                mode = ENERGY_BOTH;
            }
            else
            {
                esPIC_exit("Invalid energy mode in edf input (use full/3d/xyz/all, zr/axis/2d, or both)");
            }
        }
        spec_exist[ispec] = true;
        emin.push_back(e0);
        emax.push_back(e1);
        de.push_back(de_in);
        emode.push_back(mode);
        map_histid[ispec] = histid++;
    }
    const int nhist = histid;
    nbin.resize(nhist, 0);
    for (int h = 0; h < nhist; ++h)
    {
        nbin[h] = (int)((emax[h] - emin[h]) / de[h]) + 1;
    }
    ncount = new int *[ncells_sample];
    hist_array_full = new vector<int> *[ncells_sample];
    hist_array_axis = new vector<int> *[ncells_sample];
    for (int icell = 0; icell < ncells_sample; ++icell)
    {
        ncount[icell] = new int[nhist];
        hist_array_full[icell] = new vector<int>[nhist];
        hist_array_axis[icell] = new vector<int>[nhist];
        for (int h = 0; h < nhist; ++h)
        {
            hist_array_full[icell][h].assign(nbin[h], 0);
            hist_array_axis[icell][h].assign(nbin[h], 0);
            ncount[icell][h] = 0;
        }
    }
    return;
}

/* -----------------------------------------------------------------------*/
void EnergyDist::write_restart()
{
    FILE *ofp = fopen("esPICEFDF.restart", "wb");
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
    fwrite(&maxnspec, sizeof(int), 1, ofp);
    fwrite(&map_histid[0], sizeof(int), maxnspec, ofp);
    int nhist = (int)emin.size();
    fwrite(&nhist, sizeof(int), 1, ofp);
    for (int i = 0; i < nhist; ++i)
    {
        fwrite(&emin[i], sizeof(double), 1, ofp);
    }
    for (int i = 0; i < nhist; ++i)
    {
        fwrite(&emax[i], sizeof(double), 1, ofp);
    }
    for (int i = 0; i < nhist; ++i)
    {
        fwrite(&de[i], sizeof(double), 1, ofp);
    }
    for (int i = 0; i < nhist; ++i)
    {
        int m = emode[i];
        fwrite(&m, sizeof(int), 1, ofp);
    }
    for (int icell = 0; icell < ncells_sample; ++icell)
    {
        for (int h = 0; h < nhist; ++h)
        {
            auto write_bins = [&](const vector<int> &bins) {
                for (int ibin = 0; ibin < nbin[h]; ++ibin)
                {
                    int n = bins[ibin];
                    fwrite(&n, sizeof(int), 1, ofp);
                }
            };
            const int mode = emode[h];
            if (mode == ENERGY_AXIS_ZR)
            {
                write_bins(hist_array_axis[icell][h]);
            }
            else
            {
                write_bins(hist_array_full[icell][h]);
            }
            if (mode == ENERGY_BOTH)
            {
                write_bins(hist_array_axis[icell][h]);
            }
        }
    }
    for (int icell = 0; icell < ncells_sample; ++icell)
    {
        fwrite(&ncount[icell][0], sizeof(int), nhist, ofp);
    }
    fclose(ofp);
    return;
}

/* -----------------------------------------------------------------------*/
void EnergyDist::read_restart()
{
    if (!is_file_exist("esPICEFDF.restart"))
    {
        return;
    }
    FILE *ifp = fopen("esPICEFDF.restart", "rb");
    if (ifp == NULL)
    {
        esPIC_exit("Failed to open edf restart file[esPICEFDF.restart]");
    }
    long size_read;
    int step_old;
    size_read = fread(&step_old, sizeof(int), 1, ifp);
    bool flag_reset = step_now < step_old;
    if (flag_reset)
    {
        cout << "Warning: reset data saved for edf.\n";
        fclose(ifp);
        return;
    }
    int nspec_old;
    size_read = fread(&nspec_old, sizeof(int), 1, ifp);
    if (nspec_old != nspec)
    {
        esPIC_exit("Failed to read [esPICEFDF.restart]: # of species does not match");
    }
    bool flag_spec[MAXNSPECTYPES];
    size_read = fread(flag_spec, sizeof(bool), nspec, ifp);
    if (size_read != nspec)
    {
        esPIC_exit("Bad restart file [esPICEFDF.restart]");
    }
    for (int ispec = 0; ispec < nspec; ++ispec)
    {
        if (flag_spec[ispec] != spec_exist[ispec])
        {
            esPIC_exit("Failed to read [esPICEFDF.restart]: species to be processed do not match");
        }
    }
    int dim_old[3];
    size_read = fread(&dim_old[0], sizeof(int), 3, ifp);
    if (size_read != 3)
        esPIC_exit("Bad restart file [esPICEFDF.restart]");

    for (int i = 0; i < 3; ++i)
    {
        if (dim_old[i] != dim_cell[i])
            esPIC_exit("Failed to read [esPICEFDF.restart]: # of cells does not match");
    }

    int ncells_sample_old;
    size_read = fread(&ncells_sample_old, sizeof(int), 1, ifp);
    if (ncells_sample_old != ncells_sample)
    {
        esPIC_exit("Failed to read [esPICEFDF.restart]: # of cells to be processed does not match");
    }
    int *cell_old = new int[ncells];
    size_read = fread(&cell_old[0], sizeof(int), ncells, ifp);
    if (size_read != ncells)
    {
        esPIC_exit("Bad restart file [esPICEFDF.restart]");
    }
    for (int icell = 0; icell < ncells; ++icell)
    {
        if (cell_old[icell] != cell[icell])
        {
            esPIC_exit("Failed to read [esPICEFDF.restart]: cell mapping does not match");
        }
    }
    delete[] cell_old;
    int maxnspec;
    size_read = fread(&maxnspec, sizeof(int), 1, ifp);
    if (maxnspec != MAXNSPECTYPES)
    {
        esPIC_exit("Failed to read [esPICEFDF.restart]: max # of species does not match");
    }
    int map_histid_old[MAXNSPECTYPES];
    size_read = fread(&map_histid_old[0], sizeof(int), maxnspec, ifp);
    if (size_read != maxnspec)
    {
        esPIC_exit("Bad restart file [esPICEFDF.restart]");
    }
    for (int ispec = 0; ispec < maxnspec; ++ispec)
    {
        if (map_histid_old[ispec] != map_histid[ispec])
        {
            esPIC_exit("Failed to read [esPICEFDF.restart]: map_histid mismatch");
        }
    }
    int nhist;
    size_read = fread(&nhist, sizeof(int), 1, ifp);
    if (nhist != (int)emin.size())
    {
        esPIC_exit("Failed to read [esPICEFDF.restart]: # of edf to be processed does not match");
    }
    double *tmp = new double[nhist];
    size_read = fread(&tmp[0], sizeof(double), nhist, ifp);
    if (size_read != nhist)
    {
        esPIC_exit("Bad restart file [esPICEFDF.restart]");
    }
    for (int i = 0; i < nhist; ++i)
    {
        if (tmp[i] != emin[i])
        {
            esPIC_exit("Failed to read [esPICEFDF.restart]: emin does not match");
        }
    }
    size_read = fread(&tmp[0], sizeof(double), nhist, ifp);
    if (size_read != nhist)
    {
        esPIC_exit("Bad restart file [esPICEFDF.restart]");
    }
    for (int i = 0; i < nhist; ++i)
    {
        if (tmp[i] != emax[i])
        {
            esPIC_exit("Failed to read [esPICEFDF.restart]: emax does not match");
        }
    }
    size_read = fread(&tmp[0], sizeof(double), nhist, ifp);
    if (size_read != nhist)
    {
        esPIC_exit("Bad restart file [esPICEFDF.restart]");
    }
    for (int i = 0; i < nhist; ++i)
    {
        if (tmp[i] != de[i])
        {
            esPIC_exit("Failed to read [esPICEFDF.restart]: de does not match");
        }
    }
    delete[] tmp;
    // read emode and verify
    int *mode_old = new int[nhist];
    size_read = fread(&mode_old[0], sizeof(int), nhist, ifp);
    if (size_read != nhist)
        esPIC_exit("Bad restart file [esPICEFDF.restart] (emode)");

    for (int i = 0; i < nhist; ++i)
    {
        if (mode_old[i] != emode[i])
            esPIC_exit("Failed to read [esPICEFDF.restart]: energy mode does not match");
    }
    delete[] mode_old;
    // Read histograms
    int *nbin_old = new int[nhist];
    int maxnbin = 0;
    for (int h = 0; h < nhist; ++h)
    {
        nbin_old[h] = (int)((emax[h] - emin[h]) / de[h]) + 1;
        maxnbin = (maxnbin > nbin_old[h] ? maxnbin : nbin_old[h]);
    }
    int *hist = new int[maxnbin];
    for (int icell = 0; icell < ncells_sample; ++icell)
    {
        for (int h = 0; h < nhist; ++h)
        {
            size_read = fread(&hist[0], sizeof(int), nbin_old[h], ifp);
            if (size_read != nbin_old[h])
            {
                esPIC_exit("Failed to read [esPICEFDF.restart]: hist_array data");
            }
            for (int ibin = 0; ibin < nbin_old[h]; ++ibin)
            {
                if (emode[h] == ENERGY_AXIS_ZR)
                {
                    hist_array_axis[icell][h][ibin] = hist[ibin];
                }
                else
                {
                    hist_array_full[icell][h][ibin] = hist[ibin];
                }
            }
            if (emode[h] == ENERGY_BOTH)
            {
                size_read = fread(&hist[0], sizeof(int), nbin_old[h], ifp);
                if (size_read != nbin_old[h])
                {
                    esPIC_exit("Failed to read [esPICEFDF.restart]: hist_array data (axis)");
                }
                for (int ibin = 0; ibin < nbin_old[h]; ++ibin)
                {
                    hist_array_axis[icell][h][ibin] = hist[ibin];
                }
            }
        }
    }
    for (int icell = 0; icell < ncells_sample; ++icell)
    {
        size_read = fread(&ncount[icell][0], sizeof(int), nhist, ifp);
        if (size_read != nhist)
        {
            esPIC_exit("Bad restart file [esPICEFDF.restart]");
        }
    }
    delete[] hist;
    delete[] nbin_old;
    fclose(ifp);
    return;
}

void EnergyDist::save_edf_epf()
{
    int nx = dim_cell[0], ny = dim_cell[1];
    int icell, cellid, histid, ibin;
    const double epsE = 1e-30;
    for (int ispec = 0; ispec < nspec; ++ispec)
    {
        if (!spec_exist[ispec])
        {
            continue;
        }
        histid = map_histid[ispec];
        if (histid < 0)
        {
            continue;
        }

        struct HistOutput
        {
            vector<int> **bins;
            std::string label;
        };
        vector<HistOutput> outputs;
        const int mode = emode[histid];
        if (mode == ENERGY_AXIS_ZR)
        {
            outputs.push_back({hist_array_axis, ""});
        }
        else if (mode == ENERGY_FULL3)
        {
            outputs.push_back({hist_array_full, ""});
        }
        else
        {
            outputs.push_back({hist_array_full, ""});     // default names for full energy
            outputs.push_back({hist_array_axis, "axis"}); // explicit axis output
        }

        for (const auto &out : outputs)
        {
            const string tag_suffix = (out.label.empty() ? "" : "_" + out.label);
            ostringstream oss;
            oss << "ebin_s" << ispec << tag_suffix << ".dat";
            ofstream outfile;
            outfile.open(oss.str().c_str(), ofstream::out);
            for (ibin = 0; ibin + 1 < nbin[histid]; ++ibin)
            {
                outfile << emin[histid] + ibin * de[histid] << "\n";
            }
            outfile << emax[histid];
            outfile.close();
            ostringstream oss_edf_cnt, oss_epf_cnt, oss_edf_norm, oss_epf_norm;
            oss_edf_cnt << "edf_hist_s" << ispec << tag_suffix;
            oss_epf_cnt << "epf_hist_s" << ispec << tag_suffix;
            oss_edf_norm << "edf_hist_norm_s" << ispec << tag_suffix;
            oss_epf_norm << "epf_hist_norm_s" << ispec << tag_suffix;
            if (step_start >= 0)
            {
                oss_edf_cnt << ".dat";
                oss_epf_cnt << ".dat";
                oss_edf_norm << ".dat";
                oss_epf_norm << ".dat";
            }
            else
            {
                oss_edf_cnt << "_step" << step_now << ".dat";
                oss_epf_cnt << "_step" << step_now << ".dat";
                oss_edf_norm << "_step" << step_now << ".dat";
                oss_epf_norm << "_step" << step_now << ".dat";
            }
            ofstream out_edf_cnt, out_epf_cnt, out_edf_norm, out_epf_norm;
            out_edf_cnt.open(oss_edf_cnt.str().c_str(), ofstream::out);
            out_epf_cnt.open(oss_epf_cnt.str().c_str(), ofstream::out);
            out_edf_norm.open(oss_edf_norm.str().c_str(), ofstream::out);
            out_epf_norm.open(oss_epf_norm.str().c_str(), ofstream::out);
            // -----------------------------
            // Each line:
            // ix iy Ntot  [bin0 bin1 ...]
            // For *_norm files:
            // bins store f_E(E) and g(E)
            // -----------------------------
            icell = 0;
            for (int iy = 0; iy < ny; ++iy)
            {
                for (int ix = 0; ix < nx; ++ix)
                {
                    cellid = cell[icell];

                    if (cellid >= 0)
                    {
                        const int Ntot = ncount[cellid][histid];
                        const double dE = de[histid];
                        out_edf_cnt << ix << " " << iy << " " << Ntot;
                        out_epf_cnt << ix << " " << iy << " " << Ntot;
                        out_edf_norm << ix << " " << iy << " " << Ntot;
                        out_epf_norm << ix << " " << iy << " " << Ntot;
                        for (ibin = 0; ibin < nbin[histid]; ++ibin)
                        {
                            const int cnt = out.bins[cellid][histid][ibin];
                            out_edf_cnt << " " << cnt;
                            const double Eedge = emin[histid] + ibin * dE;
                            const double denom_legacy = (Eedge > epsE ? std::sqrt(Eedge) : 1.0);
                            out_epf_cnt << " " << (double)cnt / denom_legacy;
                            double fE = 0.0;
                            if (Ntot > 0 && dE > 0.0)
                            {
                                fE = (double)cnt / ((double)Ntot * dE);
                            }
                            out_edf_norm << " " << fE;
                            const double Ecenter = emin[histid] + (ibin + 0.5) * dE;
                            const double denom = (Ecenter > epsE ? std::sqrt(Ecenter) : 0.0);
                            double gE = (denom > 0.0 ? fE / denom : 0.0);
                            out_epf_norm << " " << gE;
                        }
                        out_edf_cnt << "\n";
                        out_epf_cnt << "\n";
                        out_edf_norm << "\n";
                        out_epf_norm << "\n";
                    }
                    ++icell;
                }
            }
            out_edf_cnt.close();
            out_epf_cnt.close();
            out_edf_norm.close();
            out_epf_norm.close();
        }
    }
    return;
}

void EnergyDist::clear_edf()
{
    auto nhist = nbin.size();
    decltype(nhist) h = 0;
    for (int icell = 0; icell < ncells_sample; ++icell)
    {
        for (h = 0; h < nhist; ++h)
        {
            hist_array_full[icell][h].assign(nbin[h], 0);
            hist_array_axis[icell][h].assign(nbin[h], 0);
            ncount[icell][h] = 0;
        }
    }
    return;
}
