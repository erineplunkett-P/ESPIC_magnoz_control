#include "Eigen/Sparse"
typedef Eigen::VectorXd Vec;

#include <iostream>
#include <sstream>
#include <cstdio>
#include <cstring>

#include "source.h"
#include "parse.h"
#include "global.h"
#include "pic_io.h"
#include "pic_part.h"

using std::cout;
using std::string;
using std::vector;

Sources::Sources(const string &filename)
    : source(),
      np_col(),
      np_res(),
      _coll_exist(false)
{
    FILE *fp = fopen(filename.c_str(), "r");
    if (NULL == fp)
    {
        cout << "No source initialized." << std::endl;
        for (int specid = 0; specid < MAXNSPECTYPES; ++specid)
            _np_col_spec[specid] = 0;

        return;
    }

    cout << "read source from[" << filename << "]." << std::endl;

    vector<string> word;
    while (parse_line(word, fp))
    {
        if (word.empty())
            continue; // this is a comment line

        if ("flux_dist" == word[0])
            init_flux(word, filename, fp);
        else if ("maxw_dist" == word[0])
            init_maxw(word, filename, fp);
        else
        {
            std::ostringstream oss;
            oss << "Illegal input setting in" << filename;
            esPIC_exit(oss.str());
        }
    }
    fclose(fp);

    np_col.resize(source.size(), 0);
    np_res.resize(source.size(), 0.);

    return;
}

Sources::Sources(const Sources &orig)
    : source(orig.source),
      np_col(orig.np_col),
      np_res(orig.np_res),
      _coll_exist(orig._coll_exist)
{
}

void Sources::init_flux(vector<string> &word,
                        const string &filename,
                        FILE *fp)
{
    int sid, flag_collector = 0;
    int flx_d[3][2] = {{0, 0}, {0, 0}, {0, 0}};
    double T, nflw, w, v[3], slo[3], shi[3], clo[3], chi[3];

    // process flux direction
    if (word.size() > 7)
    {
        std::ostringstream oss;
        oss << "Illegal input setting in" << filename;
        esPIC_exit(oss.str());
    }
    for (auto it = word.cbegin() + 1; it != word.cend(); ++it)
    {
        if ("-x" == *it)
            flx_d[0][0] = -1;
        else if ("+x" == *it)
            flx_d[0][1] = 1;
        else if ("-y" == *it)
            flx_d[1][0] = -1;
        else if ("+y" == *it)
            flx_d[1][1] = 1;
        else if ("-z" == *it)
            flx_d[2][0] = -1;
        else if ("+z" == *it)
            flx_d[2][1] = 1;
    }

    // reading: specid, nflow_rate, T, vx, vy, vz
    parse_line(word, fp);
    if (word.size() != 6)
    {
        std::ostringstream oss;
        oss << "Illegal input setting in" << filename;
        esPIC_exit(oss.str());
    }

    sid = (int)atoi(word[0].c_str()) - 1; // [0, nspec)
    nflw = (double)atof(word[1].c_str());
    T = (double)atof(word[2].c_str());
    for (int a = 0; a < 3; a++)
        v[a] = (double)atof(word[3 + a].c_str());
    w = partWeight[sid];

    // reading: slo[3], shi[3], (clo[3], chi[3])
    parse_line(word, fp);
    if (word.size() == 6)
        flag_collector = 0;
    else if (word.size() == 12)
        flag_collector = 1;
    else
    {
        std::ostringstream oss;
        oss << "Illegal input setting in" << filename;
        esPIC_exit(oss.str());
    }

    for (int a = 0; a < 3; ++a)
    {
        slo[a] = (double)atof(word[a].c_str());
        shi[a] = (double)atof(word[a + 3].c_str());
    }

    if (flag_collector)
    {
        for (int a = 0; a < 3; ++a)
        {
            clo[a] = (double)atof(word[a + 6].c_str());
            chi[a] = (double)atof(word[a + 9].c_str());
        }
        source.push_back(
            OneSource(FLUX, sid, nflw, T, v, slo, shi, clo, chi, w, flx_d));
    }
    else
    {
        source.push_back(
            OneSource(FLUX, sid, nflw, T, v, slo, shi, w, flx_d));
    }

    if (!_coll_exist && flag_collector)
        _coll_exist = true;
}

void Sources::init_maxw(vector<string> &word,
                        const string &filename,
                        FILE *fp)
{
    int sid, flag_collector = 0;
    int flx_d[3][2] = {{0, 0}, {0, 0}, {0, 0}};
    double T, nflw, w, v[3], slo[3], shi[3], clo[3], chi[3];

    // process flux direction
    if (word.size() > 1)
    {
        std::ostringstream oss;
        oss << "Illegal input setting in" << filename;
        esPIC_exit(oss.str());
    }

    // reading: specid, nflow_rate, T, vx, vy, vz
    parse_line(word, fp);
    if (word.size() != 6)
    {
        std::ostringstream oss;
        oss << "Illegal input setting in" << filename;
        esPIC_exit(oss.str());
    }

    sid = (int)atoi(word[0].c_str()) - 1; // [0, nspec)
    nflw = (double)atof(word[1].c_str());
    T = (double)atof(word[2].c_str());
    for (int a = 0; a < 3; a++)
        v[a] = (double)atof(word[3 + a].c_str());
    w = partWeight[sid];

    // reading: slo[3], shi[3], (clo[3], chi[3])
    parse_line(word, fp);
    if (word.size() == 6)
        flag_collector = 0;
    else if (word.size() == 12)
        flag_collector = 1;
    else
    {
        std::ostringstream oss;
        oss << "Illegal input setting in" << filename;
        esPIC_exit(oss.str());
    }

    for (int a = 0; a < 3; ++a)
    {
        slo[a] = (double)atof(word[a].c_str());
        shi[a] = (double)atof(word[a + 3].c_str());
    }

    if (flag_collector)
    {
        for (int a = 0; a < 3; ++a)
        {
            clo[a] = (double)atof(word[a + 6].c_str());
            chi[a] = (double)atof(word[a + 9].c_str());
        }
        source.push_back(
            OneSource(MAXW, sid, nflw, T, v, slo, shi, clo, chi, w, flx_d));
    }
    else
    {
        source.push_back(
            OneSource(MAXW, sid, nflw, T, v, slo, shi, w, flx_d));
    }

    if (!_coll_exist && flag_collector)
        _coll_exist = true;
}

void Sources::inject_particles(double dt, int istep, long *np, const Vec &pot)
{
    if (empty())
        return;

    int isrc = 0;
    for (auto it = source.cbegin(); it != source.cend(); ++it, ++isrc)
    {
        long np_tot = np[MAXNSPECTYPES];
        double fnp_gen = (it->nflowrate * dt) / it->wgt;
        long np_gen = long(fnp_gen + np_res[isrc]);

        np_res[isrc] += (fnp_gen - np_gen); // adjust # of residual particle to be injected
        np_gen += np_col[isrc];

        if (0 == np_gen)
            continue;

        np[it->specid] += np_gen;    // # of particles of this species
        np[MAXNSPECTYPES] += np_gen; // # of particles of all species
        while (np[MAXNSPECTYPES] > MAXNPARTS)
            adjustPartMem();

        if (FLUX == it->type)
            inject_flux_dist(istep, np_tot, np_gen, it, pot);
        else if (MAXW == it->type)
            inject_maxw_dist(istep, np_tot, np_gen, it, pot);
    }

    // count # of particles hitted coll for each species
    for (int specid = 0; specid < MAXNSPECTYPES; ++specid)
        _np_col_spec[specid] = 0;

    isrc = 0;
    for (auto it = source.cbegin(); it != source.cend(); ++it, ++isrc)
    {
        if (!it->collector_exist())
            continue;
        _np_col_spec[it->specid] += np_col[isrc];
    }

    // reset # of particle needs to be gen in this step
    np_col.assign(np_col.size(), 0);
    //  for (auto it = np_col.begin(); it != np_col.end(); ++it) *it = 0;
}

void Sources::inject_maxw_dist(int istep, long np_tot, long np_gen,
                               vector<OneSource>::const_iterator &it,
                               const Vec &pot)
{
    double p[MAXNPARTPROPS];
    double x0[3], dx[3], vth;
    double x[3] = {0., 0., 0.}, v[3] = {0., 0., 0.};

    vth = sqrt(it->temp / spec[it->specid][specmass]);
    for (int a = 0; a < 3; a++)
    {
        x0[a] = it->sxlo[a];
        dx[a] = it->sxhi[a] - x0[a];
    }

    if (_dim == 2)
    { // 2d simulation
        x0[2] = 0.0;
        dx[2] = 0.0;

        int nx = (int)((domain[1] - domain[0]) * dxinv[0] + 0.00001) + 1;
        for (long ip = 0; ip < np_gen; ++ip)
        {
            for (int a = 0; a < _dim; a++)
                x[a] = x0[a] + dx[a] * ranf();

            gen_maxw_vel(vth, it->vd, v);

            make_one_particle(x, v, it->specid, istep, p);

            // initial info of a particle
            p[xp0] = p[xp];
            p[yp0] = p[yp];
            p[up0] = p[up];
            p[vp0] = p[vp];

            // get potenial of nearest nodes
            double xd[2] = {p[xp0] - domain[0], p[yp0] - domain[2]};
            int ix = (int)(xd[0] * dxinv[0] + 0.5), iy = (int)(xd[1] * dxinv[1] + 0.5);
            p[Epot] = spec[it->specid][specchg] * pot(iy * nx + ix);

            for (int a = 0; a < partPropN; a++)
                part[np_tot][a] = p[a];
            ++np_tot;
        }
    } // if (_dim == 2)
    else
    {
        esPIC_exit("simulations other than 2d are not currently available");
    }

    return;
}

void Sources::inject_flux_dist(int istep, long np_tot, long np_gen,
                               vector<OneSource>::const_iterator &it,
                               const Vec &pot)
{
    double p[MAXNPARTPROPS];
    double x0[3], dx[3], vth;
    double x[3] = {0., 0., 0.}, v[3] = {0., 0., 0.};

    vth = sqrt(it->temp / spec[it->specid][specmass]);
    for (int a = 0; a < 3; a++)
    {
        x0[a] = it->sxlo[a];
        dx[a] = it->sxhi[a] - x0[a];
    }

    if (_dim == 2)
    { // 2d simulation
        x0[2] = 0.0;
        dx[2] = 0.0;

        int nx = (int)((domain[1] - domain[0]) * dxinv[0] + 0.00001) + 1;
        for (long ip = 0; ip < np_gen; ++ip)
        {
            for (int a = 0; a < _dim; a++)
                x[a] = x0[a] + dx[a] * ranf();

            for (int d = 0; d < _dim; ++d)
            { // 2D
                if (it->flux_dir[d].empty())
                    gen_maxw_vel(vth, it->vd, v, d);
                else
                    gen_flux_vel(vth, it->vd, v, d, it->flux_dir[d]);
            }

            make_one_particle(x, v, it->specid, istep, p);

            // initial info of a particle
            p[xp0] = p[xp];
            p[yp0] = p[yp];
            p[up0] = p[up];
            p[vp0] = p[vp];

            // get potenial of nearest nodes
            double xd[2] = {p[xp0] - domain[0], p[yp0] - domain[2]};
            int ix = (int)(xd[0] * dxinv[0] + 0.5), iy = (int)(xd[1] * dxinv[1] + 0.5);
            p[Epot] = spec[it->specid][specchg] * pot(iy * nx + ix);

            for (int a = 0; a < partPropN; a++)
                part[np_tot][a] = p[a];
            ++np_tot;
        }
    } // if (_dim == 2)
    else
    {
        esPIC_exit("simulations other than 2d are not currently available");
    }

    return;
}

void Sources::gen_maxw_vel(const double vth,
                           const double vd[3],
                           double v[3],
                           int d)
{
    double rf1, rf2;

    if (-1 == d)
    { // default, v in all direction follows maxwellian
        rf1 = vth * sqrt(-2.0 * log(ranf() + _small_value));
        rf2 = PI2 * ranf();
        v[0] = rf1 * cos(rf2) + vd[0];
        v[1] = rf1 * sin(rf2) + vd[1];

        if (_dim == 3)
        {
            rf1 = vth * sqrt(-2.0 * log(ranf() + _small_value));
            rf2 = PI2 * ranf();
            v[2] = rf1 * cos(rf2) + vd[2];
        }
    }
    else if (d >= 0 && d < _dim)
    {
        rf1 = vth * sqrt(-2.0 * log(ranf() + _small_value));
        rf2 = PI2 * ranf();
        v[d] = rf1 * cos(rf2) + vd[d];
    }
    else
    {
        esPIC_exit("Exceeds the dimension when generating Maxwellian velocity");
    }

    return;
}

void Sources::gen_flux_vel(const double vth,
                           const double vd[3],
                           double v[3],
                           int d,
                           const vector<int> &flx_d)
{
    int nd = flx_d.size();

    if (d >= 0 && d < _dim)
    {
        v[d] = vth * sqrt(-2 * log(ranf() + _small_value));
        v[d] *= flx_d[(int)(nd * ranf() - _small_value)]; // v[d] * (-1 or 1)
        v[d] += vd[d];
    }
    else
    {
        esPIC_exit("Exceeds the dimension when generating Maxwellian velocity");
    }

    return;
}

void Sources::make_one_particle(double x[3],
                                double v[3],
                                int specid,
                                int istep,
                                double p[partPropN])
{
    for (int a = 0; a < 3; a++)
    {
        p[xp + a] = x[a];
        p[up + a] = v[a];
    }
    p[ptype] = specid;
    p[pid] = istep;

    return;
}

int Sources::collect_particle(const double r[3], const int specid)
{
    if (!collector_exist())
        return (-1);

    bool flag = false;
    int isrc = 0;
    for (auto it = source.cbegin(); it != source.cend(); ++it, ++isrc)
    {
        if (!it->collector_exist())
            continue;
        if (specid != it->specid)
            continue;

        flag = true;
        for (int a = 0; a < 3 && flag; ++a)
        {
            flag = flag && (r[a] >= it->cxlo[a] && r[a] <= it->cxhi[a]);
        }

        // hit this collector
        if (flag)
            break;
    }

    if (flag == false)
        isrc = -1;

    return isrc;
}
