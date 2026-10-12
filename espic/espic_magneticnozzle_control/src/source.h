#ifndef SOURCE_H
#define SOURCE_H

#include <vector>
#include <string>
#include <limits>
#include <cstdlib>
#include <cmath>

#include "part_spec.h"
#include "solve_field.h"

enum SOURCETYPE
{
    MAXW,
    FLUX
};

struct OneSource
{
    int type;                // source type, MAXW or FLUX
    int coll_exist;          // true if collector exist
    int specid;              // species id of this source
    double nflowrate;        // number flow rate, npart per time
    double temp;             // temperatue of source
    double vd[3];            // drifting velocity
    double sxlo[3], sxhi[3]; // source box lower and upper bounds
    double cxlo[3], cxhi[3]; // collector box lower and upper bounds
    std::vector<int> flux_dir[3];
    double wgt; // particle weight of species for this beam

    OneSource(int tp, int sid, double nflw, double T, double v[3],
              double slo[3], double shi[3], double w,
              int flx_d[3][2])
        : type(tp),
          coll_exist(0),
          specid(sid),
          nflowrate(nflw),
          temp(T),
          wgt(w)
    {
        vd[0] = v[0];
        vd[1] = v[1];
        vd[2] = v[2];
        sxlo[0] = slo[0];
        sxlo[1] = slo[1];
        sxlo[2] = slo[2];
        sxhi[0] = shi[0];
        sxhi[1] = shi[1];
        sxhi[2] = shi[2];
        for (int i = 0; i < 3; ++i)
        {
            for (int j = 0; j < 2; ++j)
                if (abs(flx_d[i][j]) == 1)
                    flux_dir[i].push_back(flx_d[i][j]);
        }
    }

    OneSource(int tp, int sid, double nflw, double T, double v[3],
              double slo[3], double shi[3], double clo[3], double chi[3], double w,
              int flx_d[3][2])
        : type(tp),
          coll_exist(1),
          specid(sid),
          nflowrate(nflw),
          temp(T),
          wgt(w)
    {
        vd[0] = v[0];
        vd[1] = v[1];
        vd[2] = v[2];
        sxlo[0] = slo[0];
        sxlo[1] = slo[1];
        sxlo[2] = slo[2];
        sxhi[0] = shi[0];
        sxhi[1] = shi[1];
        sxhi[2] = shi[2];
        cxlo[0] = clo[0];
        cxlo[1] = clo[1];
        cxlo[2] = clo[2];
        cxhi[0] = chi[0];
        cxhi[1] = chi[1];
        cxhi[2] = chi[2];

        for (int i = 0; i < 3; ++i)
        {
            for (int j = 0; j < 2; ++j)
                if (abs(flx_d[i][j]) == 1)
                    flux_dir[i].push_back(flx_d[i][j]);
        }
    }

    const bool collector_exist() const { return coll_exist == 1; }
};

class Sources
{
public:
    typedef std::vector<OneSource>::size_type size_type;

    /* public methods */
    // constructor
    explicit Sources(const std::string & = "");

    // copy constructor
    Sources(const Sources &);

    // desctructor
    ~Sources() {}

    void inject_particles(double, int, long *, const Vec &); // particle injection
    int collect_particle(const double[], const int);         // return src id if hit one collector
    void count_np_collect(int isrc) { ++np_col[isrc]; }
    const bool empty() const { return source.empty(); } // true if no source/collector defined
    const bool collector_exist() const { return _coll_exist; }
    const int np_col_spec(int specid) const { return _np_col_spec[specid]; }
    size_type size() const { return source.size(); }

private:
    std::vector<OneSource> source;
    std::vector<long> np_col;
    std::vector<double> np_res;
    int _np_col_spec[MAXNSPECTYPES];
    bool _coll_exist;
    static constexpr int _dim = 2; // dimension of simulation, hard coded here
    static constexpr double PI2 = 6.2831853071796;
    static constexpr double _small_value = std::numeric_limits<double>::min();

    /* private method */
    double ranf() { return ((double)rand() / (double)RAND_MAX); }

    void init_flux(std::vector<std::string> &, const std::string &, FILE *);
    void init_maxw(std::vector<std::string> &, const std::string &, FILE *);

    void inject_maxw_dist(int, long, long, std::vector<OneSource>::const_iterator &, const Vec &);
    void inject_flux_dist(int, long, long, std::vector<OneSource>::const_iterator &, const Vec &);

    void gen_maxw_vel(const double, const double[], double[], int = -1);
    void gen_flux_vel(const double, const double[], double[], int, const std::vector<int> &);
    void make_one_particle(double[], double[], int, int, double[]);
};

#endif
