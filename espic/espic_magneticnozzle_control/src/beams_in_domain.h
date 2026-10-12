#ifndef BEAMS_IN_DOMAIN_H
#define BEAMS_IN_DOMAIN_H

#include <string>
#include <cstdlib>
#include "part_spec.h"
#include "solve_field.h"

struct OneBeam
{
    double xc[3];      // position of beam center
    int dim;           // dimension: 2 - 2D, 3 - 3D
    double r;          // beam radius
    double area;       // beam injection plane area: 2*r - 2D, pi*r*r - 3D
    double nvec[3];    // normal vector of injection plane
    double tvec[2][3]; // 2 tangential vectors of inj plane
    int specid;        // specid id
    double ndens;      // number density
    double temp;       // temperature
    double vd[3];      // drifting velocity
    double wgt;        // particle weight of species for this beam
    double sn;         // speed ratio of v_normal to beta(=sqrt(2kT/m))
    double p_norm;     // normalization factor for inject probability
};

class Beams
{
public:
    typedef int size_type;

    /* public methods */
    // constructor (default & read from file)
    Beams(const std::string &);

    // copy constructor
    Beams(const Beams &);

    // desctructor
    ~Beams();

    // beam injection
    void inject_particles(double, int, long *, const Vec &);
    //    int inject_particles(double (*)[MAXNPARTPROPS], long *);

    size_type size() const { return nbeams; }

private:
    OneBeam *beam;
    size_type nbeams;
    double fnparts[MAXNSPECTYPES];

    /* private methods */
    // beam injection initialization
    void inject_init();

    double ranf() { return ((double)rand() / (double)RAND_MAX); }
};

#endif
