#ifndef PART_SPEC_H
#define PART_SPEC_H

#define MAXNPARTPROPS 14
#define MAXNSPECTYPES 5
#define MAXNSPECPROPS 6

enum partProp
{
    xp,
    yp,
    zp,
    up,
    vp,
    wp,
    ptype,
    pid,
    xp0,
    yp0,
    up0,
    vp0,
    Epot,
    partPropN
};

enum specProp
{
    specmass,
    specchg,
    specQM,    // chg/mass
    specFixed, // flag to indicate a spec is fixed (1) or not (0)
    specPropN
};

#endif
