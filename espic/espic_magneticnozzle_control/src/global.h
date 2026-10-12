#ifndef GLOBAL_H
#define GLOBAL_H

#include "part_spec.h"

const double smallValue = 1e-10;
const int MAXNOBJECTS = 10;

/* flags for particle BC */
const int PBC_ABSORPT = 1;
const int PBC_REFLECT = 2;
const int PBC_PERIODIC = 4;

/* flags for field BC */
const int FBC_NEUMANN = 1;
const int FBC_DIRICHLET = 2;
const int FBC_PERIODIC = 4;

/* inflow & outflow */
const int MAXNFLOWS = 10;

extern double E0[3];
extern double B0[3];
extern bool EextOn;
extern bool BextOn;

extern double dx[3];
extern double dxinv[3];
extern double domain[6];
extern int fieldBC[6];
extern int partBC[6];
extern int partBC_new[6];
extern double fbcValue[6];

extern double objectCoord[MAXNOBJECTS][6];
extern double objectBCvalue[MAXNOBJECTS];
extern int nobjects;

extern int MAXNPARTS;
extern int MAXNPARTSINCR;

extern double (*part)[MAXNPARTPROPS];
extern double spec[MAXNSPECTYPES][MAXNSPECPROPS];
extern double partWeight[MAXNSPECTYPES];

extern double *volinv; // for axisymmetric
#endif
