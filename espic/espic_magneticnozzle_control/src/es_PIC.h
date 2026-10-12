#ifndef ES_PIC_H
#define ES_PIC_H

#include <iostream>
#include <vector>
#include <string>
#include <sstream>
#include <fstream>

#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <cstddef>
#include <sys/stat.h>

#include <memory.h>

#include "solve_field.h"
#include "parallel.h"
#include "global.h"
#include "part_spec.h"
#include "pic_part.h"
#include "pic_io.h"
#include "pic_math.h"
#include "objects_in_domain.h"
#include "beams_in_domain.h"
#include "source.h"
#include "sample_field.h"
#include "vdf.h"
#include "timer.h"
#include "comp_ellint.h"

using std ::cerr;
using std ::cin;
using std ::cout;
using std ::endl;

using std ::vector;

using std ::istringstream;
using std ::ostringstream;
using std ::string;

using std ::ifstream;
using std ::ofstream;

const int MB = 1024 * 1024;
const int KB = 1024;

void injectFlowInOut(int nspec,
                     long *nparts,
                     int ngpt[],
                     double *partWeight,
                     double dt, int istep,
                     const Vec &p);

void injectBeam(int nspec,
                long *nparts,
                double *partWeight,
                double dt, int istep);

void evalHistVelField(int nspec,
                      const long *nparts,
                      const int ngpt[3],
                      const double (*part)[MAXNPARTPROPS],
                      vector<double> &histVelField);

void evalNewHistVelField(int nspec,
                         const long *nparts,
                         const int ngpt[3],
                         const double (*part)[MAXNPARTPROPS],
                         vector<double> &new_histVelField);

void evalNewNewHistVelField(int nspec,
                            const long *nparts,
                            const int ngpt[3],
                            const double (*part)[MAXNPARTPROPS],
                            vector<double> &new_new_histVelField);

void evalAzimuthal(
    int nspec,
    const long *nparts,
    const int ngpt[3],
    const double (*part)[MAXNPARTPROPS],
    const double *particleWeight,
    std::vector<double> &number_density,
    std::vector<double> &vel_azimuthal,
    std::vector<double> &momentum_azimuthal,                // M_{rr}, M_{rz}, M_{theta,theta} (convective)
    std::vector<double> &stress_azimuthal,                  // S_{rr}, S_{rz}, S_{theta,theta} (thermal)
    std::vector<double> &ExB_component,                     // e n_e * (E_r/B_z) or -e n_e * (E_z/B_r)
    std::vector<double> &diamagnetic_component,             // (1/B) * \partial p_e
    std::vector<double> &inertial_component,                // (1/B) * (\nabla\cdot M_e terms)
    std::vector<double> &stressdiv_component,               // -(1/B) * (\nabla\cdot \pi_e terms)
    std::vector<double> &Je_theta_total_fluid_fullmomentum, // sum of the four
    std::vector<double> &Je_theta_PIC,                      // from particle data directly
    const std::vector<double> &histVelField,                // 6*nspec*nxyz; using u_z,u_r,u_theta
    double (*E)[3],
    double (*B)[3],
    int eqSelect = 0, // 1: Eq.20; 2: Eq.21
    int volSelect = 1 // 1: node volume, volinv; 2: cell volume, calculated locally
);

#endif
