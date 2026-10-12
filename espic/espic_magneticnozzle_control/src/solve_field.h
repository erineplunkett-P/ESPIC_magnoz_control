#ifndef POISSON_SOLVER_H
#define POISSON_SOLVER_H

#include <iostream>
#include <iomanip>
#include <fstream>
#include <vector>
#include <string>
#include <sstream>
#include <map>
#include <utility>

#include <cstdlib>
#include <cstdio>
#include <cstring>
#include <cmath>

#include <memory.h>

#include "Eigen/Sparse"
#include "objects_in_domain.h"

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

using std ::fixed;
using std ::setprecision;

typedef Eigen::SparseMatrix<double> SpMat;
typedef Eigen::Triplet<double> Tp;
typedef Eigen::VectorXd Vec;

typedef std::map<long, double> MapF;
using std::make_pair;

// PtrPoissonSolver is the type name of ptr that points to a poisson solver function
typedef int (*PtrPoissonSolver)(SpMat &, const Vec &, Vec &, const int[], const MapF &);

extern int lu_solver(SpMat &, const Vec &, Vec &, const int[], const MapF &);
extern int bicgstab_solver(SpMat &, const Vec &, Vec &, const int[], const MapF &);
extern int lu_nonlinear(SpMat &, const Vec &, Vec &, const int[], const MapF &);
extern int bicgstab_nonlinear(SpMat &, const Vec &, Vec &, const int[], const MapF &);

extern void solve_field_linear(const int ngpt[3],
                               const int nspec,
                               SpMat &Mat_A,
                               const MapF &fixedBV,
                               double **chgDenst,
                               Vec &p,
                               Vec &b,
                               double (*E)[3], double (*Eext)[3],
                               PtrPoissonSolver poisson_solver);

extern void solve_field_nonlinear(const int ngpt[3],
                                  const int nspec,
                                  SpMat &A,
                                  const MapF &fixedBV,
                                  double **chgDenst,
                                  Vec &p,
                                  Vec &b,
                                  double (*E)[3], double (*Eext)[3],
                                  PtrPoissonSolver poisson_solver);

extern void setup_matrix_linearFDS(const int ngpt[3],
                                   const double *dx,
                                   SpMat &Mat_A,
                                   const ObjectList &objlist,
                                   MapF &fixedBV);

extern void setup_rhs_linearFDS(const int ngpt[3], Vec &b);

extern void update_rhs_linearFDS(const int ngpt[3],
                                 const int nspec,
                                 const MapF &fixedBV,
                                 double **chgDenst,
                                 Vec &b);

extern void update_matrix_nonlinearFDS(const int ngpt[3],
                                       const MapF &fixedBV,
                                       Vec &p,
                                       SpMat &A,
                                       SpMat &J);

extern void update_rhs_nonlinearFDS(const int ngpt[3],
                                    const MapF &fixedBV,
                                    const Vec &b,
                                    Vec &p,
                                    SpMat &A,
                                    Vec &f);

const double PI = 2.0 * asin(1.0);

/* define solver type */
const int CG_SOLVER = 1;
const int LU_SOLVER = 4;
const int OUTPUT_PRECS = 12;

extern double RTOL;
extern int MAXIT;
extern double RTOLNT;
extern int MAXITNT;

#endif /* POISSON_SOLVER_H */
