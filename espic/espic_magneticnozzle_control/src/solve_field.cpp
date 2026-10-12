#include "solve_field.h"
#include "pic_io.h"
#include "timer.h"
#include "global.h"

void phi2Efield_perd(const int ngpt[3],
                     const double dx[3],
                     const Vec &p,
                     double (*E)[3])
{
    int nx = ngpt[0], ny = ngpt[1], nz = ngpt[2];

    double dx2_inv[3];
    for (int a = 0; a < 3; a++)
        dx2_inv[a] = 1.0 / (2.0 * dx[a]);

#define ID(i, j, k) (((k) * (ny) + (j)) * (nx) + (i))

    for (int k = 0; k < nz; k++)
        for (int j = 0; j < ny; j++)
            for (int i = 0; i < nx; i++)
            {
                long nid = ID(i, j, k);
                int im = i - 1, ip = i + 1;
                int jm = j - 1, jp = j + 1;
                int km = k - 1, kp = k + 1;

                /* periodic boundaries */
                //        if(i == 0) im += (nx-1);
                //        else if(i == nx-1) ip += (1-nx);
                /* fix E = 0 at left and right BC */
                if (i == 0)
                    im = ip;
                else if (i == nx - 1)
                    ip = im;

                /* periodic boundaries */
                if (j == 0)
                    jm += (ny - 1);
                else if (j == ny - 1)
                    jp += (1 - ny);
                if (k == 0)
                    km += (nz - 1);
                else if (k == nz - 1)
                    kp += (1 - nz);

                /* evaluate E field by -grad(phi) by central difference */
                /* E(x) = -grad_x(phi) = -(p(i+1,j,k)-p(i-1,j,k))/(2*dx) */
                E[nid][0] = (p(ID(im, j, k)) - p(ID(ip, j, k))) * dx2_inv[0];
                E[nid][1] = (p(ID(i, jm, k)) - p(ID(i, jp, k))) * dx2_inv[1];
                if (nz > 1)
                    E[nid][2] = (p(ID(i, j, km)) - p(ID(i, j, kp))) * dx2_inv[2];
                else
                    E[nid][2] = 0.0;
            }

#undef ID
}

void phi2Efield(const int ngpt[3],
                const double dx[3],
                const Vec &p,
                double (*E)[3])
{
    int nx = ngpt[0], ny = ngpt[1], nz = ngpt[2];

    double dx2_inv[3];
    for (int a = 0; a < 3; a++)
        dx2_inv[a] = 1.0 / (2.0 * dx[a]);

#define ID(i, j, k) (((k) * (ny) + (j)) * (nx) + (i))

    for (int k = 0; k < nz; k++)
    {
        for (int j = 0; j < ny; j++)
        {
            for (int i = 0; i < nx; i++)
            {
                long nid = ID(i, j, k);
                int im = i - 1, ip = i + 1;
                int jm = j - 1, jp = j + 1;
                int km = k - 1, kp = k + 1;

                /* evaluate E field by -grad(phi) by 2nd order discretization */
                /* E(x) = -grad_x(phi) = -(p(i+1,j,k)-p(i-1,j,k))/(2*dx) */
                if (i == 0)
                {
                    if (fieldBC[0] == FBC_PERIODIC)
                    {
                        im += (nx - 1);
                        E[nid][0] = (p(ID(im, j, k)) - p(ID(ip, j, k))) * dx2_inv[0];
                    }
                    else
                    {
                        int i0 = i, i1 = i + 1, i2 = i + 2;
                        E[nid][0] = (3 * p(ID(i0, j, k)) - 4 * p(ID(i1, j, k)) + p(ID(i2, j, k))) * dx2_inv[0];
                    }
                }
                else if (i == nx - 1)
                {
                    if (fieldBC[1] == FBC_PERIODIC)
                    {
                        ip += (1 - nx);
                        E[nid][0] = (p(ID(im, j, k)) - p(ID(ip, j, k))) * dx2_inv[0];
                    }
                    else
                    {
                        int i0 = i, i1 = i - 1, i2 = i - 2;
                        E[nid][0] = -(3 * p(ID(i0, j, k)) - 4 * p(ID(i1, j, k)) + p(ID(i2, j, k))) * dx2_inv[0];
                    }
                }
                else
                {
                    E[nid][0] = (p(ID(im, j, k)) - p(ID(ip, j, k))) * dx2_inv[0];
                }

                if (j == 0)
                {
                    if (fieldBC[2] == FBC_PERIODIC)
                    {
                        jm += (ny - 1);
                        E[nid][1] = (p(ID(i, jm, k)) - p(ID(i, jp, k))) * dx2_inv[1];
                    }
                    else
                    {
                        int j0 = j, j1 = j + 1, j2 = j + 2;
                        E[nid][1] = (3 * p(ID(i, j0, k)) - 4 * p(ID(i, j1, k)) + p(ID(i, j2, k))) * dx2_inv[1];
                    }
                }
                else if (j == ny - 1)
                {
                    if (fieldBC[3] == FBC_PERIODIC)
                    {
                        jp += (1 - ny);
                        E[nid][1] = (p(ID(i, jm, k)) - p(ID(i, jp, k))) * dx2_inv[1];
                    }
                    else
                    {
                        int j0 = j, j1 = j - 1, j2 = j - 2;
                        E[nid][1] = -(3 * p(ID(i, j0, k)) - 4 * p(ID(i, j1, k)) + p(ID(i, j2, k))) * dx2_inv[1];
                    }
                }
                else
                {
                    E[nid][1] = (p(ID(i, jm, k)) - p(ID(i, jp, k))) * dx2_inv[1];
                }

                if (nz <= 1)
                    E[nid][2] = 0.0;
                else
                {
                    if (k == 0)
                    {
                        if (fieldBC[4] == FBC_PERIODIC)
                        {
                            km += (nz - 1);
                            E[nid][2] = (p(ID(i, j, km)) - p(ID(i, j, kp))) * dx2_inv[2];
                        }
                        else
                        {
                            int k0 = k, k1 = k + 1, k2 = k + 2;
                            E[nid][2] = (3 * p(ID(i, j, k0)) - 4 * p(ID(i, j, k1)) + p(ID(i, j, k2))) * dx2_inv[2];
                        }
                    }
                    if (k == nz - 1)
                    {
                        if (fieldBC[5] == FBC_PERIODIC)
                        {
                            kp += (1 - nz);
                            E[nid][2] = (p(ID(i, j, km)) - p(ID(i, j, kp))) * dx2_inv[2];
                        }
                        else
                        {
                            int k0 = k, k1 = k - 1, k2 = k - 2;
                            E[nid][2] = -(3 * p(ID(i, j, k0)) - 4 * p(ID(i, j, k1)) + p(ID(i, j, k2))) * dx2_inv[2];
                        }
                    }
                    else
                        E[nid][2] = (p(ID(i, j, km)) - p(ID(i, j, kp))) * dx2_inv[2];
                }
            }
        }
    }

#undef ID
}

int bicgstab_solver(SpMat &A,
                    const Vec &b,
                    Vec &p,
                    const int ngpt[3],
                    const MapF &fixedBV)
{
    static bool first = true;
    static Eigen::BiCGSTAB<SpMat, Eigen::IncompleteLUT<double>> bicgSolver;

    if (first)
    {
        bicgSolver.compute(A);
        if (bicgSolver.info() != Eigen::Success)
        {
            cerr << "bicgSolver compute failed!\n";
            exit(-1);
        }
        bicgSolver.setTolerance(RTOL);
        bicgSolver.setMaxIterations(MAXIT);
        first = false;
    }

    //  p = bicgSolver.solve(b);
    p = bicgSolver.solveWithGuess(b, p);

    cout << ", iter/max_iter = " << bicgSolver.iterations() << "/" << bicgSolver.maxIterations()
         << ", error/tol = " << bicgSolver.error() << "/" << bicgSolver.tolerance();

    return bicgSolver.info() == Eigen::Success;
}

int lu_solver(SpMat &A,
              const Vec &b,
              Vec &p,
              const int ngpt[3],
              const MapF &fixedBV)
{
    static bool first = true;
    static Eigen::SparseLU<SpMat> luSolver;

    if (first)
    {
        luSolver.compute(A);
        if (luSolver.info() != Eigen::Success)
        {
            cerr << "LU decomposition failed!\n";
            exit(-1);
        }
        first = false;
    }

    p = luSolver.solve(b);

    return luSolver.info() == Eigen::Success;
}

void solve_field_linear(const int ngpt[3],
                        const int nspec,
                        SpMat &Mat_A,
                        const MapF &fixedBV,
                        double **chgDenst,
                        Vec &p,
                        Vec &b,
                        double (*E)[3], double (*Eext)[3],
                        PtrPoissonSolver poisson_solver)
{
    long nxyz = ngpt[0] * ngpt[1] * ngpt[2];

    /* update rhs for field solver (b[nid] = -rho_tot[nid]) */
    update_rhs_linearFDS(ngpt, nspec, fixedBV, chgDenst, b);

    /* solver for phi with selected poisson solver */
    poisson_solver(Mat_A, b, p, ngpt, fixedBV);

    /* compute electric field from potential field */
    //  phi2Efield_perd(ngpt, dx, p, E);  // periodic in ymin and ymax
    phi2Efield(ngpt, dx, p, E); // 2nd order discretization; one-sided discret at BC

    if (EextOn)
    {
        for (long nid = 0; nid < nxyz; nid++)
            for (int a = 0; a < 3; a++)
                E[nid][a] += Eext[nid][a];
    }
    //  for(long nid = 0; nid < nxyz; nid++) p[nid] = 0.0;
}

int lu_nonlinear(SpMat &A,
                 const Vec &b,
                 Vec &p,
                 const int ngpt[3],
                 const MapF &fixedBV)
{
    Eigen::SparseLU<SpMat> luSolver;
    SpMat J(A);
    Vec f(b);
    int iter_tot = 0, iter_nt;
    bool is_converged = false;
    double residule;

    /* Newton's method to solve non-linear Poisson's equation */
    /* Within each Newton's iteration, LU solver is used */
    for (iter_nt = 0; (iter_nt < MAXITNT) && (!is_converged); iter_nt++)
    {
        update_rhs_nonlinearFDS(ngpt, fixedBV, b, p, A, f);
        update_matrix_nonlinearFDS(ngpt, fixedBV, p, A, J);

        luSolver.compute(J);
        if (luSolver.info() != Eigen::Success)
        {
            cerr << "LU decomposition failed!\n";
            exit(-1);
        }

        f = luSolver.solve(f);
        p = p - f;

        residule = f.norm();
        is_converged = residule < RTOLNT;
        // cout << "LU: " << "iter_nt = " << iter_nt << ", residual = " << residule << "\n";
    }

    cout << ", iter_nt/maxiter_nt = " << iter_nt << "/" << MAXITNT << ", iter_tot = " << iter_tot
         << ", error/tol = " << residule << "/" << RTOLNT;

    return (int)is_converged;
}

int bicgstab_nonlinear(SpMat &A,
                       const Vec &b,
                       Vec &p,
                       const int ngpt[3],
                       const MapF &fixedBV)
{
    Eigen::BiCGSTAB<SpMat, Eigen::IncompleteLUT<double>> bicgSolver;
    SpMat J(A);
    Vec f(b);
    int iter_tot = 0, iter_nt;
    bool is_converged = false;
    double residule;

    /* Newton's method to solve non-linear Poisson's equation */
    /* Within each Newton's iteration, bi-conjugate iterative solver is used */
    for (iter_nt = 0; (iter_nt < MAXITNT) && (!is_converged); iter_nt++)
    {
        update_rhs_nonlinearFDS(ngpt, fixedBV, b, p, A, f);
        update_matrix_nonlinearFDS(ngpt, fixedBV, p, A, J);

        bicgSolver.compute(J);
        if (bicgSolver.info() != Eigen::Success)
        {
            cerr << "bicgSolver compute failed!\n";
            exit(-1);
        }
        bicgSolver.setTolerance(RTOL);
        bicgSolver.setMaxIterations(MAXIT);

        //    f = bicgSolver.solve(f);
        f = bicgSolver.solveWithGuess(f, f);
        p = p - f;

        residule = f.norm();
        is_converged = residule < RTOLNT;

        iter_tot += bicgSolver.iterations();
        // cout << "BICGSTAB: " << "iter_nt = " << iter_nt << ", residual = " << residule << "\n";
    }

    cout << ", iter_nt/maxiter_nt = " << iter_nt << "/" << MAXITNT << ", iter_tot = " << iter_tot
         << ", error/tol = " << residule << "/" << RTOLNT;

    return (int)is_converged;
}

void solve_field_nonlinear(const int ngpt[3],
                           const int nspec,
                           SpMat &Mat_A,
                           const MapF &fixedBV,
                           double **chgDenst,
                           Vec &p,
                           Vec &b,
                           double (*E)[3], double (*Eext)[3],
                           PtrPoissonSolver poisson_solver)
{
    long nxyz = ngpt[0] * ngpt[1] * ngpt[2];
    //  Eigen::SparseLU< SpMat > luSolver;

    /* store charge density & FBC value into b */
    update_rhs_linearFDS(ngpt, nspec, fixedBV, chgDenst, b);

    /* solver for phi with selected poisson solver */
    poisson_solver(Mat_A, b, p, ngpt, fixedBV);

    /* compute electric field from potential field */
    phi2Efield(ngpt, dx, p, E); // 2nd order discretization; one-sided discret at BC

    if (EextOn)
    {
        for (long nid = 0; nid < nxyz; nid++)
            for (int a = 0; a < 3; a++)
                E[nid][a] += Eext[nid][a];
    }
    //  for(long nid = 0; nid < nxyz; nid++) p[nid] = 0.0;
}
