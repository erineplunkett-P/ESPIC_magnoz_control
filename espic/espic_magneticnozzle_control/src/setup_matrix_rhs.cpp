#include "es_PIC.h"
#include "global.h"

inline bool is_fixedBC(int a)
{
    return a == 1;
}

double Radius(double j, double dr)
{
    return j * dr;
}

void innerNodes(const int i, const int j, const int k,
                const int ngpt[3], const double *dx,
                long col[5], double f[5])
{
    int nx = ngpt[0], ny = ngpt[1];
    double b = -1.0 / (dx[0] * dx[0]), c = -1.0 / (dx[1] * dx[1]);
    double a, d;

    a = -2.0 * (b + c);
    d = (j == 0 ? 0.0 : 1 / (Radius(j, dx[1]) * 2. * dx[1])); // symmetric axis if (j == 0)

#define ID(i, j, k) (((k) * (ny) + (j)) * (nx) + (i))

    col[0] = ID(i, j - 1, k);
    f[0] = c + d;
    col[1] = ID(i - 1, j, k);
    f[1] = b;
    col[2] = ID(i, j, k);
    f[2] = a;
    col[3] = ID(i + 1, j, k);
    f[3] = b;
    col[4] = ID(i, j + 1, k);
    f[4] = c - d;

#undef ID
}

int bcNeumannPeriodic(const int i, const int j, const int k,
                      const int ngpt[3], const double *dx,
                      int meshFunc,
                      long col[5], double f[5])
{
    int nx = ngpt[0], ny = ngpt[1];
    double b = -1.0 / (dx[0] * dx[0]), c = -1.0 / (dx[1] * dx[1]);

#define ID(i, j, k) (((k) * (ny) + (j)) * (nx) + (i))

    /* lower boundary */
    if (j == 0)
    {
        if (fieldBC[2] == FBC_NEUMANN)
        {
            f[4] += c;   // correct coeff (i, j+1, k)
            col[0] = -1; // void node (i, j-1, k)
            meshFunc += 2;
        }
        else
        {
            esPIC_exit("Lower boundary can only be zero-Neumann");
        }
    }

    /* left boundary */
    if (i == 0)
    {
        if (fieldBC[0] == FBC_NEUMANN)
        {
            f[3] += b;   // correct coeff (i+1, j, k)
            col[1] = -1; // void node (i-1, j, k)
            meshFunc += 2;
        }
        else if (fieldBC[0] == FBC_PERIODIC)
        {
            col[1] = ID(i + nx - 2, j, k);
            meshFunc += 8;
        }
    }

    /* right boundary */
    if (i == nx - 1)
    {
        if (fieldBC[1] == FBC_NEUMANN)
        {
            f[1] += b;   // correct coeff (i-1, j, k)
            col[3] = -1; // void node (i+1, j, k)
            meshFunc += 2;
        }
        else if (fieldBC[1] == FBC_PERIODIC)
        {
            col[3] = ID(i - nx + 2, j, k);
            meshFunc += 8;
        }
    }

    /* upper boundary */
    if (j == ny - 1)
    {
        if (fieldBC[3] == FBC_NEUMANN)
        {
            f[0] += c;                                   // correct coeff (i, j-1, k)
            f[0] -= 1 / (2. * dx[1] * Radius(j, dx[1])); // axisymmetric
            col[4] = -1;                                 // void node (i, j+1, k)
            meshFunc += 2;
        }
        else
        {
            esPIC_exit("Upper boundary can not be periodic");
        }
    }

#undef ID

    return meshFunc;
}

/* Two-dimensional case with periodic BC at all boundaries for linear field solver. */
void setup_matrix_linearFDS(const int ngpt[3],
                            const double *dx,
                            SpMat &A,
                            const ObjectList &objlist,
                            MapF &fixedBV)
{
    int nx = ngpt[0], ny = ngpt[1], nz = ngpt[2];
    int nxyz = nx * ny * nz;
    ObjectList::size_type nobjs = objlist.size();

    int *meshFunc = new int[nxyz];
    memset(meshFunc, 0, nxyz * sizeof(int));

#define ID(i, j, k) (((k) * (ny) + (j)) * (nx) + (i))

    int i, j, k = 0;

    /* objects in the simulation domain, fixed BC */
    for (j = 0; j < ny; j++)
    {
        for (i = 0; i < nx; i++)
        {
            double x[3] = {i * dx[0], j * dx[1], k * dx[2]};

            for (ObjectList::size_type iobj = 0; iobj < nobjs; iobj++)
                if (is_inside_object(x, objlist[iobj]))
                {
                    meshFunc[ID(i, j, k)] = 1;
                    fixedBV.insert(make_pair(ID(i, j, k), objlist[iobj].potential(x)));
                }
        }
    }

    /* left boundary */
    i = 0;
    if (fieldBC[0] == FBC_DIRICHLET)
    {
        for (j = 0; j < ny; j++)
        {
            meshFunc[ID(i, j, k)] = 1;
            fixedBV.insert(make_pair(ID(i, j, k), fbcValue[0]));
        }
    }
    else if (fieldBC[0] == FBC_NEUMANN)
    {
        for (j = 0; j < ny; j++)
        {
            long nid = ID(i, j, k);

            if (is_fixedBC(meshFunc[nid]))
                continue;

            meshFunc[nid] -= 2;
        }
    }
    else if (fieldBC[0] == FBC_PERIODIC)
    {
        for (j = 0; j < ny; j++)
        {
            long nid = ID(i, j, k);

            if (is_fixedBC(meshFunc[nid]))
                continue;

            meshFunc[nid] -= 8;
        }
    }

    /* right boundary */
    i = nx - 1;
    if (fieldBC[1] == FBC_DIRICHLET)
    {
        for (j = 0; j < ny; j++)
        {
            meshFunc[ID(i, j, k)] = 1;
            fixedBV.insert(make_pair(ID(i, j, k), fbcValue[1]));
        }
    }
    else if (fieldBC[1] == FBC_NEUMANN)
    {
        for (j = 0; j < ny; j++)
        {
            long nid = ID(i, j, k);

            if (is_fixedBC(meshFunc[nid]))
                continue;

            meshFunc[nid] -= 2;
        }
    }
    else if (fieldBC[1] == FBC_PERIODIC)
    {
        for (j = 0; j < ny; j++)
        {
            long nid = ID(i, j, k);

            if (is_fixedBC(meshFunc[nid]))
                continue;

            meshFunc[nid] -= 8;
        }
    }

    /* lower boundary */
    j = 0;
    if (fieldBC[2] == FBC_NEUMANN)
    {

        fbcValue[2] = 0.0; // symmetric BC, zero-Neumann

        for (i = 0; i < nx; i++)
        {
            long nid = ID(i, j, k);

            if (is_fixedBC(meshFunc[nid]))
                continue;

            meshFunc[nid] -= 2;
        }
    }
    else
    {
        esPIC_exit("Lower boundary can only be zero-Neumann");
    }

    /* upper boundary */
    j = ny - 1;
    if (fieldBC[3] == FBC_DIRICHLET)
    {
        for (i = 0; i < nx; i++)
        {
            meshFunc[ID(i, j, k)] = 1;
            fixedBV.insert(make_pair(ID(i, j, k), fbcValue[3]));
        }
    }
    else if (fieldBC[3] == FBC_NEUMANN)
    {
        for (i = 0; i < nx; i++)
        {
            long nid = ID(i, j, k);

            if (is_fixedBC(meshFunc[nid]))
                continue;

            meshFunc[nid] -= 2;
        }
    }
    else if (fieldBC[3] == FBC_PERIODIC)
    {
        esPIC_exit("Upper boundary can not be periodic");
    }
    cout << "# of grid poits in (x, y, z) direction = (" << nx << ", " << ny << ", " << nz << ")" << endl;

    vector<Tp> stiff_coeff;

    for (k = 0; k < nz; k++)
    {
        for (j = 0; j < ny; j++)
        {
            for (i = 0; i < nx; i++)
            {
                long row = ID(i, j, k), col[5] = {0, 0, 0, 0, 0};
                double f[5] = {0.0, 0.0, 0.0, 0.0, 0.0};

                if (is_fixedBC(meshFunc[row]))
                {
                    stiff_coeff.push_back(Tp(row, row, 1.0));
                    continue;
                }
                else
                {
                    innerNodes(i, j, k, ngpt, dx, col, f);
                }

                if (meshFunc[row] < 0)
                {
                    int flag = bcNeumannPeriodic(i, j, k, ngpt, dx, meshFunc[row], col, f);

                    if (flag != 0)
                        esPIC_exit("Boundary assemble incorrect");
                }

                for (int a = 0; a < 5; a++)
                {
                    if (col[a] < 0)
                        continue;

                    stiff_coeff.push_back(Tp(row, col[a], f[a]));
                }
            }
        }
    }
    //  k = 0;
    //  for (j = ny-1; j >= 0; j--) {
    //    i = 0;
    //    cout << meshFunc[ID(i,j,k)];
    //    for (i = 1; i < nx; i++)
    //      printf(" %5d", meshFunc[ID(i,j,k)]);
    //    cout << "\n";
    //  }

#undef ID

    A.setFromTriplets(stiff_coeff.begin(), stiff_coeff.end());

    /*
      ofstream ofile("matA.dat");

      for (int k = 0; k < A.outerSize(); k++)
        for (SpMat::InnerIterator it(A, k); it; ++it) {
          ofile << "(" << it.row() << ", " << it.col() << ") = "
            << it.value() << "\n";
        }
      ofile.close();
      ofile.clear();
    */

    delete[] meshFunc;
}

/* Setup the rhs vector for linear field solver */
void setup_rhs_linearFDS(const int ngpt[3], Vec &b)
{
}
