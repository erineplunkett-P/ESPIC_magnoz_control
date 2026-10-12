#include "es_PIC.h"
#include "global.h"

/* Update the rhs vector for linear poisson field solver */
void update_rhs_linearFDS(const int ngpt[3],
                          const int nspec,
                          const MapF &fixedBV,
                          double **chgDenst,
                          Vec &b)
{
    int nx = ngpt[0], ny = ngpt[1], nz = ngpt[2];
    size_t nxy = nx * ny;
    size_t nxyz = nxy * nz;

    for (size_t nid = 0; nid != nxyz; nid++)
        b(nid) = chgDenst[nid][nspec];

#define ID(i, j, k) (((k) * (ny) + (j)) * (nx) + (i))

    int i, j, k;
    for (k = 0; k < nz; k++)
    {
        /*----------------------------------------------*/
        /*-------- handling boundary conditions --------*/
        /*----------------------------------------------*/

        /*--------- left-side boundary (i = 0) ---------*/
        i = 0;
        if (fieldBC[0] == FBC_NEUMANN)
        {
            for (j = 0; j < ny; j++)
            {
                size_t nid = ID(i, j, k);
                b(nid) -= 2.0 * dxinv[0] * fbcValue[0];
            }
        }

        /*------- right-side boundary (i = nx-1) -------*/
        i = nx - 1;
        if (fieldBC[1] == FBC_NEUMANN)
        {
            for (j = 0; j < ny; j++)
            {
                size_t nid = ID(i, j, k);
                b(nid) += 2.0 * dxinv[0] * fbcValue[1];
            }
        }

        /*-------- lower-side boundary (j = 0) --------*/
        j = 0;
        // axisymmetric - zero Neumann boundary (fbcValue[2] = 0.0) so do not update
        //    if(fieldBC[2] == FBC_NEUMANN) {
        //      for(i = 0; i < nx; i++) {
        //        size_t nid = ID(i,j,k);
        //        b(nid) -= 2.0*dxinv[1]*fbcValue[2];
        //      }
        //    }

        /*------- upper-side boundary (j = ny-1) -------*/
        j = ny - 1;
        if (fieldBC[3] == FBC_NEUMANN)
        {
            for (i = 0; i < nx; i++)
            {
                size_t nid = ID(i, j, k);
                // axisymmetric
                b(nid) += 2.0 * dxinv[1] * fbcValue[3] + fbcValue[3] / (j * dx[1]);
            }
        }

        /*-------- handling fixed phi points --------*/
        for (MapF::const_iterator it = fixedBV.begin(); it != fixedBV.end(); it++)
        {
            b(it->first) = it->second;
        }
    } // end for(k = 0; k < nz; k++)

#undef ID
}

const double phi_ref = 0.0;

/* Update the rhs vector for non-linear poisson field solver */
void update_matrix_nonlinearFDS(const int ngpt[3],
                                const MapF &fixedBV,
                                Vec &p,
                                SpMat &A,
                                SpMat &J)
{
    int nx = ngpt[0], ny = ngpt[1], nz = ngpt[2];
    size_t nxy = nx * ny;
    size_t nxyz = nxy * nz;
    MapF::const_iterator it;

    for (size_t nid = 0; nid < nxyz; nid++)
    {
        /*-------- handling unfixed phi points --------*/
        it = fixedBV.find(nid);
        if (it != fixedBV.end())
            continue;
        // phi_ref hard coded here
        J.coeffRef(nid, nid) = A.coeffRef(nid, nid) + exp(p(nid) - phi_ref);
    }
}

/* Update the rhs vector for non-linear poisson field solver */
void update_rhs_nonlinearFDS(const int ngpt[3],
                             const MapF &fixedBV,
                             const Vec &b,
                             Vec &p,
                             SpMat &A,
                             Vec &f)
{
    int nx = ngpt[0], ny = ngpt[1], nz = ngpt[2];
    size_t nxy = nx * ny;
    size_t nxyz = nxy * nz;
    MapF::const_iterator it;

    f = A * p - b;
    for (size_t nid = 0; nid < nxyz; nid++)
    {
        /*-------- handling unfixed phi points --------*/
        it = fixedBV.find(nid);
        if (it != fixedBV.end())
            continue;
        f(nid) += exp(p(nid) - phi_ref); // phi_ref hard coded here
    }
}
