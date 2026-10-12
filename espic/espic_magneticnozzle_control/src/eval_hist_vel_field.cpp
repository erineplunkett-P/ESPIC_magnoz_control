#include "es_PIC.h"

void evalHistVelField(int nspec,
                      const long *nparts,
                      const int ngpt[3],
                      const double (*part)[MAXNPARTPROPS],
                      vector<double> &histVelField)
{
    int nx = ngpt[0], ny = ngpt[1], nz = ngpt[2];
    int nxy = nx * ny, nxyz = nxy * nz;
    double *vel_field = new double[nspec * nxyz * 6];
    double *nsum = new double[nspec * nxyz];

    memset(vel_field, 0, nspec * nxyz * 6 * sizeof(double));
    memset(nsum, 0, nspec * nxyz * sizeof(double));

    long npartsTot = nparts[MAXNSPECTYPES];

#define igpt(i, j, k) ((k) * (nxy) + (j) * (nx) + (i))

    double partX[3], partV[3];
    double xpm, ypm;
    double ddx0, ddy0, ddx1, ddy1, w;
    double r2, rj2, rjp2, vinv; // axisymmetric
    int nid, ienty;
    int ix, iy, iz = 0;

    for (long ipart = 0; ipart < npartsTot; ipart++)
    {
        int itype = (int)(part[ipart][ptype] + smallValue);

        for (int ip = 0; ip < 3; ip++)
        {
            partX[ip] = part[ipart][xp + ip];
            partV[ip] = part[ipart][up + ip];
        }

        /* in case very few particles move more than one cell per step. */
        /* in this case the time step dt should be reduced. */
        if (partX[0] >= domain[1])
            partX[0] = domain[1] - smallValue;
        if (partX[1] >= domain[3])
            partX[1] = domain[3] - smallValue;
        if (partX[2] >= domain[5])
            partX[2] = domain[5] - smallValue;
        if (partX[0] <= domain[0])
            partX[0] = domain[0] + smallValue;
        if (partX[1] <= domain[2])
            partX[1] = domain[2] + smallValue;
        if (partX[2] <= domain[4])
            partX[2] = domain[4] + smallValue;

        xpm = (partX[0] - domain[0]) * dxinv[0];
        ypm = (partX[1] - domain[2]) * dxinv[1];
        ix = (int)xpm;
        iy = (int)ypm;

        r2 = SQR(partX[1]);
        rj2 = iy * dx[1] + domain[2];
        rjp2 = rj2 + dx[1];
        rj2 *= rj2;
        rjp2 *= rjp2;
        vinv = PI * (rjp2 - rj2) * dx[0];
        vinv = 1. / vinv;

        ddx0 = xpm - ix;
        ddx1 = 1.0 - ddx0;
        ddx0 *= dx[0];
        ddx1 *= dx[0];
        ddy0 = PI * (r2 - rj2);
        ddy1 = PI * (rjp2 - r2);

        nid = igpt(ix, iy, iz);
        w = ddx1 * ddy1 * vinv;
        nsum[nid * nspec + itype] += w;
        ienty = 6 * nspec * nid + 6 * itype;
        for (int ip = 0; ip < 3; ip++)
        {
            vel_field[ienty + ip] += (w * partV[ip]); // review for electron axial and total velo. 
            vel_field[ienty + ip + 3] += (w * partV[ip] * partV[ip]);
        }

        nid = igpt(ix + 1, iy, iz);
        w = ddx0 * ddy1 * vinv;
        nsum[nid * nspec + itype] += w;
        ienty = 6 * nspec * nid + 6 * itype;
        for (int ip = 0; ip < 3; ip++)
        {
            vel_field[ienty + ip] += (w * partV[ip]);
            vel_field[ienty + ip + 3] += (w * partV[ip] * partV[ip]);
        }

        nid = igpt(ix, iy + 1, iz);
        w = ddx1 * ddy0 * vinv;
        nsum[nid * nspec + itype] += w;
        ienty = 6 * nspec * nid + 6 * itype;
        for (int ip = 0; ip < 3; ip++)
        {
            vel_field[ienty + ip] += (w * partV[ip]);
            vel_field[ienty + ip + 3] += (w * partV[ip] * partV[ip]);
        }

        nid = igpt(ix + 1, iy + 1, iz);
        w = ddx0 * ddy0 * vinv;
        nsum[nid * nspec + itype] += w;
        ienty = 6 * nspec * nid + 6 * itype;
        for (int ip = 0; ip < 3; ip++)
        {
            vel_field[ienty + ip] += (w * partV[ip]);
            vel_field[ienty + ip + 3] += (w * partV[ip] * partV[ip]);
        }
    } // for(long ipart = 0; ipart < npartsTot; ipart++)

#undef igpt

    for (int nid = 0; nid < nxyz; nid++)
    {
        for (int ispec = 0; ispec < nspec; ispec++)
        {
            double ns = nsum[nid * nspec + ispec];

            if (ns > 0.0)
            {
                for (int i = 0; i < 6; i++)
                {
                    int ienty = 6 * nspec * nid + 6 * ispec + i;
                    vel_field[ienty] /= ns;
                }
            }
        }
    }

    for (int ienty = 0; ienty < 6 * nspec * nxyz; ienty++)
        histVelField.push_back(vel_field[ienty]);

    delete[] nsum;
    delete[] vel_field;
}

void evalNewHistVelField(int nspec,
                         const long *nparts,
                         const int ngpt[3],
                         const double (*part)[MAXNPARTPROPS],
                         vector<double> &new_histVelField)
{
    int nx = ngpt[0], ny = ngpt[1], nz = ngpt[2];
    int nxy = nx * ny, nxyz = nxy * nz;
    double *vel_field = new double[nspec * nxyz * 5];
    double *nsum = new double[nspec * nxyz];

    memset(vel_field, 0, nspec * nxyz * 5 * sizeof(double));
    memset(nsum, 0, nspec * nxyz * sizeof(double));

    long npartsTot = nparts[MAXNSPECTYPES];

#define igpt(i, j, k) ((k) * (nxy) + (j) * (nx) + (i))

    double partX[3], partV[3];
    double xpm, ypm;
    double ddx0, ddy0, ddx1, ddy1, w;
    double r2, rj2, rjp2, vinv; // axisymmetric
    int nid, ienty;
    int ix, iy, iz = 0;
    for (long ipart = 0; ipart < npartsTot; ipart++)
    {
        int itype = (int)(part[ipart][ptype] + smallValue);
        for (int ip = 0; ip < 3; ip++)
        {
            partX[ip] = part[ipart][xp + ip];
            partV[ip] = part[ipart][up + ip];
        }

        /* in case very few particles move more than one cell per step. */
        /* in this case the time step dt should be reduced. */
        if (partX[0] >= domain[1])
            partX[0] = domain[1] - smallValue;
        if (partX[1] >= domain[3])
            partX[1] = domain[3] - smallValue;
        if (partX[2] >= domain[5])
            partX[2] = domain[5] - smallValue;
        if (partX[0] <= domain[0])
            partX[0] = domain[0] + smallValue;
        if (partX[1] <= domain[2])
            partX[1] = domain[2] + smallValue;
        if (partX[2] <= domain[4])
            partX[2] = domain[4] + smallValue;

        xpm = (partX[0] - domain[0]) * dxinv[0];
        ypm = (partX[1] - domain[2]) * dxinv[1];
        ix = (int)xpm;
        iy = (int)ypm;

        r2 = SQR(partX[1]);
        rj2 = iy * dx[1] + domain[2];
        rjp2 = rj2 + dx[1];
        rj2 *= rj2;
        rjp2 *= rjp2;
        vinv = PI * (rjp2 - rj2) * dx[0];
        vinv = 1. / vinv;

        ddx0 = xpm - ix;
        ddx1 = 1.0 - ddx0;
        ddx0 *= dx[0];
        ddx1 *= dx[0];
        ddy0 = PI * (r2 - rj2);
        ddy1 = PI * (rjp2 - r2);

        nid = igpt(ix, iy, iz);
        w = ddx1 * ddy1 * vinv;
        nsum[nid * nspec + itype] += w;
        ienty = 5 * nspec * nid + 5 * itype;
        // for (int ip = 0; ip < 5; ip++)
        // {
        // vel_field[ienty + ip] += (w * partV[ip]);
        // vel_field[ienty + ip + 3] += (w * partV[ip] * partV[ip]);
        // }

        vel_field[ienty + 0] += (w * partV[0] * partV[0] * partV[0]);
        vel_field[ienty + 1] += (w * partV[1] * partV[1] * partV[1]);
        vel_field[ienty + 2] += (w * partV[0] * partV[0] * partV[1]);
        vel_field[ienty + 3] += (w * partV[0] * partV[1] * partV[1]);
        vel_field[ienty + 4] += (w * partV[0] * partV[1]);

        nid = igpt(ix + 1, iy, iz);
        w = ddx0 * ddy1 * vinv;
        nsum[nid * nspec + itype] += w;
        ienty = 5 * nspec * nid + 5 * itype;
        // for (int ip = 0; ip < 3; ip++)
        // {
        // vel_field[ienty + ip] += (w * partV[ip]);
        // vel_field[ienty + ip + 3] += (w * partV[ip] * partV[ip]);
        // }

        vel_field[ienty + 0] += (w * partV[0] * partV[0] * partV[0]);
        vel_field[ienty + 1] += (w * partV[1] * partV[1] * partV[1]);
        vel_field[ienty + 2] += (w * partV[0] * partV[0] * partV[1]);
        vel_field[ienty + 3] += (w * partV[0] * partV[1] * partV[1]);
        vel_field[ienty + 4] += (w * partV[0] * partV[1]);

        nid = igpt(ix, iy + 1, iz);
        w = ddx1 * ddy0 * vinv;
        nsum[nid * nspec + itype] += w;
        ienty = 5 * nspec * nid + 5 * itype;
        // for (int ip = 0; ip < 3; ip++)
        // {
        // vel_field[ienty + ip] += (w * partV[ip]);
        // vel_field[ienty + ip + 3] += (w * partV[ip] * partV[ip]);
        // }

        vel_field[ienty + 0] += (w * partV[0] * partV[0] * partV[0]);
        vel_field[ienty + 1] += (w * partV[1] * partV[1] * partV[1]);
        vel_field[ienty + 2] += (w * partV[0] * partV[0] * partV[1]);
        vel_field[ienty + 3] += (w * partV[0] * partV[1] * partV[1]);
        vel_field[ienty + 4] += (w * partV[0] * partV[1]);

        nid = igpt(ix + 1, iy + 1, iz);
        w = ddx0 * ddy0 * vinv;
        nsum[nid * nspec + itype] += w;
        ienty = 5 * nspec * nid + 5 * itype;
        // for (int ip = 0; ip < 3; ip++)
        // {
        // vel_field[ienty + ip] += (w * partV[ip]);
        // vel_field[ienty + ip + 3] += (w * partV[ip] * partV[ip]);
        // }

        vel_field[ienty + 0] += (w * partV[0] * partV[0] * partV[0]);
        vel_field[ienty + 1] += (w * partV[1] * partV[1] * partV[1]);
        vel_field[ienty + 2] += (w * partV[0] * partV[0] * partV[1]);
        vel_field[ienty + 3] += (w * partV[0] * partV[1] * partV[1]);
        vel_field[ienty + 4] += (w * partV[0] * partV[1]);
    }

#undef igpt

    for (int nid = 0; nid < nxyz; nid++)
    {
        for (int ispec = 0; ispec < nspec; ispec++)
        {
            double ns = nsum[nid * nspec + ispec];

            if (ns > 0.0)
            {
                for (int i = 0; i < 5; i++)
                {
                    int ienty = 5 * nspec * nid + 5 * ispec + i;
                    vel_field[ienty] /= ns;
                }
            }
        }
    }

    for (int ienty = 0; ienty < 5 * nspec * nxyz; ienty++)
        new_histVelField.push_back(vel_field[ienty]);

    delete[] nsum;
    delete[] vel_field;
}

void evalNewNewHistVelField(int nspec,
                            const long *nparts,
                            const int ngpt[3],
                            const double (*part)[MAXNPARTPROPS],
                            vector<double> &new_new_histVelField)
{
    int nx = ngpt[0], ny = ngpt[1], nz = ngpt[2];
    int nxy = nx * ny, nxyz = nxy * nz;
    double *vel_field = new double[nspec * nxyz * 3];
    double *nsum = new double[nspec * nxyz];

    memset(vel_field, 0, nspec * nxyz * 3 * sizeof(double));
    memset(nsum, 0, nspec * nxyz * sizeof(double));

    long npartsTot = nparts[MAXNSPECTYPES];

#define igpt(i, j, k) ((k) * (nxy) + (j) * (nx) + (i))

    double partX[3], partV[3];
    double xpm, ypm;
    double ddx0, ddy0, ddx1, ddy1, w;
    double r2, rj2, rjp2, vinv; // axisymmetric
    int nid, ienty;
    int ix, iy, iz = 0;

    for (long ipart = 0; ipart < npartsTot; ipart++)
    {
        int itype = (int)(part[ipart][ptype] + smallValue);
        for (int ip = 0; ip < 3; ip++)
        {
            partX[ip] = part[ipart][xp + ip];
            partV[ip] = part[ipart][up + ip];
        }

        /* in case very few particles move more than one cell per step. */
        /* in this case the time step dt should be reduced. */
        if (partX[0] >= domain[1])
            partX[0] = domain[1] - smallValue;
        if (partX[1] >= domain[3])
            partX[1] = domain[3] - smallValue;
        if (partX[2] >= domain[5])
            partX[2] = domain[5] - smallValue;
        if (partX[0] <= domain[0])
            partX[0] = domain[0] + smallValue;
        if (partX[1] <= domain[2])
            partX[1] = domain[2] + smallValue;
        if (partX[2] <= domain[4])
            partX[2] = domain[4] + smallValue;

        xpm = (partX[0] - domain[0]) * dxinv[0];
        ypm = (partX[1] - domain[2]) * dxinv[1];
        ix = (int)xpm;
        iy = (int)ypm;

        r2 = SQR(partX[1]);
        rj2 = iy * dx[1] + domain[2];
        rjp2 = rj2 + dx[1];
        rj2 *= rj2;
        rjp2 *= rjp2;
        vinv = PI * (rjp2 - rj2) * dx[0];
        vinv = 1. / vinv;

        ddx0 = xpm - ix;
        ddx1 = 1.0 - ddx0;
        ddx0 *= dx[0];
        ddx1 *= dx[0];
        ddy0 = PI * (r2 - rj2);
        ddy1 = PI * (rjp2 - r2);

        nid = igpt(ix, iy, iz);
        w = ddx1 * ddy1 * vinv;
        nsum[nid * nspec + itype] += w;
        ienty = 3 * nspec * nid + 3 * itype;
        // for (int ip = 0; ip < 3; ip++)
        // {
        vel_field[ienty + 0] += (w * partV[0] * partV[0] * partV[0] * partV[0]);
        vel_field[ienty + 1] += (w * partV[1] * partV[1] * partV[1] * partV[1]);
        vel_field[ienty + 2] += (w * partV[0] * partV[0] * partV[1] * partV[1]);
        // vel_field[ienty + ip + 3] += (ddx1 * ddy1 * partV[ip] * partV[ip]);
        // }

        nid = igpt(ix + 1, iy, iz);
        w = ddx0 * ddy1 * vinv;
        nsum[nid * nspec + itype] += w;
        ienty = 3 * nspec * nid + 3 * itype;
        // for (int ip = 0; ip < 3; ip++)
        // {
        vel_field[ienty + 0] += (w * partV[0] * partV[0] * partV[0] * partV[0]);
        vel_field[ienty + 1] += (w * partV[1] * partV[1] * partV[1] * partV[1]);
        vel_field[ienty + 2] += (w * partV[0] * partV[0] * partV[1] * partV[1]);
        // vel_field[ienty + ip + 3] += (ddx0 * ddy1 * partV[ip] * partV[ip]);
        // }

        nid = igpt(ix, iy + 1, iz);
        w = ddx1 * ddy0 * vinv;
        nsum[nid * nspec + itype] += w;
        ienty = 3 * nspec * nid + 3 * itype;
        // for (int ip = 0; ip < 3; ip++)
        // {
        vel_field[ienty + 0] += (w * partV[0] * partV[0] * partV[0] * partV[0]);
        vel_field[ienty + 1] += (w * partV[1] * partV[1] * partV[1] * partV[1]);
        vel_field[ienty + 2] += (w * partV[0] * partV[0] * partV[1] * partV[1]);
        // vel_field[ienty + ip + 3] += (ddx1 * ddy0 * partV[ip] * partV[ip]);
        // }

        nid = igpt(ix + 1, iy + 1, iz);
        w = ddx0 * ddy0 * vinv;
        nsum[nid * nspec + itype] += w;
        ienty = 3 * nspec * nid + 3 * itype;

        vel_field[ienty + 0] += (w * partV[0] * partV[0] * partV[0] * partV[0]);
        vel_field[ienty + 1] += (w * partV[1] * partV[1] * partV[1] * partV[1]);
        vel_field[ienty + 2] += (w * partV[0] * partV[0] * partV[1] * partV[1]);
        // vel_field[ienty + ip + 3] += (ddx0 * ddy0 * partV[ip] * partV[ip]);
        // }

    } // for(long ipart = 0; ipart < npartsTot; ipart++)

#undef igpt

    for (int nid = 0; nid < nxyz; nid++)
    {
        for (int ispec = 0; ispec < nspec; ispec++)
        {
            double ns = nsum[nid * nspec + ispec];

            if (ns > 0.0)
            {
                for (int i = 0; i < 3; i++)
                {
                    int ienty = 3 * nspec * nid + 3 * ispec + i;
                    vel_field[ienty] /= ns;
                }
            }
        }
    }

    for (int ienty = 0; ienty < 3 * nspec * nxyz; ienty++)
    {
        new_new_histVelField.push_back(vel_field[ienty]);
    }

    delete[] nsum;
    delete[] vel_field;
}