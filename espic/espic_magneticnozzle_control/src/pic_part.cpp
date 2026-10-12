#include "es_PIC.h"
#include <new>
#include <cstring>

using namespace std;

extern vector<double> fnp_res[MAXNFLOWS][6]; // (fractional) residual number of particles

int *vacMap;
double ***chgDenst_thrd;
double *volinv;

/* compare function used for qsort */
int compare(const void *a, const void *b)
{
    return *(int *)a - *(int *)b;
}

void cross_prod(double vecC[3], double vecA[3], double vecB[3])
{
    double temp[3];

    temp[0] = vecA[1] * vecB[2] - vecA[2] * vecB[1];
    temp[1] = vecA[2] * vecB[0] - vecA[0] * vecB[2];
    temp[2] = vecA[0] * vecB[1] - vecA[1] * vecB[0];

    for (int a = 0; a < 3; a++)
        vecC[a] = temp[a];
}

void adjustPartMem()
{
    const int oldMaxParts = MAXNPARTS;
    MAXNPARTS += MAXNPARTSINCR;
    cout << "Adjusting particle memory to (MAXNPARTS = " << MAXNPARTS << ")" << endl;

    double(*ptr_part_new)[MAXNPARTPROPS] = nullptr;
    int *ptr_vacMap_new = nullptr;

    try
    {
        ptr_part_new = new double[MAXNPARTS][MAXNPARTPROPS];
        ptr_vacMap_new = new int[MAXNPARTS];
    }
    catch (const std::bad_alloc &e)
    {
        esPIC_exit(e.what());
    }

    // copy old contents and zero-fill the newly added region
    std::memcpy(ptr_part_new, part, static_cast<size_t>(oldMaxParts) * MAXNPARTPROPS * sizeof(double));
    std::memcpy(ptr_vacMap_new, vacMap, static_cast<size_t>(oldMaxParts) * sizeof(int));
    std::memset(ptr_part_new + oldMaxParts, 0, static_cast<size_t>(MAXNPARTS - oldMaxParts) * MAXNPARTPROPS * sizeof(double));
    std::memset(ptr_vacMap_new + oldMaxParts, 0, static_cast<size_t>(MAXNPARTS - oldMaxParts) * sizeof(int));

    delete[] part;
    delete[] vacMap;
    part = ptr_part_new;
    vacMap = ptr_vacMap_new;
}

/* Allocate memory for part; initialize particle # counter (*parts); load ambient particles if input file [ambient.part] exists. */
bool initPart(const int flag_read_restart, long *nparts, double *partWeight)
{
    memset(nparts, 0, (MAXNSPECTYPES + 1) * sizeof(long));

    // c++ style
    // note: it usually does NOT throw exception on (debian) linux due to memory overcommit mechanism
    try
    {
        part = new double[MAXNPARTS][MAXNPARTPROPS];
        vacMap = new int[MAXNPARTS];
    }
    catch (const std::bad_alloc &e)
    {
        esPIC_exit(e.what());
    }
    //  if(part   == nullptr) esPIC_exit("not enough memory for particle allocation");
    //  if(vacMap == nullptr) esPIC_exit("not enough memory for particle allocation");

    // c style
    //  part   = (double (*)[MAXNPARTPROPS])malloc(MAXNPARTS*MAXNPARTPROPS*sizeof(double));
    //  vacMap = (int *)malloc(MAXNPARTS*sizeof(int));
    //  part   = (double (*)[MAXNPARTPROPS])calloc(MAXNPARTS*MAXNPARTPROPS, sizeof(double));
    //  vacMap = (int *)calloc(MAXNPARTS, sizeof(int));

    //  if(part   == nullptr) esPIC_exit("not enough memory for particle allocation");
    //  if(vacMap == nullptr) esPIC_exit("not enough memory for particle allocation");
    //  memset(part,   0, MAXNPARTS*MAXNPARTPROPS*sizeof(double));

    if (!is_file_exist("weight.part"))
    {
        esPIC_exit("[weight.part] does NOT exist");
    }
    else
    {
        ifstream infile("weight.part");
        if (!infile)
            esPIC_exit("Cannot read [weight.part]");

        cout << "Reading weights of ambient particles from [weight.part]..." << endl;

        string line;
        int itype = 0;
        getline(infile, line); // skip 1st line which is the comment
        while (getline(infile, line))
        {
            istringstream iss(line);

            if (!(iss >> partWeight[itype]))
                break;

            ++itype;
        }
        cout << "Weights of " << itype << " species have been read...\n";

        infile.close();
        infile.clear();
    }

    if (flag_read_restart)
        return false;
    if (!is_file_exist("ambient.part"))
        return false;

    ifstream infile("ambient.part");
    cout << "Loading ambient particles from [ambient.part]..." << endl;

    if (!infile)
        esPIC_exit("Cannot read [ambient.part]");

    /* Give particles their initial positions, velocities and species types */
    int ipart = 0;
    string line;
    getline(infile, line); // skip 1st line
    while (getline(infile, line))
    {
        istringstream iss(line);

        for (int i = 0; i < partPropN; ++i)
            part[ipart][i] = 0.0;

        if (!(iss >> part[ipart][xp] >> part[ipart][yp] >> part[ipart][zp] >> part[ipart][up] >> part[ipart][vp] >> part[ipart][wp] >> part[ipart][ptype]))
            break;
        part[ipart][pid] = -1; // group ID number

        ++nparts[(int)(part[ipart][ptype] + 1e-6)];
        if (++ipart == MAXNPARTS)
            adjustPartMem();
    }
    nparts[MAXNSPECTYPES] = ipart;
    //  cout << "nparts[0, 1, end] = " << "(" << nparts[0] << ", " << nparts[1] << ", " << nparts[MAXNSPECTYPES] << ")" << endl;

    infile.close();
    infile.clear();

    return true;
}

void bilinear_weighting(double px[3], double domain[6], int &ix, int &iy, double f[4])
{
    double xpm, ypm;
    double ddx0, ddx1, ddy0, ddy1;

    /* in case very few particles move more than one cell per step. */
    /* in this case the time step dt should be reduced. */
    if (px[0] >= domain[1])
        px[0] = domain[1] - smallValue;
    if (px[1] >= domain[3])
        px[1] = domain[3] - smallValue;
    if (px[2] >= domain[5])
        px[2] = domain[5] - smallValue;
    if (px[0] < domain[0])
        px[0] = domain[0];
    if (px[1] < domain[2])
        px[1] = domain[2];
    if (px[2] < domain[4])
        px[2] = domain[4];

    xpm = (px[0] - domain[0]) * dxinv[0];
    ypm = (px[1] - domain[2]) * dxinv[1];
    ix = (int)xpm;
    iy = (int)ypm;

    ddx0 = xpm - ix;
    ddx1 = 1.0 - ddx0;
    ddy0 = ypm - iy;
    ddy1 = 1.0 - ddy0;

    f[0] = ddx1 * ddy1;
    f[1] = ddx0 * ddy1;
    f[2] = ddx1 * ddy0;
    f[3] = ddx0 * ddy0;
}

void area_weighting(double px[3], double domain[6], int &ix, int &iy, double f[4])
{
    double xpm, ypm;
    double ddx0, ddx1, ddy0, ddy1;
    double rj2, rjp2, r2, vinv; // axisymmetric

    /* in case very few particles move more than one cell per step. */
    /* in this case the time step dt should be reduced. */
    if (px[0] >= domain[1])
        px[0] = domain[1] - smallValue;
    if (px[1] >= domain[3])
        px[1] = domain[3] - smallValue;
    if (px[2] >= domain[5])
        px[2] = domain[5] - smallValue;
    if (px[0] < domain[0])
        px[0] = domain[0];
    if (px[1] < domain[2])
        px[1] = domain[2];
    if (px[2] < domain[4])
        px[2] = domain[4];

    xpm = (px[0] - domain[0]) * dxinv[0];
    ypm = (px[1] - domain[2]) * dxinv[1];
    ix = (int)xpm;
    iy = (int)ypm;

    r2 = SQR(px[1]);
    rj2 = iy * dx[1] + domain[2];
    rjp2 = rj2 + dx[1];
    rj2 *= rj2;
    rjp2 *= rjp2;
    vinv = PI * (rjp2 - rj2) * dx[0];
    vinv = 1.0 / vinv;

    ddx0 = xpm - ix;
    ddx1 = 1.0 - ddx0;
    ddx0 *= dx[0];
    ddx1 *= dx[0];
    ddy0 = PI * (r2 - rj2);
    ddy1 = PI * (rjp2 - r2);

    f[0] = ddx1 * ddy1 * vinv;
    f[1] = ddx0 * ddy1 * vinv;
    f[2] = ddx1 * ddy0 * vinv;
    f[3] = ddx0 * ddy0 * vinv;
}

void field2part(double partp[MAXNPARTPROPS], const int ngpt[3],
                double Ep[3], double Bp[3],
                double (*E)[3], double (*B)[3])
{
    int nx = ngpt[0], ny = ngpt[1];
    int nxy = nx * ny;
    int ix, iy;
    double f[4];
    //  double xpm, ypm;    // position of particle within a mesh(cell)
    //  double ddx0, ddy0;
    //  double ddx1, ddy1;
    //
#define igpt(i, j, k) ((k) * (nxy) + (j) * (nx) + (i))
    //
    //  if(partp[xp] >= domain[1]) partp[xp] = domain[1] - 1e-8;
    //  if(partp[yp] >= domain[3]) partp[yp] = domain[3] - 1e-8;
    //
    //  xpm = (partp[xp] - domain[0])*dxinv[0];
    //  ypm = (partp[yp] - domain[2])*dxinv[1];
    //  ix  = (int)xpm; ddx0 = xpm - ix;
    //  iy  = (int)ypm; ddy0 = ypm - iy;
    //  ddx1 = 1.0 - ddx0;
    //  ddy1 = 1.0 - ddy0;
    //
    //  /* 2D interpolation */
    //  /* fields from mesh nodes to particles */
    //  for(int a = 0; a < 2; a++) {
    //    Ep[a] = E[igpt(ix,iy,0)][a]    *ddx1*ddy1
    //          + E[igpt(ix+1,iy,0)][a]  *ddx0*ddy1
    //          + E[igpt(ix,iy+1,0)][a]  *ddx1*ddy0
    //          + E[igpt(ix+1,iy+1,0)][a]*ddx0*ddy0;
    //
    //    Bp[a] = B[igpt(ix,iy,0)][a]    *ddx1*ddy1
    //          + B[igpt(ix+1,iy,0)][a]  *ddx0*ddy1
    //          + B[igpt(ix,iy+1,0)][a]  *ddx1*ddy0
    //          + B[igpt(ix+1,iy+1,0)][a]*ddx0*ddy0;
    //  }
    //  Ep[2] = Bp[2] = 0.0;
    //

    area_weighting(partp, domain, ix, iy, f);

    /* 2D interpolation */
    /* fields from mesh nodes to particles */
    for (int a = 0; a < 2; a++)
    {
        Ep[a] = E[igpt(ix, iy, 0)][a] * f[0] + E[igpt(ix + 1, iy, 0)][a] * f[1] + E[igpt(ix, iy + 1, 0)][a] * f[2] + E[igpt(ix + 1, iy + 1, 0)][a] * f[3];

        Bp[a] = B[igpt(ix, iy, 0)][a] * f[0] + B[igpt(ix + 1, iy, 0)][a] * f[1] + B[igpt(ix, iy + 1, 0)][a] * f[2] + B[igpt(ix + 1, iy + 1, 0)][a] * f[3];
    }
    Ep[2] = Bp[2] = 0.0;

#undef igpt
}

//{
//  int nx = ngpt[0], ny = ngpt[1];
//  int nxy = nx*ny;
//  int ix, iy, iz;
//  double xpm, ypm, zpm;   // position of particle within a mesh(cell)
//  double ddx, ddy, ddz;
//  double f, g;
//
//  #define igpt(i,j,k) ((k)*(nxy) + (j)*(nx) + (i))
//
//  xpm = (partp[xp] - domain[0])*dxinv[0];
//  ypm = (partp[yp] - domain[2])*dxinv[1];
//  zpm = (partp[zp] - domain[4])*dxinv[2];
//  ix  = (int)xpm; ddx = xpm - ix;
//  iy  = (int)ypm; ddy = ypm - iy;
//  iz  = (int)zpm; ddz = zpm - iz;
//
//  /* tri-linear field interpolations (linear in x by linear in y by linear in z) */
//
//  /* E field interpolation */
//
//  f = E[igpt(ix,iy,iz)][0] + ddx*(E[igpt(ix+1,iy,iz)][0] - E[igpt(ix,iy,iz)][0]);
//  f = f + ddy*(E[igpt(ix,iy+1,iz)][0]
//        + ddx*(E[igpt(ix+1,iy+1,iz)][0] - E[igpt(ix,iy+1,iz)][0]) - f);
//  g = E[igpt(ix,iy,iz+1)][0] + ddx*(E[igpt(ix+1,iy,iz+1)][0] - E[igpt(ix,iy,iz+1)][0]);
//  g = g + ddy*(E[igpt(ix,iy+1,iz+1)][0]
//        + ddx*(E[igpt(ix+1,iy+1,iz+1)][0] - E[igpt(ix,iy+1,iz+1)][0]) - g);
//  Ep[0] = f + ddz*(g - f);
//
//  f = E[igpt(ix,iy,iz)][1] + ddy*(E[igpt(ix,iy+1,iz)][1] - E[igpt(ix,iy,iz)][1]);
//  f = f + ddz*(E[igpt(ix,iy,iz+1)][1]
//        + ddy*(E[igpt(ix,iy+1,iz+1)][1] - E[igpt(ix,iy,iz+1)][1]) - f);
//  g = E[igpt(ix+1,iy,iz)][1] + ddy*(E[igpt(ix+1,iy+1,iz)][1] - E[igpt(ix+1,iy,iz)][1]);
//  g = g + ddz*(E[igpt(ix+1,iy,iz+1)][1]
//        + ddy*(E[igpt(ix+1,iy+1,iz+1)][1] - E[igpt(ix+1,iy,iz+1)][1]) - g);
//  Ep[1] = f + ddx*(g - f);
//
//  f = E[igpt(ix,iy,iz)][2] + ddz*(E[igpt(ix,iy,iz+1)][2] - E[igpt(ix,iy,iz)][2]);
//  f = f + ddx*(E[igpt(ix+1,iy,iz)][2]
//        + ddz*(E[igpt(ix+1,iy,iz+1)][2] - E[igpt(ix+1,iy,iz)][2]) - f);
//  g = E[igpt(ix,iy+1,iz)][2] + ddz*(E[igpt(ix,iy+1,iz+1)][2] - E[igpt(ix,iy+1,iz)][2]);
//  g = g + ddx*(E[igpt(ix+1,iy+1,iz)][2]
//        + ddz*(E[igpt(ix+1,iy+1,iz+1)][2] - E[igpt(ix+1,iy+1,iz)][2]) - g);
//  Ep[2] = f + ddy*(g - f);
//
//  /* B field interpolation */
//
//  f = B[igpt(ix,iy,iz)][0] + ddx*(B[igpt(ix+1,iy,iz)][0] - B[igpt(ix,iy,iz)][0]);
//  f = f + ddy*(B[igpt(ix,iy+1,iz)][0]
//        + ddx*(B[igpt(ix+1,iy+1,iz)][0] - B[igpt(ix,iy+1,iz)][0]) - f);
//  g = B[igpt(ix,iy,iz+1)][0] + ddx*(B[igpt(ix+1,iy,iz+1)][0] - B[igpt(ix,iy,iz+1)][0]);
//  g = g + ddy*(B[igpt(ix,iy+1,iz+1)][0]
//        + ddx*(B[igpt(ix+1,iy+1,iz+1)][0] - B[igpt(ix,iy+1,iz+1)][0]) - g);
//  Bp[0] = f + ddz*(g - f);
//
//  f = B[igpt(ix,iy,iz)][1] + ddy*(B[igpt(ix,iy+1,iz)][1] - B[igpt(ix,iy,iz)][1]);
//  f = f + ddz*(B[igpt(ix,iy,iz+1)][1]
//        + ddy*(B[igpt(ix,iy+1,iz+1)][1] - B[igpt(ix,iy,iz+1)][1]) - f);
//  g = B[igpt(ix+1,iy,iz)][1] + ddy*(B[igpt(ix+1,iy+1,iz)][1] - B[igpt(ix+1,iy,iz)][1]);
//  g = g + ddz*(B[igpt(ix+1,iy,iz+1)][1]
//        + ddy*(B[igpt(ix+1,iy+1,iz+1)][1] - B[igpt(ix+1,iy,iz+1)][1]) - g);
//  Bp[1] = f + ddx*(g - f);
//
//  f = B[igpt(ix,iy,iz)][2] + ddz*(B[igpt(ix,iy,iz+1)][2] - B[igpt(ix,iy,iz)][2]);
//  f = f + ddx*(B[igpt(ix+1,iy,iz)][2]
//        + ddz*(B[igpt(ix+1,iy,iz+1)][2] - B[igpt(ix+1,iy,iz)][2]) - f);
//  g = B[igpt(ix,iy+1,iz)][2] + ddz*(B[igpt(ix,iy+1,iz+1)][2] - B[igpt(ix,iy+1,iz)][2]);
//  g = g + ddx*(B[igpt(ix+1,iy+1,iz)][2]
//        + ddz*(B[igpt(ix+1,iy+1,iz+1)][2] - B[igpt(ix+1,iy+1,iz)][2]) - g);
//  Bp[2] = f + ddy*(g - f);
//
//  #undef  igpt
//}

bool appPartBC(double partX[3], double partV[3], int ispec, int ibound_out)
{
    bool partAbsorbed = false;
    double domainL = 0.0;
    int dim = -1, dim2 = -1;

    // int pbcType = partBC[ibound_out];
    int pbcType = -1000;

    if (ispec == 0)
    {
        pbcType = partBC[ibound_out];
    }
    else if (ispec == 1)
    {
        pbcType = partBC_new[ibound_out];
    }
    else
    {
        esPIC_exit("We only support 2 species for now");
    }

    switch (pbcType)
    {
    case PBC_ABSORPT: // absorption BC
        partAbsorbed = true;
        break;

    case PBC_REFLECT: // reflect BC
        dim = ibound_out / 2;
        partX[dim] = 2 * domain[ibound_out] - partX[dim];
        partV[dim] = -partV[dim];
        break;

    case PBC_PERIODIC: // perodic BC
        dim = ibound_out / 2;
        dim2 = dim * 2;
        domainL = domain[dim2 + 1] - domain[dim2];
        if (dim2 == ibound_out)
            partX[dim] += domainL; // partX[dim] < domain_max
        else
            partX[dim] -= domainL;
        break;
    }

    return partAbsorbed;
}

bool appPartObjBC(double partX[3], const ObjectList &objlist)
{
    bool flag = false;

    ObjectList::size_type nobjs = objlist.size();

    for (ObjectList::size_type iobj = 0; iobj < nobjs; iobj++)
    {
        if (!objlist[iobj].realObj)
            continue;

        flag = is_inside_object(partX, objlist[iobj]);

        if (flag)
            return flag;
    }

    return flag;
}

void partMove(const double half_dt,
              const int ngpt[3],
              long *nparts,
              double (*part)[MAXNPARTPROPS],
              double (*E)[3], double (*B)[3],
              const ObjectList &objlist,
              Sources &sourcelist,
              const int nthrds)
{
    double Ep[3], Bp[3];
    double QM;
    double t0, t[3], s[3], QME_hfdt[3], Uhf[3];
    double dt = 2.0 * half_dt;
    long npartsTot = nparts[MAXNSPECTYPES];
    long nparts_per_thrd = npartsTot / nthrds;
    long thrd_start[MAXNTHRDS], thrd_end[MAXNTHRDS];
    long nvacs = 0, nvacs_thrd[MAXNTHRDS] = {0};

    for (int ithrd = 0; ithrd < nthrds; ithrd++)
    {
        thrd_start[ithrd] = nparts_per_thrd * ithrd;
        thrd_end[ithrd] = nparts_per_thrd * (ithrd + 1);
    }
    thrd_end[nthrds - 1] = npartsTot;

#pragma omp parallel default(none) private(QM, Ep, Bp, t0, t, s, QME_hfdt, Uhf)      \
    shared(thrd_start, thrd_end, part, spec, ngpt, domain,                           \
               objlist, partBC, partBC_new, fnp_res, sourcelist, vacMap, nvacs_thrd, \
               E, B, dx, dt, half_dt)
    {
        int tid = omp_get_thread_num();
        double y2, z2, r1, r2, sina, cosa, vr1, vt1, vr2, vt2;
        double old_r, old_x;

        for (long ipart = thrd_start[tid]; ipart < thrd_end[tid]; ipart++)
        {
            int ispec = (long)(part[ipart][ptype] + 1e-6);

            /* fixed particles, do not move */
            if (spec[ispec][specFixed] > 0.0)
                continue;

            QM = spec[ispec][specQM];
            old_x = (part[ipart] + xp)[0];
            old_r = (part[ipart] + xp)[1];

            field2part(part[ipart], ngpt, Ep, Bp, E, B); // field interpolation to part

            /* compute t and s vectors needed for rotation */
            //      Bm = sqrt(SQR(Bp[0]) + SQR(Bp[1]) + SQR(Bp[2]));
            //      t0 = QM*Bm*half_dt;
            //      s0 = 2.0*t0/(1.0 + SQR(t0));
            //      for(int a = 0; a < 3; a++) {
            //        double Bp_unit = Bp[a]/(Bm + smallValue);
            //
            //        t[a] = t0*Bp_unit;
            //        s[a] = s0*Bp_unit;
            //        QME_hfdt[a] = QM*half_dt*Ep[a];
            //      }
            t0 = 0.0;
            for (int a = 0; a < 3; a++)
            {
                QME_hfdt[a] = QM * half_dt * Ep[a];
                t[a] = QM * half_dt * Bp[a];
                t0 += SQR(t[a]);
            }
            t0 = 2.0 / (1.0 + t0);
            for (int a = 0; a < 3; a++)
                s[a] = t0 * t[a];

            /* Update particle velocities using leap-frog scheme */
            for (int a = 0; a < 3; a++)
                Uhf[a] = part[ipart][a + up] + QME_hfdt[a];
            cross_prod(part[ipart] + up, Uhf, t);
            for (int a = 0; a < 3; a++)
                part[ipart][a + up] += Uhf[a];
            cross_prod(part[ipart] + up, part[ipart] + up, s);

            for (int a = 0; a < 3; a++)
            {
                part[ipart][a + up] += (Uhf[a] + QME_hfdt[a]);

                //        double dl = part[ipart][a+up]*dt;
                //        if(dl > dx[a]) {
                //          printf("->: Warning: part[%ld]'s V(%d) = %e (*dt = %e > dx), possibly moves > 1 cell <--\n", ipart, a, part[ipart][a+up], dl);
                //          cout << "-> Warning: part[" << ipart << "]V(" << a << ") = " << part[ipart][a+up] << " (*dt = " << dl << "), possibly move more than 1 cell <-" << endl;
                //        }
            }

            /* Update particle positions */
            //      for(int a = 0; a < 3; a++) part[ipart][a+xp] += dt*part[ipart][a+up];
            // axisymmetric
            r1 = part[ipart][yp];
            vr1 = part[ipart][vp];
            vt1 = part[ipart][wp];
            y2 = dt * vr1 + r1;
            z2 = dt * vt1;
            r2 = sqrt(SQR(y2) + SQR(z2));
            if (r2 > 0)
            {
                sina = z2 / r2;
                cosa = y2 / r2;
            }
            else
            {
                sina = 0.0;
                cosa = 1.0;
            }
            vr2 = cosa * vr1 + sina * vt1;
            vt2 = cosa * vt1 - sina * vr1;
            part[ipart][vp] = vr2;
            part[ipart][wp] = vt2;

            part[ipart][xp] += dt * part[ipart][up];
            part[ipart][yp] = r2;
            part[ipart][zp] = 0.0;

            //      /* 2D simulation - do not move particles in z direction */
            //      if(ngpt[2] == 1) part[ipart][zp] = 0.0;

            /* Judge whether need to apply particle BC */
            bool poutFlag[6];
            poutFlag[0] = part[ipart][xp] < domain[0];
            poutFlag[1] = part[ipart][xp] > domain[1];
            poutFlag[2] = part[ipart][yp] < domain[2];
            poutFlag[3] = part[ipart][yp] > domain[3];
            if (ngpt[2] == 1)
            { // 2D simulation - no particle out from z direction
                poutFlag[4] = poutFlag[5] = false;
            }
            else
            {
                poutFlag[4] = part[ipart][zp] < domain[4];
                poutFlag[5] = part[ipart][zp] > domain[5];
            }

            bool outofDomain = false, partAbsorbed = false, partObjAbsorbed = false;
            int ibound_out;

            for (ibound_out = 0; ibound_out < 6; ibound_out++)
            {
                outofDomain = outofDomain || poutFlag[ibound_out];
                if (outofDomain)
                    break;
            }

            if (outofDomain)
            // partAbsorbed = appPartBC(part[ipart] + xp, part[ipart] + up, ispec, ibound_out);
            {
                if (ibound_out == 0)
                {
                    int temp_pbctype = -1000;
                    int dim = -1;
                    int dim2 = -1;
                    double domainL = 0.0;
                    double *partX = part[ipart] + xp;
                    double *partV = part[ipart] + up;
                    partAbsorbed = false;    // reset the flag to make it safe
                    partObjAbsorbed = false; // reset the flag to make it safe
                    if (ispec == 0)
                    {
                        temp_pbctype = partBC[ibound_out];
                    }
                    else if (ispec == 1)
                    {
                        temp_pbctype = partBC_new[ibound_out];
                    }
                    else
                    {
                        esPIC_exit("We only support 2 species for now");
                    }
                    switch (temp_pbctype)
                    {
                    case PBC_ABSORPT: // absorption BC
                    {
                        partAbsorbed = true;
                        break;
                    }
                    case PBC_REFLECT: // reflect BC, revise to hard code the reflection plane, only those particles actually hit on the emission plane should be absorbed, assuming only has one conducting object sitting at the emission plane
                    {
                        // dim = ibound_out / 2;
                        // partX[dim] = 2 * domain[ibound_out] - partX[dim];
                        // partV[dim] = -partV[dim];
                        double emit_plane_upperbound = objlist[0].xmax[1];
                        double emit_plane_lowerbound = objlist[0].xmin[1];
                        double intersection = (old_r * partX[0] - partX[1] * old_x) / (partX[0] - old_x); // y=kx+b, b = (y0x1-y1x0)/(x1-x0)
                        // if (((partX[1] + old_y) / 2.0 < emit_plane_lowerbound) || ((partX[1] + old_y) / 2.0 > emit_plane_upperbound)) // the reflection point is not on the emission plane
                        if ((intersection < emit_plane_lowerbound) || (intersection > emit_plane_upperbound))
                        {
                            // reflect the particle
                            dim = ibound_out / 2;
                            partX[dim] = 2 * domain[ibound_out] - partX[dim];
                            partV[dim] = -partV[dim];
                        }
                        else // need to apply Obj boundary
                        {
                            if (appPartObjBC(part[ipart] + xp, objlist))
                            {
                                vacMap[thrd_start[tid] + (nvacs_thrd[tid]++)] = ipart;
                                partObjAbsorbed = true;
                            }
                            else // these particles are not in the conductor, but still techniqually needs to be absorbed as they hit the emission plane
                            {
                                partAbsorbed = true;
                            }
                        }
                        break;
                    }
                    case PBC_PERIODIC: // perodic BC
                    {
                        dim = ibound_out / 2;
                        dim2 = dim * 2;
                        domainL = domain[dim2 + 1] - domain[dim2];
                        if (dim2 == ibound_out)
                            partX[dim] += domainL; // partX[dim] < domain_max
                        else
                            partX[dim] -= domainL;
                        break;
                    }
                    }
                }
            }

            if (partAbsorbed)
                vacMap[thrd_start[tid] + (nvacs_thrd[tid]++)] = ipart;

            if (!partAbsorbed)
            {
                int isrc = sourcelist.collect_particle(part[ipart] + xp, (int)(part[ipart][ptype] + 0.0001));
                if (isrc >= 0)
                {
                    vacMap[thrd_start[tid] + (nvacs_thrd[tid]++)] = ipart;
#pragma omp critical
                    sourcelist.count_np_collect(isrc);
                }
                else if (appPartObjBC(part[ipart] + xp, objlist))
                {
                    if (partObjAbsorbed == false)
                    {
                        vacMap[thrd_start[tid] + (nvacs_thrd[tid]++)] = ipart;
                    }
                }
            } // end if(!partAbsorbed)

        } // for(int ipart = 0; ipart < npartsTot; ipart++)
    } // end #pragma omp parallel

    for (int ithrd = 0; ithrd < nthrds; ithrd++)
        nvacs += nvacs_thrd[ithrd];
    //  cout << npartsTot << ", " << nvacs << endl;

    if (nvacs > 0)
    {
        nparts[MAXNSPECTYPES] = npartsTot - nvacs;
        long insert_i = nvacs_thrd[0];
        for (int ithrd = 1; ithrd < nthrds; ithrd++)
        {
            long vac_start = thrd_start[ithrd];
            if (insert_i == vac_start)
                continue;
            long vac_end = vac_start + nvacs_thrd[ithrd];
            for (long ipart = vac_start; ipart < vac_end; ipart++)
                vacMap[insert_i++] = vacMap[ipart];
        }

        int vac_end_id = vacMap[nvacs - 1];
        for (int ivac = 0; ivac < nvacs; ivac++)
        {
            int vac_id = vacMap[ivac];
            while (--npartsTot == vac_end_id)
            {
                --nvacs;
                --nparts[(int)(part[vac_end_id][ptype] + 0.0001)];
                if (ivac == nvacs)
                    break;
                vac_end_id = vacMap[nvacs - 1];
            }
            if (ivac < nvacs)
            {
                --nparts[(int)(part[vac_id][ptype] + 0.0001)];
                for (int a = 0; a < partPropN; a++)
                    part[vac_id][a] = part[npartsTot][a];
            }
        }
        if (nparts[MAXNSPECTYPES] != npartsTot)
        {
            esPIC_exit("move vacancy");
        }
    } // if(nvacs > 0)
}

void meshChgBoundary(const int nid[2],
                     const int nspec,
                     int pbc_out_0[2],
                     int pbc_out_1[2],
                     double *partWeight,
                     double **chgDenst)
{
    int pbc[2];
    for (int ispec = 0; ispec < nspec; ispec++)
    {
        if (ispec == 0)
        {
            pbc[0] = pbc_out_0[0];
            pbc[1] = pbc_out_0[1];
        }
        else if (ispec == 1)
        {
            pbc[0] = pbc_out_1[0];
            pbc[1] = pbc_out_1[1];
        }
        else
        {
            esPIC_exit("Only support 2 species for now");
        }
        if ((pbc[0] == PBC_PERIODIC) || (pbc[1] == PBC_PERIODIC))
        {

            chgDenst[nid[0]][ispec] += chgDenst[nid[1]][ispec];
            chgDenst[nid[1]][ispec] = chgDenst[nid[0]][ispec];
        }
        else
        {
            if (pbc[0] == PBC_REFLECT || pbc[0] == PBC_ABSORPT)
                chgDenst[nid[0]][ispec] *= 2;

            if (pbc[1] == PBC_REFLECT || pbc[1] == PBC_ABSORPT)
                chgDenst[nid[1]][ispec] *= 2;
        }
    } // end for(int ispec = 0; ispec < nspec; ispec++)
}

/*
void meshChgObject(const int ngpt[3],
                   const int nspec,
                   double *partWeight,
                   double **chgDenst)
{
  int nx = ngpt[0], ny = ngpt[1];

  #define ID(i,j,k) (((k)*(ny) + (j))*(nx) + (i))

  for(int iobj = 0; iobj < nobjects; iobj++) {
    int istrt, iend, jstrt, jend;
    int i, j, k = 0;

    istrt = (int)((objectCoord[iobj][0] - domain[0])*dxinv[0] + 1e-6);
    iend  = (int)((objectCoord[iobj][1] - domain[0])*dxinv[0] + 1e-6);
    jstrt = (int)((objectCoord[iobj][2] - domain[2])*dxinv[1] + 1e-6);
    jend  = (int)((objectCoord[iobj][3] - domain[2])*dxinv[1] + 1e-6);

    j = jstrt;
    for(i = istrt+1; i < iend; i++) {
      long nid = ID(i,j,k);
      for(int ispec = 0; ispec < nspec; ispec++) {
        chgDenst[nid][ispec] *= 2;
      }
    }

    j = jend;
    for(i = istrt+1; i < iend; i++) {
      long nid = ID(i,j,k);
      for(int ispec = 0; ispec < nspec; ispec++) {
        chgDenst[nid][ispec] *= 2;
      }
    }

    i = istrt;
    for(j = jstrt+1; j < jend; j++) {
      long nid = ID(i,j,k);
      for(int ispec = 0; ispec < nspec; ispec++) {
        chgDenst[nid][ispec] *= 2;
      }
    }

    i = iend;
    for(j = jstrt; j < jend; j++) {
      long nid = ID(i,j,k);
      for(int ispec = 0; ispec < nspec; ispec++) {
        chgDenst[nid][ispec] *= 2;
      }
    }

  }   // end for(int iobj = 0; iobj < nobjects; iobj++)

  #undef ID
}
*/

void part2meshChg(int nspec, long *nparts,
                  const int ngpt[3],
                  double *partWeight,
                  double **chgDenst,
                  const int nthrds)
{
    int nx = ngpt[0], ny = ngpt[1], nz = ngpt[2];
    int nxy = nx * ny, nxyz = nxy * nz;

    static int first(1);
    if (first)
    {
        double Rj, Rjj;

        chgDenst_thrd = new double **[nthrds];
        for (int ithrd = 0; ithrd < nthrds; ithrd++)
        {
            chgDenst_thrd[ithrd] = new double *[nxyz];
            for (int id = 0; id < nxyz; id++)
                chgDenst_thrd[ithrd][id] = new double[nspec + 1];
        }

        volinv = new double[nxy];
        for (int j = 1; j + 1 < ny; ++j)
        {
            Rj = j * dx[1];

            for (int i = 0; i < nx; ++i)
            {
                double a = 0.0, dx_half = 0.5 * dx[0];

                for (int jj = -1; jj <= 1; jj += 2)
                {
                    if (j + jj < 0 || j + jj >= ny)
                        continue;

                    for (int ii = -1; ii <= 1; ii += 2)
                    {
                        if (i + ii < 0 || i + ii >= nx)
                            continue;

                        Rjj = 0.5 * (Rj + (j + jj) * dx[1]);
                        Rjj = Rj > Rjj ? (SQR(Rj) - SQR(Rjj)) : (SQR(Rjj) - SQR(Rj));
                        a += dx_half * Rjj;
                    }
                }

                a *= PI;

                volinv[j * nx + i] = 1.0 / a;
            }
        }

        // lower bound: j = 0, Rj = 0.0
        Rj = 0.0;
        Rjj = SQR(dx[1]);
        for (int i = 0; i < nx; ++i)
        {
            double a = 0.0, dx_half = 0.5 * dx[0];

            for (int ii = -1; ii <= 1; ii += 2)
            {
                if (i + ii < 0 || i + ii >= nx)
                    continue;
                a += dx_half * (Rjj - Rj);
            }

            a *= 0.5 * PI;

            volinv[0 * nx + i] = 1.0 / a;
        }

        // upper bound: j = ny-1
        Rj = SQR((ny - 1) * dx[1]);
        Rjj = SQR((ny - 2) * dx[1]);
        for (int i = 0; i < nx; ++i)
        {
            double a = 0.0, dx_half = 0.5 * dx[0];

            for (int ii = -1; ii <= 1; ii += 2)
            {
                if (i + ii < 0 || i + ii >= nx)
                    continue;
                a += dx_half * (Rj - Rjj);
            }

            a *= 0.5 * PI;

            volinv[(ny - 1) * nx + i] = 1.0 / a;
        }

        first = 0;
    }

    for (int ithrd = 0; ithrd < nthrds; ithrd++)
        for (int id = 0; id < nxyz; id++)
        {
            for (int ispec = 0; ispec <= nspec; ispec++)
                chgDenst_thrd[ithrd][id][ispec] = 0.0;
            //      memset(chgDenst_thrd[ithrd][id], 0, (nspec+1)*sizeof(double));
        }

    for (int id = 0; id < nxyz; id++)
    {
        for (int ispec = 0; ispec <= nspec; ispec++)
            chgDenst[id][ispec] = 0.0;
        //    memset(chgDenst[id], 0, (nspec+1)*sizeof(double));
    }

    long npartsTot = nparts[MAXNSPECTYPES];
    long nparts_per_thrd = npartsTot / nthrds;
    long thrd_start[MAXNTHRDS], thrd_end[MAXNTHRDS];

    for (int ithrd = 0; ithrd < nthrds; ithrd++)
    {
        thrd_start[ithrd] = nparts_per_thrd * ithrd;
        thrd_end[ithrd] = nparts_per_thrd * (ithrd + 1);
    }
    thrd_end[nthrds - 1] = npartsTot;
    for (long ipart = 0; ipart < npartsTot; ipart++)
        vacMap[ipart] = MAXNPARTS;

#pragma omp parallel default(none)              \
    shared(part, nparts, chgDenst_thrd, domain, \
               ngpt, nx, nxy, dx, dxinv, nspec, thrd_start, thrd_end, smallValue)
    {
#define igpt(i, j, k) ((k) * (nxy) + (j) * (nx) + (i))

        int tid = omp_get_thread_num();
        int ix, iy;
        double partX[3];
        double f[4];

        for (int ipart = thrd_start[tid]; ipart < thrd_end[tid]; ipart++)
        {
            int itype = (int)(part[ipart][ptype] + smallValue);
            partX[0] = part[ipart][xp];
            partX[1] = part[ipart][yp];
            partX[2] = part[ipart][zp];

            //      bilinear_weighting(partX, domain, ix, iy, f);

            area_weighting(partX, domain, ix, iy, f);

            chgDenst_thrd[tid][igpt(ix, iy, 0)][itype] += f[0];
            chgDenst_thrd[tid][igpt(ix + 1, iy, 0)][itype] += f[1];
            chgDenst_thrd[tid][igpt(ix, iy + 1, 0)][itype] += f[2];
            chgDenst_thrd[tid][igpt(ix + 1, iy + 1, 0)][itype] += f[3];

            for (int a = 0; a < 3; a++)
                part[ipart][xp + a] = partX[a];
            if (ngpt[2] == 1)
                part[ipart][zp] = 0.0;
        }
#undef igpt
    } // end #pragma omp parallel

#define igpt(i, j, k) ((k) * (nxy) + (j) * (nx) + (i))

    int nid = 0;
    for (int ithrd = 0; ithrd < nthrds; ithrd++)
        for (int k = 0; k < nz; k++)
            for (int j = 0; j < ny; j++)
                for (int i = 0; i < nx; i++)
                {
                    for (int itype = 0; itype < nspec; itype++)
                    {
                        nid = igpt(i, j, k);
                        chgDenst[nid][itype] += chgDenst_thrd[ithrd][nid][itype] * partWeight[itype] * spec[itype][specchg];
                        //            if(j == 0) {
                        //              printf("%e ", chgDenst[igpt(i,j,k)][itype]);
                        //              printf("%e ", chgDenst_thrd[1][igpt(i+1,j+1,k+1)][itype]);
                        //            }
                    }
                    //          if(j == 0) printf("\n");
                }

    // --------------------------------------------------------------------------
    //  axisymmetric
    // --------------------------------------------------------------------------
    //  /* boundaries at jmin and jmax */
    //  for(int i = 0; i < nx; i++) {
    //    int nid[2];
    //    nid[0] = igpt(i,0,0);
    //    nid[1] = igpt(i,ny-1,0);
    //
    //    meshChgBoundary(nid, nspec, &partBC[2], partWeight, chgDenst);
    //  }
    //
    //  /* boundaries at imin and imax */
    //  for(int j = 0; j < ny; j++) {
    //    int nid[2];
    //    nid[0] = igpt(0,j,0);
    //    nid[1] = igpt(nx-1,j,0);
    //
    //    meshChgBoundary(nid, nspec, &partBC[0], partWeight, chgDenst);
    //
    //  }

    //  meshChgObject(ngpt, nspec, partWeight, chgDenst);

    /* compute total charge density by summing over all species */
    for (int j = 0; j < ny; j++)
    {
        for (int i = 0; i < nx; i++)
        {
            int nid = igpt(i, j, 0);

            for (int ispec = 0; ispec < nspec; ispec++)
            {
                chgDenst[nid][ispec] *= volinv[nid];
                chgDenst[nid][nspec] += chgDenst[nid][ispec];
            }
        }
    }

#undef igpt
}
