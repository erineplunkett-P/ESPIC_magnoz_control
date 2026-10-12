#include "es_PIC.h"
#include "edf.h"

double E0[3] = {0., 0., 0.};
double B0[3] = {0., 0., 0.};
bool EextOn = false; // flag whether external E exists
bool BextOn = true;  // flag whether external B exists
int numBextCoil;
std::vector<std::vector<double>> BextCoil_Location;
std::vector<double> BextCoil_ThroatStrength;
double dx[3] = {1., 1., 1.};
double dxinv[3];
double domain[6] = {0.0, 1.0, 0.0, 1.0, 0.0, 1.0};
int fieldBC[6] = {0, 0, 0, 0, 0, 0};
int partBC[6] = {0, 0, 0, 0, 0, 0};
int partBC_new[6] = {0, 0, 0, 0, 0, 0};
double fbcValue[6];

int vthrd[3] = {1, 1, 1};

int MAXNPARTS = 10000;     // estimated max # of particles, read from input
int MAXNPARTSINCR = 10000; // increment of max # of particles per re-allocate part

double (*part)[MAXNPARTPROPS];
double spec[MAXNSPECTYPES][MAXNSPECPROPS];
double partWeight[MAXNSPECTYPES];

/* output data of particle information*/
void dump_particles(const int istep, long *nparts,
                    const vector<vector<double>> &dump_part_box,
                    const Vec &pot)
{
    const char odir[] = "dumped_particles";
    if (1)
    { // make a directory to put dumped particles if not exist
        struct stat st;
        int ierr = stat(odir, &st);

        if (ierr == 0)
        { // "field_data" exists
            if (!S_ISDIR(st.st_mode))
            {
                ostringstream oss;
                oss << "\"" << odir << "\" exists, but is not a directory";
                esPIC_exit(oss.str());
            }
        }
        else
        {
            mkdir(odir, S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
        }
    }

    long npartsTot = nparts[MAXNSPECTYPES];
    int nbox = (int)(dump_part_box.size());

    // # of nodes in the x-direction
    int nx = (int)((domain[1] - domain[0]) * dxinv[0] + 0.0000001) + 1;
    int ix, iy;

    for (int ibox = 0; ibox < nbox; ++ibox)
    {
        double xlo[2] = {dump_part_box[ibox][0], dump_part_box[ibox][2]};
        double xhi[2] = {dump_part_box[ibox][1], dump_part_box[ibox][3]};
        int ixlo[2] = {(int)(xlo[0] + 0.5), (int)(xlo[1] + 0.5)};
        int ixhi[2] = {(int)(xhi[0] + 0.5), (int)(xhi[1] + 0.5)};
        vector<int> specid;

        for (auto it = dump_part_box[ibox].cbegin() + 6; it != dump_part_box[ibox].cend(); ++it)
            specid.push_back((int)(*it + 0.00001));

        int nspec = (int)(specid.size());
        for (int ispec = 0; ispec < nspec; ++ispec)
        {
            ostringstream oss;
            oss << "./" << odir << "/particle_spec" << specid[ispec]
                << "_x" << ixlo[0] << "_" << ixhi[0]
                << "_y" << ixlo[1] << "_" << ixhi[1]
                << "_step" << istep << ".dat";

            ofstream outfile(oss.str().c_str());

            outfile << "# x0 y0 u0 v0 Ep0 x y u v Ep t_in\n";

            int itype;
            double x[2], y[2], u[2], v[2], Ep[2], tin;

            for (long ipart = 0; ipart < npartsTot; ipart++)
            {
                itype = (int)(part[ipart][ptype] + 0.000001);
                if (itype != specid[ispec])
                    continue;

                x[1] = part[ipart][xp];
                y[1] = part[ipart][yp];

                if (x[1] < xlo[0] || x[1] > xhi[0] || y[1] < xlo[1] || y[1] > xhi[1])
                    continue;

                // initial info of particle
                tin = part[ipart][pid];
                x[0] = part[ipart][xp0];
                y[0] = part[ipart][yp0];
                u[0] = part[ipart][up0];
                v[0] = part[ipart][vp0];
                Ep[0] = part[ipart][Epot];

                // current info of particle
                u[1] = part[ipart][up];
                v[1] = part[ipart][vp];

                ix = (int)((x[1] - domain[0]) * dxinv[0] + 0.5);
                iy = (int)((y[1] - domain[2]) * dxinv[1] + 0.5);
                Ep[1] = spec[itype][specchg] * pot(iy * nx + ix);

                for (int i = 0; i < 2; ++i)
                    outfile << x[i] << " " << y[i] << " "
                            << u[i] << " " << v[i] << " " << Ep[i] << " ";
                outfile << tin << "\n";
            } // for (long ipart = 0; ipart < npartsTot; ipart++)

            outfile.close();
        }
    }
}

void write_restart(const int step_now, const long npall, const string &ofname)
{
    FILE *ofp;
    ofp = fopen(ofname.c_str(), "wb");

    fwrite(&step_now, sizeof(int), 1, ofp);
    fwrite(&npall, sizeof(long), 1, ofp);
    fwrite(part, sizeof(double), npall * MAXNPARTPROPS, ofp);
    fclose(ofp);

    //  char buffer[1024];
    //  ofstream outfile(ofname.c_str(), std::ios::binary);
    //  memcpy(buffer, &step_now, sizeof(int));
    //  outfile.write(buffer, sizeof(int));
    //  memcpy(buffer, &npall, sizeof(long));
    //  outfile.write(buffer, sizeof(long));
    //
    //  size_t ns = MAXNPARTPROPS*sizeof(double);
    //  for (long ip = 0; ip < npall; ++ip) {
    //    memcpy(buffer, part[ip], ns);
    //    outfile.write(buffer, ns);
    //  }
    //  outfile.close();
}

void read_restart(int &step_now, int nspec, long nparts[MAXNSPECTYPES + 1])
{
    string ifname;
    struct stat st0, st1;

    int ierr0 = stat("restart.0", &st0);
    int ierr1 = stat("restart.1", &st1);

    if (ierr0 && ierr1)
    {
        cout << ("--> Warning: flag of reading restart file turned on, but no restart file to read. <--\n");
        cerr << ("--> Warning: flag of reading restart file turned on, but no restart file to read. <--\n");
        return;
    }

    if (ierr0 == 0 && ierr1 == 0)
    {
        if (difftime(st0.st_mtime, st1.st_mtime) > 0.0)
            ifname = "restart.0";
        else
            ifname = "restart.1";
    }
    else
    {
        if (ierr0 == 0)
            ifname = "restart.0";
        else if (ierr1 == 0)
            ifname = "restart.1";
    }

    cout << "Read restart file from [" << ifname << "]\n";

    FILE *ifp;
    ifp = fopen(ifname.c_str(), "rb");
    if (ifp == NULL)
    {
        ostringstream oss;
        oss << "Failed to open restart file [" << ifname << "]";
        esPIC_exit(oss.str());
    }

    long size_read;
    size_read = fread(&step_now, sizeof(int), 1, ifp);
    size_read = fread(&(nparts[MAXNSPECTYPES]), sizeof(long), 1, ifp);

    while (nparts[MAXNSPECTYPES] > MAXNPARTS)
        adjustPartMem();
    size_read = fread(&(part[0][0]), sizeof(double), nparts[MAXNSPECTYPES] * MAXNPARTPROPS, ifp);

    if (size_read != nparts[MAXNSPECTYPES] * MAXNPARTPROPS)
    {
        ostringstream oss;
        oss << "Failed to read particle info from [" << ifname << "]";
        esPIC_exit(oss.str());
    }

    for (long ip = 0; ip < nparts[MAXNSPECTYPES]; ++ip)
        ++nparts[(long)part[ip][ptype]];

    printf("istep = %d", step_now);
    for (int ispec = 0; ispec < nspec; ispec++)
    {
        printf(", np[%d]/np_col[%d] = %ld/%d", ispec, ispec, nparts[ispec], 0);
    }
    printf(", np[all] = %ld", nparts[MAXNSPECTYPES]);
    printf("\n");

    fclose(ifp);
}

/* read simulation control parameters */
bool initCtrl(double &dt, int &endstep,
              int &ivevalfield, int &ivdumpPart,
              int &ivdump, int &flag_read_restart, int &ivrestart,
              bool &dumpField, bool &dumpPart)
{
    int hybrid = 0;
    int seed;

    /* read control input parameters of simulation */
    ifstream infile("simulationCtrl.in");

    if (!infile)
        esPIC_exit("Input file [simulationCtrl.in] missing");
    else
    {
        cout << "Reading [simulationCtrl.in]..." << endl;
    }

    string line;
    getline(infile, line);
    istringstream iss(line);
    iss >> MAXNPARTS;
    cout << "MAXNPARTS = " << MAXNPARTS << endl;

    getline(infile, line);
    iss.str(line);
    iss >> MAXNPARTSINCR;
    cout << "MAXNPARTSINCR = " << MAXNPARTSINCR << endl;

    getline(infile, line);
    iss.str(line);
    iss >> dt;

    getline(infile, line);
    iss.str(line);
    iss >> endstep;
    cout << "dt = " << dt << ", endstep = " << endstep << endl;

    getline(infile, line);
    iss.str(line);
    iss >> hybrid;
    if (hybrid == 1)
    {
        cout << "Conducting hybrid PIC simulation." << endl;
    }
    else
    {
        cout << "Conducting full PIC simulation." << endl;
    }

    getline(infile, line);
    iss.str(line);
    iss >> seed;
    if (seed < 0)
        srand((int)time(NULL));
    else if (seed > 0)
        srand(seed);

    getline(infile, line);
    iss.str(line);
    iss >> ivevalfield;
    cout << "ivevalfield = " << ivevalfield << endl;

    getline(infile, line);
    iss.str(line);
    iss >> ivdump;
    cout << "ivdump = " << ivdump << endl;

    getline(infile, line);
    iss.str(line);
    iss >> ivdumpPart;
    cout << "ivdumpPart = " << ivdumpPart << endl;

    getline(infile, line);
    iss.str(line);
    iss >> flag_read_restart;
    flag_read_restart = (flag_read_restart == 1) ? 1 : 0;

    getline(infile, line);
    iss.str(line);
    iss >> ivrestart;
    cout << "ivrestart = " << ivrestart << endl;

    if (ivdump >= ivevalfield && ivevalfield != 0)
    {
        dumpField = true;

        if ((ivdump % ivevalfield) != 0)
            esPIC_exit("ivdump (mod) ivevalfield != 0");
    }

    if (ivdumpPart > 0)
        dumpPart = true;

    getline(infile, line);
    iss.str(line);
    iss >> vthrd[0] >> vthrd[1] >> vthrd[2];
    cout << "virtual openmp threads vthrd[3] = {" << vthrd[0] << ", " << vthrd[1] << ", " << vthrd[2] << "}" << endl;

    getline(infile, line);
    iss.str(line);
    iss >> EextOn;
    cout << "EextOn = " << EextOn << endl;

    getline(infile, line);
    iss.str(line);
    iss >> BextOn;
    cout << "BextOn = " << BextOn << endl;

    getline(infile, line);
    iss.str(line);
    iss >> numBextCoil;
    cout << "numBextCoil = " << numBextCoil << endl;

    BextCoil_Location.resize(numBextCoil);
    BextCoil_ThroatStrength.resize(numBextCoil);
    for (int i = 0; i < numBextCoil; i++)
    {
        getline(infile, line);
        iss.str(line);
        double xa, ya, za;
        iss >> xa >> ya >> za >> BextCoil_ThroatStrength[i];
        BextCoil_Location[i].push_back(xa);
        BextCoil_Location[i].push_back(ya);
        BextCoil_Location[i].push_back(za);
    }

    infile.close();

    return (hybrid == 1);
}

void init_part_dump(vector<vector<double>> &part_dump_box,
                    long *nparts, const int istep, const int ivdumpPart,
                    const Vec &pot)
{

    /* read control input parameters of simulation */
    ifstream infile("particleOutputCtrl.in");

    if (!infile)
        esPIC_exit("Input file [particleOutputCtrl.in] missing");
    else
    {
        cout << "Reading [particleOutputCtrl.in]..." << endl;
    }

    string line;
    while (getline(infile, line))
    {
        auto sit = line.begin();
        while (sit != line.end())
        {
            if (*sit == '!')
                break;
            ++sit;
        }

        string content(line.begin(), sit);
        istringstream iss(content);

        vector<double> helper;
        double item;
        while (iss >> item)
        {
            helper.push_back(item);
        }
        part_dump_box.push_back(helper);
    }

    infile.close();

    if (istep > 0 && (istep % ivdumpPart) == 0)
        dump_particles(istep, nparts, part_dump_box, pot);

    return;
}

double RTOL = 1e-6;   // tolerance for iterative poisson solver
int MAXIT = 10000;    // max # of iterations for poisson solver
double RTOLNT = 1e-6; // tol for Newton's method (non-linear poisson)
int MAXITNT = 20;     // max # of iterations for Newton's method

void read_field_solver(int &solverType)
{
    ifstream infile("field_solver.in");
    string line;

    if (!infile)
        esPIC_exit("Input file \"field_solver.in\" does NOT exist");

    getline(infile, line);
    istringstream isstrm(line);
    isstrm >> solverType;

    getline(infile, line);
    isstrm.str(line);
    isstrm >> RTOL;

    getline(infile, line);
    isstrm.str(line);
    isstrm >> MAXIT;

    getline(infile, line);
    isstrm.str(line);
    isstrm >> RTOLNT;

    getline(infile, line);
    isstrm.str(line);
    isstrm >> MAXITNT;

    infile.close();
    infile.clear();
}

PtrPoissonSolver init_field_solver(const int ngpt[3],
                                   const bool HYBRID,
                                   const ObjectList &objlist,
                                   double **chgDenst,
                                   SpMat &A,
                                   Vec &b,
                                   Vec &p,
                                   MapF &fixedBV)
{
    /* read control parameter for field solver */
    int solver = 1;
    read_field_solver(solver);

    //  long n = ngpt[0]*ngpt[1]*ngpt[2];

    for (long r = 0; r < ngpt[0] * ngpt[1] * ngpt[2]; ++r)
        p[r] = 0.0;

    if (HYBRID)
    {
        /* Set up the matrix */
        // stiffness matrix same as linear poisson
        setup_matrix_linearFDS(ngpt, dx, A, objlist, fixedBV);

        /* Set up right hand side */
        // b used to store charge density
        setup_rhs_linearFDS(ngpt, b);
    }
    else
    {
        /* Set up the matrix */
        setup_matrix_linearFDS(ngpt, dx, A, objlist, fixedBV);
        /* Set up right hand side */
        setup_rhs_linearFDS(ngpt, b);
    }

    if (HYBRID)
    {
        switch (solver)
        {
        case CG_SOLVER:
            cout << "solver selection: non-linear CG (max iterations = " << MAXIT << ")" << endl;
            return bicgstab_nonlinear;
        case LU_SOLVER:
            cout << "solver selection: non-linear Direct Method with LU decomposition." << endl;
            return lu_nonlinear;
        default:
            esPIC_exit("solver type has not been specified");
        }
    }
    else
    {
        switch (solver)
        {
        case CG_SOLVER:
            cout << "solver selection: linear CG (max iterations = " << MAXIT << ")" << endl;
            return bicgstab_solver;
        case LU_SOLVER:
            cout << "solver selection: linear Direct Method with LU decomposition." << endl;
            return lu_solver;
        default:
            esPIC_exit("solver type has not been specified");
        }
    }

    return 0;
}

int genMesh(int ngpt[3])
{
    ifstream infile("meshGeom.in");

    if (!infile)
        return 0; // return generate mesh flag (0 - no mesh generated)
    else
    {
        cout << "Generating mesh according to [meshGeom.in]..." << endl;
    }

    string line;
    getline(infile, line);
    istringstream iss(line);
    iss >> domain[0] >> domain[1] >> domain[2] >> domain[3] >> domain[4] >> domain[5];

    getline(infile, line);
    iss.str(line);
    iss >> dx[0] >> dx[1] >> dx[2];

    getline(infile, line);
    iss.str(line);
    iss >> fieldBC[0] >> fieldBC[1] >> fieldBC[2] >> fieldBC[3] >> fieldBC[4] >> fieldBC[5];

    getline(infile, line);
    iss.str(line);
    iss >> partBC[0] >> partBC[1] >> partBC[2] >> partBC[3] >> partBC[4] >> partBC[5];

    getline(infile, line);
    iss.str(line);
    iss >> partBC_new[0] >> partBC_new[1] >> partBC_new[2] >> partBC_new[3] >> partBC_new[4] >> partBC_new[5];

    getline(infile, line);
    iss.str(line);
    iss >> fbcValue[0] >> fbcValue[1] >> fbcValue[2] >> fbcValue[3] >> fbcValue[4] >> fbcValue[5];

    infile.close();
    infile.clear();

    for (int a = 0; a < 3; a++)
        dxinv[a] = 1.0 / dx[a];
    for (int a = 0; a < 3; a++)
        ngpt[a] = (int)((domain[2 * a + 1] - domain[2 * a]) * dxinv[a] + 1e-5) + 1;

    return ngpt[0] * ngpt[1] * ngpt[2];
}

void initField(const int ngpt[3],
               double (*&Eext)[3], double (*&Bext)[3],
               double (*E)[3], double (*B)[3])
{
    int nx = ngpt[0], ny = ngpt[1], nz = ngpt[2];
    long nxyz = nx * ny * nz;

    memset(E, 0, 3 * nxyz * sizeof(double));
    memset(B, 0, 3 * nxyz * sizeof(double));

#define igpt(i, j, k) (((k) * (ny) + (j)) * (nx) + (i))

    if (EextOn)
    {
        double A = 10, a = 0.5 * PI;

        Eext = new double[nxyz][3];

        for (int k = 0; k < nz; k++)
            for (int j = 0; j < ny; j++)
                for (int i = 0; i < nx; i++)
                {
                    long nid = igpt(i, j, k);

                    Eext[nid][0] = A * sin(a * i * dx[0]);
                    Eext[nid][1] = Eext[nid][2] = 0.0;
                    //          for(int a = 0; a < 3; a++)
                    //            Eext[nid][a] = E0[a];
                    for (int a = 0; a < 3; a++)
                        E[nid][a] += Eext[nid][a];
                }
    }

    double *Bx = new double[nxyz];
    double *By = new double[nxyz];
    double *XX = new double[nxyz];
    double *YY = new double[nxyz];

    std::fill_n(Bx, nxyz, 0.0);
    std::fill_n(By, nxyz, 0.0);
    std::fill_n(XX, nxyz, 0.0);
    std::fill_n(YY, nxyz, 0.0);

    if (BextOn)
    {
        Bext = new double[nxyz][3];

        // for (int k = 0; k < nz; k++)
        //     for (int j = 0; j < ny; j++)
        //         for (int i = 0; i < nx; i++)
        //         {
        //             long nid = igpt(i, j, k);

        //             for (int a = 0; a < 3; a++)
        //                 Bext[nid][a] = B0[a];

        //             for (int a = 0; a < 3; a++)
        //                 E[nid][a] += Eext[nid][a];
        //         }
        memset(Bext, 0, 3 * nxyz * sizeof(double));

        // for (int k = 0; k < nz; k++)
        //     for (int j = 0; j < ny; j++)
        //         for (int i = 0; i < nx; i++)
        //         {
        //             long nid = igpt(i, j, k);

        //             for (int a = 0; a < 3; a++)
        //                 Bext[nid][a] = B0[a];

        //             for (int a = 0; a < 3; a++)
        //                 E[nid][a] += Eext[nid][a];
        //         }
        for (int icoil = 0; icoil < numBextCoil; icoil++)
        {
            double Xa = BextCoil_Location[icoil][0];
            double Ya = BextCoil_Location[icoil][1];
            double b_throat = BextCoil_ThroatStrength[icoil];
            // For safety, clear all of the temporary arrays
            std::fill_n(Bx, nxyz, 0.0);
            std::fill_n(By, nxyz, 0.0);
            std::fill_n(XX, nxyz, 0.0);
            std::fill_n(YY, nxyz, 0.0);
            double kk_t = sqrt(4.0 * Ya * 0.0 / (pow((0.0 - Xa), 2.0) + pow((0.0 + Ya), 2.0)));
            double tr2_t = pow(0.0 + Ya, 2.0) + pow(0.0 - Xa, 2.0);
            double tr4_t = pow(0.0 - Ya, 2.0) + pow(0.0 - Xa, 2.0);
            double tr_z_t = -pow(0.0 - Xa, 2.0) - 0.0 * 0.0 + Ya * Ya;
            // double tre1_t = M_PI / 2.0 * (1.0 + 0.25 * pow(kk_t, 2.0) + 0.140625 * pow(kk_t, 4.0) + 0.09765625 * pow(kk_t, 6.0) + 0.07476806640625 * pow(kk_t, 8.0));
            // double tre2_t = M_PI / 2.0 * (1.0 - 0.25 * pow(kk_t, 2.0) - 0.140625 * pow(kk_t, 4.0) / 3.0 - 0.09765625 * pow(kk_t, 6.0) / 5.0 - 0.07476806640625 * pow(kk_t, 8.0) / 7.0);
            //  double tre1_t = std::comp_ellint_1(kk_t);
            //  double tre2_t = std::comp_ellint_2(kk_t);
            double tre1_t = ellint::comp_ellint_1(kk_t);
            double tre2_t = ellint::comp_ellint_2(kk_t);
            double coeff = b_throat / (pow(tr2_t, -0.5) * (tr_z_t / tr4_t * tre2_t + tre1_t));
            for (int k = 0; k < nz; k++)
            {
                for (int j = 1; j < ny; j++)
                {
                    for (int i = 0; i < nx; i++)
                    {
                        long index = igpt(i, j, k);
                        XX[index] = i * dx[0];
                        YY[index] = j * dx[1];
                        double kk = sqrt(4.0 * Ya * YY[index] / (pow((XX[index] - Xa), 2.0) + pow((YY[index] + Ya), 2.0)));
                        double tr2 = pow(YY[index] + Ya, 2.0) + pow(XX[index] - Xa, 2.0);
                        double tr3 = pow(tr2, -0.5) * (XX[index] - Xa) / YY[index];
                        double tr4 = pow(YY[index] - Ya, 2.0) + pow(XX[index] - Xa, 2.0);
                        double tr_r = pow(XX[index] - Xa, 2.0) + YY[index] * YY[index] + Ya * Ya;
                        double tr_z = -pow(XX[index] - Xa, 2.0) - YY[index] * YY[index] + Ya * Ya;
                        // double tre1 = M_PI / 2.0 * (1.0 + 0.25 * pow(kk, 2.0) + 0.140625 * pow(kk, 4.0) + 0.09765625 * pow(kk, 6.0) + 0.07476806640625 * pow(kk, 8.0));
                        // double tre2 = M_PI / 2.0 * (1.0 - 0.25 * pow(kk, 2.0) - 0.140625 * pow(kk, 4.0) / 3.0 - 0.09765625 * pow(kk, 6.0) / 5.0 - 0.07476806640625 * pow(kk, 8.0) / 7.0);
                        // double tre1 = std::comp_ellint_1(kk);
                        // double tre2 = std::comp_ellint_2(kk);
                        double tre1 = ellint::comp_ellint_1(kk);
                        double tre2 = ellint::comp_ellint_2(kk);
                        By[index] = coeff * tr3 * (tr_r / tr4 * tre2 - tre1);
                        Bx[index] = coeff * pow(tr2, -0.5) * (tr_z / tr4 * tre2 + tre1);
                    }
                }
                // here j = 0
                int j = 0;
                for (int i = 0; i < nx; i++)
                {
                    long index = igpt(i, j, k);
                    XX[index] = i * dx[0];
                    YY[index] = j * dx[1];
                    double kk = sqrt(4.0 * Ya * YY[index] / (pow((XX[index] - Xa), 2.0) + pow((YY[index] + Ya), 2.0)));
                    double tr2 = pow(YY[index] + Ya, 2.0) + pow(XX[index] - Xa, 2.0);
                    double tr4 = pow(YY[index] - Ya, 2.0) + pow(XX[index] - Xa, 2.0);
                    double tr_z = -pow(XX[index] - Xa, 2.0) - YY[index] * YY[index] + Ya * Ya;
                    // double tre1 = M_PI / 2.0 * (1.0 + 0.25 * pow(kk, 2.0) + 0.140625 * pow(kk, 4.0) + 0.09765625 * pow(kk, 6.0) + 0.07476806640625 * pow(kk, 8.0));
                    // double tre2 = M_PI / 2.0 * (1.0 - 0.25 * pow(kk, 2.0) - 0.140625 * pow(kk, 4.0) / 3.0 - 0.09765625 * pow(kk, 6.0) / 5.0 - 0.07476806640625 * pow(kk, 8.0) / 7.0);
                    // double tre1 = std::comp_ellint_1(kk);
                    // double tre2 = std::comp_ellint_2(kk);
                    double tre1 = ellint::comp_ellint_1(kk);
                    double tre2 = ellint::comp_ellint_2(kk);
                    By[index] = 0.0;
                    Bx[index] = coeff * pow(tr2, -0.5) * (tr_z / tr4 * tre2 + tre1);
                }
            }

            for (int k = 0; k < nz; k++)
            {
                for (int j = 0; j < ny; j++)
                {
                    for (int i = 0; i < nx; i++)
                    {
                        long index = igpt(i, j, k);
                        Bext[index][0] += Bx[index];
                        Bext[index][1] += By[index];
                        Bext[index][2] = 0.0;
                    }
                }
            }
        }
        for (int k = 0; k < nz; k++)
        {
            for (int j = 0; j < ny; j++)
            {
                for (int i = 0; i < nx; i++)
                {
                    long index = igpt(i, j, k);
                    B[index][0] = Bext[index][0];
                    B[index][1] = Bext[index][1];
                    B[index][2] = Bext[index][2];
                }
            }
        }
        std::ofstream outfile;
        std::string outfilepath;

        outfilepath = "./Bz.dat";
        outfile.open(outfilepath, std::ofstream::out);
        for (int k = 0; k < nz; k++)
        {
            for (int j = 0; j < ny; j++)
            {
                for (int i = 0; i < nx; i++)
                {
                    long index = igpt(i, j, k);
                    outfile << B[index][0] << " ";
                }
                outfile << std::endl;
            }
            outfile << std::endl;
        }
        outfile.close();
        outfile.clear();

        outfilepath = "./Br.dat";
        outfile.open(outfilepath, std::ofstream::out);
        for (int k = 0; k < nz; k++)
        {
            for (int j = 0; j < ny; j++)
            {
                for (int i = 0; i < nx; i++)
                {
                    long index = igpt(i, j, k);
                    outfile << B[index][1] << " ";
                }
                outfile << std::endl;
            }
            outfile << std::endl;
        }
        outfile.close();
        outfile.clear();
    }

    delete[] Bx;
    delete[] By;
    delete[] XX;
    delete[] YY;

#undef igpt
}

int initSpec()
{
    int nspec = 0;
    /* read in [spec.in] to properties of different species */
    /* hard coded here - electron: type 0; ion: type 1 */

    ifstream infile("specProp.in");
    if (!infile)
        esPIC_exit("Input file [specProp.in] missing");
    else
    {
        cout << "Reading [specProp.in]..." << endl;
    }

    string line;
    getline(infile, line);
    while (getline(infile, line))
    {
        istringstream iss(line);
        iss >> spec[nspec][specmass] >> spec[nspec][specchg] >> spec[nspec][specFixed];
        spec[nspec][specQM] = spec[nspec][specchg] / spec[nspec][specmass];
        ++nspec;
    }

    return nspec;
}

void evalVel(const int nspec,
             const long *nparts,
             double (*vel)[3],
             double (*vth)[3])
{
    long npartsTot = nparts[MAXNSPECTYPES];

    for (int ispec = 0; ispec < nspec; ispec++)
        for (int a = 0; a < 3; a++)
        {
            vel[ispec][a] = 0.0;
            vth[ispec][a] = 0.0;
        }

    for (long ipart = 0; ipart < npartsTot; ipart++)
    {
        int ispec = (int)(part[ipart][ptype] + 1e-6);
        double partV[3];
        for (int a = 0; a < 3; a++)
            partV[a] = part[ipart][up + a];
        for (int a = 0; a < 3; a++)
        {
            vel[ispec][a] += partV[a];
            vth[ispec][a] += SQR(partV[a]);
        }
    }

    for (int ispec = 0; ispec < nspec; ispec++)
    {
        for (int a = 0; a < 3; a++)
        {
            vel[ispec][a] /= nparts[ispec];
            vth[ispec][a] /= nparts[ispec];

            vth[ispec][a] = sqrt(vth[ispec][a] - SQR(vel[ispec][a]));
        }
    }
}

int main(int argc, char *argv[])
{
    double tstart, tend, elapsed = 0.0;
    double t_inj_start, t_inj_end, inj_elapsed = 0.0;
    double t_pmove_start, t_pmove_end, pmove_elapsed = 0.0;
    double t_depos_start, t_depos_end, depos_elapsed = 0.0;
    double t_fieldSol_start = 0.0, t_fieldSol_end = 0.0, fieldSol_elapsed = 0.0;
    double dt;
    int endstep;
    int ivdump, ivrestart, ivevalfield, ivdumpPart;
    int rstrt_id = 1, flag_read_restart = 0;
    bool dumpField = false, dumpPart = false;

    // flag: 0 - full PIC; 1 - hybrid simulation
    bool HYBRID = initCtrl(dt, endstep,
                           ivevalfield, ivdumpPart,
                           ivdump, flag_read_restart, ivrestart,
                           dumpField, dumpPart);
    double half_dt = 0.5 * dt;

    /* # of omp threads per process */
    int nthrds = vthrd[0] * vthrd[1] * vthrd[2];
    if (nthrds > MAXNTHRDS)
    {
        ostringstream oss;
        oss << "nthrds (" << nthrds << ") requested exceeds the maxmium (" << MAXNTHRDS << ")" << endl;
        esPIC_exit(oss.str());
    }

    /* set # of omp threads launched by each process */
    omp_set_num_threads(nthrds);

    /* Get information on all species */
    int nspec = initSpec(); // # of species considered
    cout << "nspec = " << nspec << endl;

    /* Initialize mesh and field */
    int ngpt[3];
    long nxyz = genMesh(ngpt);
    double (*E)[3] = new double[nxyz][3];
    double (*B)[3] = new double[nxyz][3];
    double (*Eext)[3] = NULL, (*Bext)[3] = NULL;

    /* Initialize external E and B fields */
    initField(ngpt, Eext, Bext, E, B);
    cout << "Memory allocated for (*E)[3] = " << (double)(sizeof(double) * nxyz * 3) / MB << " MB" << endl;
    cout << "Memory allocated for (*B)[3] = " << (double)(sizeof(double) * nxyz * 3) / MB << " MB" << endl;

    /* Initialize array to store charge density at every node */
    double **chgDenst = new double *[nxyz];
    for (long nid = 0; nid < nxyz; nid++)
        chgDenst[nid] = new double[nspec + 1]; // total charge density stored in chg[id][nspec]
    cout << "Memory allocated for **chgDenst = " << (double)(sizeof(double) * (nspec + 1) * nxyz) / MB << " MB" << endl;

    int istep = 0;
    long nparts[MAXNSPECTYPES + 1];

    // initialize part[MAXNPARTS][MAXNPARTPROPS]
    initPart(flag_read_restart, nparts, partWeight);
    if (flag_read_restart)
        read_restart(istep, nspec, nparts); // read restart file

    SampleField fdsample("sampleCtrl.in", istep, nspec, ngpt);
    VelDist vdf("vdfCtrl.in", istep, nspec, ngpt);
    EnergyDist edf("edfCtrl.in", istep, nspec, ngpt);

    part2meshChg(nspec, nparts, ngpt, partWeight, chgDenst, nthrds);

    ObjectList objlist("objectsIndomain.in");
    cout << "# of objects inside domain = " << objlist.size() << endl;

    const long nrows(nxyz), ncols(nxyz);
    SpMat Mat_A(nrows, ncols); // Compressed sparse column (CSC) format matrix
    Vec b(nrows);              // rhs vector
    Vec pot(nrows);            // solution vector

    /* Initialize the field solver */
    MapF fixedBV;                    // phi value for fixed phi points: map<nid, phi>
    PtrPoissonSolver poisson_solver; // ptr to solver function
    poisson_solver = init_field_solver(ngpt, HYBRID, objlist, chgDenst, Mat_A, b, pot, fixedBV);

    cout << "step = " << istep << ": total particles = " << nparts[MAXNSPECTYPES] << endl;

    int nFieldEvals = dumpField ? ivdump / ivevalfield : 0;
    vector<double> histPhi, histChgDenst, histVelField, new_histVelField, new_new_histVelField;
    vector<double> num_density, vel_azimuthal, momentum_azimuthal, stress_azimuthal, ExB_component, diamagnetic_component, inertial_component, stressdiv_component, J_theta_total_fluid, J_theta_PIC; //, J_theta_total_fluid_20, J_theta_total_fluid_21, J_theta_PIC_20, J_theta_PIC_21;
    histPhi.reserve(nFieldEvals * nxyz);
    histChgDenst.reserve(nFieldEvals * nxyz * (nspec + 1));
    histVelField.reserve(nFieldEvals * nxyz * nspec * 6);
    new_histVelField.reserve(nFieldEvals * nxyz * nspec * 5);
    new_new_histVelField.reserve(nFieldEvals * nxyz * nspec * 3);
    num_density.reserve(nFieldEvals * 1 * nspec * nxyz);
    vel_azimuthal.reserve(nFieldEvals * 1 * nspec * nxyz);
    ExB_component.reserve(nFieldEvals * 1 * nspec * nxyz);
    diamagnetic_component.reserve(nFieldEvals * 1 * nspec * nxyz);
    inertial_component.reserve(nFieldEvals * 1 * nspec * nxyz);
    stressdiv_component.reserve(nFieldEvals * 1 * nspec * nxyz);
    J_theta_total_fluid.reserve(nFieldEvals * 1 * nspec * nxyz);
    // J_theta_total_fluid_20.reserve(nFieldEvals * 1 * nspec * nxyz);
    // J_theta_total_fluid_21.reserve(nFieldEvals * 1 * nspec * nxyz);
    J_theta_PIC.reserve(nFieldEvals * 1 * nspec * nxyz);
    // J_theta_PIC_20.reserve(nFieldEvals * 1 * nspec * nxyz);
    // J_theta_PIC_21.reserve(nFieldEvals * 1 * nspec * nxyz);
    momentum_azimuthal.reserve(nFieldEvals * 3 * nspec * nxyz);
    stress_azimuthal.reserve(nFieldEvals * 3 * nspec * nxyz);
    cout << "histPhi.size() = " << histPhi.size() << ", histPhi.capacity() = " << histPhi.capacity()
         << ", memory required for " << nFieldEvals << " evaluations = "
         << (double)(sizeof(vector<double>) * histPhi.capacity()) / MB << " MB" << endl;
    cout << "histChgDenst.size() = " << histChgDenst.size() << ", histChgDenst.capacity() = " << histChgDenst.capacity()
         << ", memory required for " << nFieldEvals << " evaluations = "
         << (double)(sizeof(vector<double>) * histChgDenst.capacity()) / MB << " MB" << endl;
    cout << "histVelField.size() = " << histVelField.size() << ", histVelField.capacity() = " << histVelField.capacity()
         << ", memory required for " << nFieldEvals << " evaluations = "
         << (double)(sizeof(vector<double>) * histVelField.capacity()) / MB << " MB" << endl;
    cout << "new_histVelField.size() = " << new_histVelField.size() << ", new_histVelField.capacity() = " << new_histVelField.capacity()
         << ", memory required for " << nFieldEvals << " evaluations = "
         << (double)(sizeof(vector<double>) * new_histVelField.capacity()) / MB << " MB" << endl;
    cout << "new_new_histVelField.size() = " << new_new_histVelField.size() << ", new_new_histVelField.capacity() = " << new_new_histVelField.capacity()
         << ", memory required for " << nFieldEvals << " evaluations = "
         << (double)(sizeof(vector<double>) * new_new_histVelField.capacity()) / MB << " MB" << endl;
    cout << "num_density.size() = " << num_density.size() << ", num_density.capacity() = " << num_density.capacity()
         << ", memory required for " << nFieldEvals << " evaluations = "
         << (double)(sizeof(vector<double>) * num_density.capacity()) / MB << " MB" << endl;
    cout << "vel_azimuthal.size() = " << vel_azimuthal.size() << ", vel_azimuthal.capacity() = " << vel_azimuthal.capacity()
         << ", memory required for " << nFieldEvals << " evaluations = "
         << (double)(sizeof(vector<double>) * vel_azimuthal.capacity()) / MB << " MB" << endl;
    cout << "ExB_component.size() = " << ExB_component.size() << ", ExB_component.capacity() = " << ExB_component.capacity()
         << ", memory required for " << nFieldEvals << " evaluations = "
         << (double)(sizeof(vector<double>) * ExB_component.capacity()) / MB << " MB" << endl;
    cout << "diamagnetic_component.size() = " << diamagnetic_component.size() << ", diamagnetic_component.capacity() = " << diamagnetic_component.capacity()
         << ", memory required for " << nFieldEvals << " evaluations = "
         << (double)(sizeof(vector<double>) * diamagnetic_component.capacity()) / MB << " MB" << endl;
    cout << "inertial_component.size() = " << inertial_component.size() << ", inertial_component.capacity() = " << inertial_component.capacity()
         << ", memory required for " << nFieldEvals << " evaluations = "
         << (double)(sizeof(vector<double>) * inertial_component.capacity()) / MB << " MB" << endl;
    cout << "stressdiv_component.size() = " << stressdiv_component.size() << ", stressdiv_component.capacity() = " << stressdiv_component.capacity()
         << ", memory required for " << nFieldEvals << " evaluations = "
         << (double)(sizeof(vector<double>) * stressdiv_component.capacity()) / MB << " MB" << endl;
    cout << "J_theta_total_fluid.size() = " << J_theta_total_fluid.size() << ", J_theta_total_fluid.capacity() = " << J_theta_total_fluid.capacity()
         << ", memory required for " << nFieldEvals << " evaluations = "
         << (double)(sizeof(vector<double>) * J_theta_total_fluid.capacity()) / MB << " MB" << endl;
    // cout << "J_theta_total_fluid_20.size() = " << J_theta_total_fluid_20.size() << ", J_theta_total_fluid_20.capacity() = " << J_theta_total_fluid_20.capacity()
    //      << ", memory required for " << nFieldEvals << " evaluations = "
    //      << (double)(sizeof(vector<double>) * J_theta_total_fluid_20.capacity()) / MB << " MB" << endl;
    // cout << "J_theta_total_fluid_21.size() = " << J_theta_total_fluid_21.size() << ", J_theta_total_fluid_21.capacity() = " << J_theta_total_fluid_21.capacity()
    //      << ", memory required for " << nFieldEvals << " evaluations = "
    //      << (double)(sizeof(vector<double>) * J_theta_total_fluid_21.capacity()) / MB << " MB" << endl;
    cout << "J_theta_PIC.size() = " << J_theta_PIC.size() << ", J_theta_PIC.capacity() = " << J_theta_PIC.capacity()
         << ", memory required for " << nFieldEvals << " evaluations = "
         << (double)(sizeof(vector<double>) * J_theta_PIC.capacity()) / MB << " MB" << endl;
    // cout << "J_theta_PIC_20.size() = " << J_theta_PIC_20.size() << ", J_theta_PIC_20.capacity() = " << J_theta_PIC_20.capacity()
    //      << ", memory required for " << nFieldEvals << " evaluations = "
    //      << (double)(sizeof(vector<double>) * J_theta_PIC_20.capacity()) / MB << " MB" << endl;
    // cout << "J_theta_PIC_21.size() = " << J_theta_PIC_21.size() << ", J_theta_PIC_21.capacity() = " << J_theta_PIC_21.capacity()
    //      << ", memory required for " << nFieldEvals << " evaluations = "
    //      << (double)(sizeof(vector<double>) * J_theta_PIC_21.capacity()) / MB << " MB" << endl;
    cout << "momentum_azimuthal.size() = " << momentum_azimuthal.size() << ", momentum_azimuthal.capacity() = " << momentum_azimuthal.capacity()
         << ", memory required for " << nFieldEvals << " evaluations = "
         << (double)(sizeof(vector<double>) * momentum_azimuthal.capacity()) / MB << " MB" << endl;
    cout << "stress_azimuthal.size() = " << stress_azimuthal.size() << ", stress_azimuthal.capacity() = " << stress_azimuthal.capacity()
         << ", memory required for " << nFieldEvals << " evaluations = "
         << (double)(sizeof(vector<double>) * stress_azimuthal.capacity()) / MB << " MB" << endl;

    /* initialize potential field */
    if (nparts[MAXNSPECTYPES] > 0)
    {
#ifndef NOPOISSON
        if (HYBRID) // if conducting hybrid simulations, nonlinear field solver is chosen
            solve_field_nonlinear(ngpt, nspec, Mat_A, fixedBV, chgDenst, pot, b, E, Eext, poisson_solver);
        else
            solve_field_linear(ngpt, nspec, Mat_A, fixedBV, chgDenst, pot, b, E, Eext, poisson_solver);
#endif
        printf("\n");
    }

    vector<vector<double>> part_dump_box;
    if (dumpPart)
        init_part_dump(part_dump_box, nparts, istep, ivdumpPart, pot);

    /* Start PIC simulation loop till max time reached */
    GET_TIME(tstart);

    Beams beamlist("beams_in_domain.in");
    Sources sourcelist("sources_in_domain.in");

    while (istep < endstep)
    {
        ++istep;

        /* move particles by leap-frog scheme */
        GET_TIME(t_pmove_start);
        partMove(half_dt, ngpt, nparts, part, E, B, objlist, sourcelist, nthrds);
        GET_TIME(t_pmove_end);
        pmove_elapsed += (t_pmove_end - t_pmove_start);

        /* inject particles at in-/outflow boundaries */
        GET_TIME(t_inj_start);
        injectFlowInOut(nspec, nparts, ngpt, partWeight, dt, istep, pot);
        beamlist.inject_particles(dt, istep, nparts, pot);
        sourcelist.inject_particles(dt, istep, nparts, pot);
        GET_TIME(t_inj_end);
        inj_elapsed += (t_inj_end - t_inj_start);

        /* scatter - deposit particle charges to mesh */
        GET_TIME(t_depos_start);
        part2meshChg(nspec, nparts, ngpt, partWeight, chgDenst, nthrds);
        GET_TIME(t_depos_end);
        depos_elapsed += (t_depos_end - t_depos_start);

        /* solve linear Poisson equation for field */
        printf("istep = %d", istep);

        GET_TIME(t_fieldSol_start);
#ifndef NOPOISSON
        if (HYBRID) // if conducting hybrid simulations, nonlinear field solver is chosen
            solve_field_nonlinear(ngpt, nspec, Mat_A, fixedBV, chgDenst, pot, b, E, Eext, poisson_solver);
        else
            solve_field_linear(ngpt, nspec, Mat_A, fixedBV, chgDenst, pot, b, E, Eext, poisson_solver);
#endif
        if (HYBRID) // if conducting hybrid simulations, nonlinear field solver is chosen
            GET_TIME(t_fieldSol_end);
        fieldSol_elapsed += (t_fieldSol_end - t_fieldSol_start);

        for (int ispec = 0; ispec < nspec; ispec++)
        {
            printf(", np[%d]/np_col[%d] = %ld/%d", ispec, ispec, nparts[ispec], sourcelist.np_col_spec(ispec));
        }
        printf(", np[all] = %ld", nparts[MAXNSPECTYPES]);
        printf("\n");

        if (dumpField && !(istep % ivevalfield))
        {
            for (long nid = 0; nid < nxyz; nid++)
            {
                histPhi.push_back(pot[nid]);
                for (int ispec = 0; ispec <= nspec; ispec++)
                    histChgDenst.push_back(chgDenst[nid][ispec]);
            }

            evalHistVelField(nspec, nparts, ngpt, part, histVelField);
            evalNewHistVelField(nspec, nparts, ngpt, part, new_histVelField);
            evalNewNewHistVelField(nspec, nparts, ngpt, part, new_new_histVelField);
            evalAzimuthal(nspec, nparts, ngpt, part, partWeight, num_density, vel_azimuthal, momentum_azimuthal, stress_azimuthal,
                          ExB_component, diamagnetic_component, inertial_component,
                          stressdiv_component, J_theta_total_fluid, J_theta_PIC,
                          histVelField, E, B, 0, 1);
            // evalAzimuthal(nspec, nparts, ngpt, part, partWeight, num_density, vel_azimuthal, momentum_azimuthal, stress_azimuthal,
            //               ExB_component, diamagnetic_component, inertial_component,
            //               stressdiv_component, J_theta_total_fluid_20, J_theta_PIC_20,
            //               histVelField, E, B, 1, 1);
            // evalAzimuthal(nspec, nparts, ngpt, part, partWeight, num_density, vel_azimuthal, momentum_azimuthal, stress_azimuthal,
            //   ExB_component, diamagnetic_component, inertial_component,
            //   stressdiv_component, J_theta_total_fluid_21, J_theta_PIC_21,
            //   histVelField, E, B, 2, 1);
        }

        /* dump data */
        if (!(istep % ivdump))
        {
            ofstream outfile;

            int idump = istep / ivdump - 1;
            if (dumpField)
            {
                // electric potential field
                outfile.open("hist_phi.dat", ofstream ::out | ofstream ::app);
                for (int ieval = 0; ieval != nFieldEvals;)
                {
                    long vec_start = ieval * nxyz;
                    ++ieval;
                    outfile << (idump * ivdump + ieval * ivevalfield) * dt << " " << ngpt[0] << " " << ngpt[1] << " " << ngpt[2];
                    for (long nid = 0; nid < nxyz; nid++)
                        outfile << " " << histPhi.at(vec_start + nid);
                    outfile << "\n";
                }
                outfile.close();
                outfile.clear();
                histPhi.clear();

                // charge density field
                outfile.open("hist_chgDenst.dat", ofstream ::out | ofstream ::app);
                for (int ieval = 0; ieval != nFieldEvals;)
                {
                    long vec_start = ieval * nxyz * (nspec + 1);
                    ++ieval;
                    outfile << (idump * ivdump + ieval * ivevalfield) * dt << " " << ngpt[0] << " " << ngpt[1] << " " << ngpt[2] << " " << nspec + 1;
                    for (long nid = 0; nid < nxyz; nid++)
                        for (int ispec = 0; ispec <= nspec; ispec++)
                            outfile << " " << histChgDenst.at(vec_start + nid * (nspec + 1) + ispec);
                    outfile << "\n";
                }
                outfile.close();
                outfile.clear();
                histChgDenst.clear();

                // v and v^2 field
                outfile.open("hist_velField.dat", ofstream ::out | ofstream ::app);
                for (int ieval = 0; ieval != nFieldEvals;)
                {
                    long vec_start = ieval * nxyz * nspec * 6;
                    ++ieval;
                    outfile << (idump * ivdump + ieval * ivevalfield) * dt << " " << ngpt[0] << " " << ngpt[1] << " " << ngpt[2] << " " << nspec;
                    for (int ienty = 0; ienty < 6 * nspec * nxyz; ienty++)
                        outfile << " " << histVelField.at(vec_start + ienty);
                    outfile << "\n";
                }
                outfile.close();
                outfile.clear();
                histVelField.clear();

                // v^3 field
                outfile.open("new_hist_velField.dat", ofstream ::out | ofstream ::app);
                for (int ieval = 0; ieval != nFieldEvals;)
                {
                    long vec_start = ieval * nxyz * nspec * 5;
                    ++ieval;
                    outfile << (idump * ivdump + ieval * ivevalfield) * dt << " " << ngpt[0] << " " << ngpt[1] << " " << ngpt[2] << " " << nspec;
                    for (int ienty = 0; ienty < 5 * nspec * nxyz; ienty++)
                    {
                        outfile << " " << new_histVelField.at(vec_start + ienty);
                    }
                    outfile << "\n";
                }
                outfile.close();
                outfile.clear();
                new_histVelField.clear();

                // v^4 field
                outfile.open("new_new_hist_velField.dat", ofstream ::out | ofstream ::app);
                for (int ieval = 0; ieval != nFieldEvals;)
                {
                    long vec_start = ieval * nxyz * nspec * 3;
                    ++ieval;
                    outfile << (idump * ivdump + ieval * ivevalfield) * dt << " " << ngpt[0] << " " << ngpt[1] << " " << ngpt[2] << " " << nspec;
                    for (int ienty = 0; ienty < 3 * nspec * nxyz; ienty++)
                    {
                        outfile << " " << new_new_histVelField.at(vec_start + ienty);
                    }
                    outfile << "\n";
                }
                outfile.close();
                outfile.clear();
                new_new_histVelField.clear();

                // azimuthal components
                outfile.open("vel_azimuthal.dat", ofstream ::out | ofstream ::app);
                for (int ieval = 0; ieval != nFieldEvals;)
                {
                    long vec_start = ieval * nxyz * nspec;
                    ++ieval;
                    outfile << (idump * ivdump + ieval * ivevalfield) * dt << " " << ngpt[0] << " " << ngpt[1] << " " << ngpt[2] << " " << nspec;
                    for (int ienty = 0; ienty < nspec * nxyz; ienty++)
                    {
                        outfile << " " << vel_azimuthal.at(vec_start + ienty);
                    }
                    outfile << "\n";
                }
                outfile.close();
                outfile.clear();
                vel_azimuthal.clear();

                outfile.open("ExB_component.dat", ofstream ::out | ofstream ::app);
                for (int ieval = 0; ieval != nFieldEvals;)
                {
                    long vec_start = ieval * nxyz * nspec;
                    ++ieval;
                    outfile << (idump * ivdump + ieval * ivevalfield) * dt << " " << ngpt[0] << " " << ngpt[1] << " " << ngpt[2] << " " << nspec;
                    for (int ienty = 0; ienty < nspec * nxyz; ienty++)
                    {
                        outfile << " " << ExB_component.at(vec_start + ienty);
                    }
                    outfile << "\n";
                }
                outfile.close();
                outfile.clear();
                ExB_component.clear();

                outfile.open("diamagnetic_component.dat", ofstream ::out | ofstream ::app);
                for (int ieval = 0; ieval != nFieldEvals;)
                {
                    long vec_start = ieval * nxyz * nspec;
                    ++ieval;
                    outfile << (idump * ivdump + ieval * ivevalfield) * dt << " " << ngpt[0] << " " << ngpt[1] << " " << ngpt[2] << " " << nspec;
                    for (int ienty = 0; ienty < nspec * nxyz; ienty++)
                    {
                        outfile << " " << diamagnetic_component.at(vec_start + ienty);
                    }
                    outfile << "\n";
                }
                outfile.close();
                outfile.clear();
                diamagnetic_component.clear();

                outfile.open("inertial_component.dat", ofstream ::out | ofstream ::app);
                for (int ieval = 0; ieval != nFieldEvals;)
                {
                    long vec_start = ieval * nxyz * nspec;
                    ++ieval;
                    outfile << (idump * ivdump + ieval * ivevalfield) * dt << " " << ngpt[0] << " " << ngpt[1] << " " << ngpt[2] << " " << nspec;
                    for (int ienty = 0; ienty < nspec * nxyz; ienty++)
                    {
                        outfile << " " << inertial_component.at(vec_start + ienty);
                    }
                    outfile << "\n";
                }
                outfile.close();
                outfile.clear();
                inertial_component.clear();

                outfile.open("stressdiv_component.dat", ofstream ::out | ofstream ::app);
                for (int ieval = 0; ieval != nFieldEvals;)
                {
                    long vec_start = ieval * nxyz * nspec;
                    ++ieval;
                    outfile << (idump * ivdump + ieval * ivevalfield) * dt << " " << ngpt[0] << " " << ngpt[1] << " " << ngpt[2] << " " << nspec;
                    for (int ienty = 0; ienty < nspec * nxyz; ienty++)
                    {
                        outfile << " " << stressdiv_component.at(vec_start + ienty);
                    }
                    outfile << "\n";
                }
                outfile.close();
                outfile.clear();
                stressdiv_component.clear();

                outfile.open("J_theta_total_fluid.dat", ofstream ::out | ofstream ::app);
                for (int ieval = 0; ieval != nFieldEvals;)
                {
                    long vec_start = ieval * nxyz * nspec;
                    ++ieval;
                    outfile << (idump * ivdump + ieval * ivevalfield) * dt << " " << ngpt[0] << " " << ngpt[1] << " " << ngpt[2] << " " << nspec;
                    for (int ienty = 0; ienty < nspec * nxyz; ienty++)
                    {
                        outfile << " " << J_theta_total_fluid.at(vec_start + ienty);
                    }
                    outfile << "\n";
                }
                outfile.close();
                outfile.clear();
                J_theta_total_fluid.clear();

                outfile.open("J_theta_PIC.dat", ofstream ::out | ofstream ::app);
                for (int ieval = 0; ieval != nFieldEvals;)
                {
                    long vec_start = ieval * nxyz * nspec;
                    ++ieval;
                    outfile << (idump * ivdump + ieval * ivevalfield) * dt << " " << ngpt[0] << " " << ngpt[1] << " " << ngpt[2] << " " << nspec;
                    for (int ienty = 0; ienty < nspec * nxyz; ienty++)
                    {
                        outfile << " " << J_theta_PIC.at(vec_start + ienty);
                        ;
                    }
                    outfile << "\n";
                }
                outfile.close();
                outfile.clear();
                J_theta_PIC.clear();

                // outfile.open("J_theta_total_fluid_20.dat", ofstream ::out | ofstream ::app);
                // for (int ieval = 0; ieval != nFieldEvals;)
                // {
                //     long vec_start = ieval * nxyz * nspec;
                //     ++ieval;
                //     outfile << (idump * ivdump + ieval * ivevalfield) * dt << " " << ngpt[0] << " " << ngpt[1] << " " << ngpt[2] << " " << nspec;
                //     for (int ienty = 0; ienty < nspec * nxyz; ienty++)
                //     {
                //         outfile << " " << J_theta_total_fluid_20.at(vec_start + ienty);
                //     }
                //     outfile << "\n";
                // }
                // outfile.close();
                // outfile.clear();
                // J_theta_total_fluid_20.clear();

                // outfile.open("J_theta_PIC_20.dat", ofstream ::out | ofstream ::app);
                // for (int ieval = 0; ieval != nFieldEvals;)
                // {
                //     long vec_start = ieval * nxyz * nspec;
                //     ++ieval;
                //     outfile << (idump * ivdump + ieval * ivevalfield) * dt << " " << ngpt[0] << " " << ngpt[1] << " " << ngpt[2] << " " << nspec;
                //     for (int ienty = 0; ienty < nspec * nxyz; ienty++)
                //     {
                //         outfile << " " << J_theta_PIC_20.at(vec_start + ienty);
                //         ;
                //     }
                //     outfile << "\n";
                // }
                // outfile.close();
                // outfile.clear();
                // J_theta_PIC_20.clear();

                // outfile.open("J_theta_total_fluid_21.dat", ofstream ::out | ofstream ::app);
                // for (int ieval = 0; ieval != nFieldEvals;)
                // {
                //     long vec_start = ieval * nxyz * nspec;
                //     ++ieval;
                //     outfile << (idump * ivdump + ieval * ivevalfield) * dt << " " << ngpt[0] << " " << ngpt[1] << " " << ngpt[2] << " " << nspec;
                //     for (int ienty = 0; ienty < nspec * nxyz; ienty++)
                //     {
                //         outfile << " " << J_theta_total_fluid_21.at(vec_start + ienty);
                //         ;
                //     }
                //     outfile << "\n";
                // }
                // outfile.close();
                // outfile.clear();
                // J_theta_total_fluid_21.clear();

                // outfile.open("J_theta_PIC_21.dat", ofstream ::out | ofstream ::app);
                // for (int ieval = 0; ieval != nFieldEvals;)
                // {
                //     long vec_start = ieval * nxyz * nspec;
                //     ++ieval;
                //     outfile << (idump * ivdump + ieval * ivevalfield) * dt << " " << ngpt[0] << " " << ngpt[1] << " " << ngpt[2] << " " << nspec;
                //     for (int ienty = 0; ienty < nspec * nxyz; ienty++)
                //     {
                //         outfile << " " << J_theta_PIC_21.at(vec_start + ienty);
                //         ;
                //     }
                //     outfile << "\n";
                // }
                // outfile.close();
                // outfile.clear();
                // J_theta_PIC_21.clear();

                // azimuthal momentum and stress
                outfile.open("momentum_azimuthal.dat", ofstream ::out | ofstream ::app);
                for (int ieval = 0; ieval != nFieldEvals;)
                {
                    long vec_start = ieval * nxyz * nspec * 4;
                    ++ieval;
                    outfile << (idump * ivdump + ieval * ivevalfield) * dt << " " << ngpt[0] << " " << ngpt[1] << " " << ngpt[2] << " " << nspec;
                    for (int ienty = 0; ienty < 3 * nspec * nxyz; ienty++)
                    {
                        outfile << " " << momentum_azimuthal.at(vec_start + ienty);
                    }
                    outfile << "\n";
                }
                outfile.close();
                outfile.clear();
                momentum_azimuthal.clear();

                outfile.open("stress_azimuthal.dat", ofstream ::out | ofstream ::app);
                for (int ieval = 0; ieval != nFieldEvals;)
                {
                    long vec_start = ieval * nxyz * nspec * 4;
                    ++ieval;
                    outfile << (idump * ivdump + ieval * ivevalfield) * dt << " " << ngpt[0] << " " << ngpt[1] << " " << ngpt[2] << " " << nspec;
                    for (int ienty = 0; ienty < 3 * nspec * nxyz; ienty++)
                    {
                        outfile << " " << stress_azimuthal.at(vec_start + ienty);
                    }
                    outfile << "\n";
                }
                outfile.close();
                outfile.clear();
                stress_azimuthal.clear();
            }
        } //  end if(!(istep % ivdump))

        if (istep > 0 && (istep % ivdumpPart) == 0)
            dump_particles(istep, nparts, part_dump_box, pot);

        if (ivrestart > 0 && !(istep % ivrestart))
        {
            ostringstream oss;
            rstrt_id = rstrt_id == 0; // toggle restart file id
            oss << "restart." << rstrt_id;
            write_restart(istep, nparts[MAXNSPECTYPES], oss.str());
        }

        if (fdsample.sample(istep, nparts[MAXNSPECTYPES], part))
        {
            fdsample.postproc_field(dt, HYBRID, pot, fixedBV, Mat_A, poisson_solver, B);
        }

        vdf.sample(istep, nparts[MAXNSPECTYPES], part);
        edf.sample(istep, nparts[MAXNSPECTYPES], part);

    } // end while(istep < endstep)

    GET_TIME(tend);
    elapsed += (tend - tstart);
    cout << "time elapsed (injection/deposit/pmove/solveField/main_loop) = ("
         << inj_elapsed << "/" << depos_elapsed << "/"
         << pmove_elapsed << "/" << fieldSol_elapsed << "/"
         << elapsed << ")s." << endl;

    /* free memory allocated */
    for (long nid = 0; nid < nxyz; nid++)
    {
        delete[] chgDenst[nid];
        chgDenst[nid] = NULL;
    }
    delete[] chgDenst;
    delete[] part;
    delete[] E;
    delete[] B;
    if (EextOn)
        delete[] Eext;
    if (BextOn)
        delete[] Bext;

    return 0;
}
