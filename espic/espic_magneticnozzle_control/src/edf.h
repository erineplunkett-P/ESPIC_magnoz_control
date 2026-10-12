#ifndef EDF_H
#define EDF_H

#include <vector>
#include <string>
#include "part_spec.h"
#include "global.h"

using std::vector;

class EnergyDist
{
public:
    // Constructor
    EnergyDist(const std::string &, int, int, int[3]);

    // Destructor
    ~EnergyDist();

    void sample(int, long, double (*)[MAXNPARTPROPS]); // count histogram for edf

private:
    int step_start;
    int step_now;
    int nevery_sample;
    int nevery_restart;
    int nevery_save;

    int nspec;               // total # of species
    vector<bool> spec_exist; // true if needs to sample this species
    int dim_cell[3];         // nx, ny, nz of cells
    int ncells;              // total # of cells
    int ncells_sample;

    int *cell; // map cell -> sampled-cellid or -1

    // One energy histogram per species (if enabled)
    int map_histid[MAXNSPECTYPES];

    // Energy mode per histogram
    enum EnergyMode
    {
        ENERGY_FULL3 = 0,  // 0.5 m (vz^2 + vr^2 + vtheta^2)
        ENERGY_AXIS_ZR = 1, // 0.5 m (vz^2 + vr^2)
        ENERGY_BOTH = 2     // output both ENERGY_FULL3 and ENERGY_AXIS_ZR
    };
    vector<int> emode; // size = nhist, values in EnergyMode

    vector<double> emin;
    vector<double> emax;
    vector<double> de;
    vector<int> nbin;

    vector<int> **hist_array_full; // hist_array_full[cellid][histid][ibin]
    vector<int> **hist_array_axis; // hist_array_axis[cellid][histid][ibin]
    int **ncount;             // ncount[cellid][histid]

    void init(const std::string &);
    void write_restart();
    void read_restart();
    void save_edf_epf();
    void clear_edf();
};

#endif
