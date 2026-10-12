#ifndef VDF_H
#define VDF_H

#include <vector>
#include <string>
#include "part_spec.h"
#include "global.h"

using std::vector;

class VelDist
{
public:
    // Constructor
    VelDist(const std::string &, int, int, int[3]);

    // Destructor
    ~VelDist();

    void sample(int, long, double (*)[MAXNPARTPROPS]); // count histogram for vdf

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

    int *cell; // true if needs to sample in this cell

    vector<int> ivel[MAXNSPECTYPES]; // ivel[ispec][i] = 0(x), 1(y) or 2(z)
    int map_histid[MAXNSPECTYPES][3];
    vector<double> vmin;
    vector<double> vmax;
    vector<double> dv;
    vector<int> nbin;
    vector<int> **hist_array;
    int **ncount;

    void init(const std::string &);
    void write_restart();
    void read_restart();
    void save_vdf();
    void clear_vdf();
};

#endif
