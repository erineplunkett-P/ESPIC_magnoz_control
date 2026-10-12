#ifndef SAMPLEFIELD_H
#define SAMPLEFIELD_H

#include <string>
#include "solve_field.h"
#include "part_spec.h"
#include "global.h"

class SampleField
{
public:
    // Constructor
    SampleField(const std::string &, int, int, int[3]);

    // Destructor
    ~SampleField();

    struct Sample
    {
        double NSUM;

        // Vz*Vr
        // double VZVRSUM;

        double VSUM[3];
        double V2SUM[3];

        // --- New raw moment accumulators (axisymmetric use: v[0], v[1]) ---
        // Matches evalNewHistVelField ordering:
        // [0] v0^3, [1] v1^3, [2] v0^2*v1, [3] v0*v1^2, [4] v0*v1
        double V3SUM[5];

        // Matches evalNewNewHistVelField ordering:
        // [0] v0^4, [1] v1^4, [2] v0^2*v1^2
        double V4SUM[3];
    };

    // public methods
    bool sample(int, long, double (*)[MAXNPARTPROPS]); // sample flow field
    void postproc_field(const double, const bool, const Vec &, const MapF &, SpMat &, PtrPoissonSolver, double (*B)[3] = nullptr, int eqSelect = 0);

private:
    int step_start;
    int step_now;
    int nevery_sample;
    int nevery_write;
    int nevery_postproc;
    uint32_t nsamples;
    int nspec;
    int ngpt[3];
    long nnodes;
    Sample *spfield;

    // private method
    void read_input(const std::string &);
    bool read_sample(int);
    void write_sample();
    void compute_moment_sum(Sample *, double *, double);
    // For the ease of the calculation of averaged azimuthal quantities
    void compute_time_averaged_azimuthal(double (*E)[3],
                                         double (*B)[3],
                                         int eqSelect,
                                         std::vector<double> &number_density,
                                         std::vector<double> &vel_azimuthal,
                                         std::vector<double> &momentum_azimuthal,
                                         std::vector<double> &stress_azimuthal,
                                         std::vector<double> &ExB_component,
                                         std::vector<double> &diamagnetic_component,
                                         std::vector<double> &inertial_component,
                                         std::vector<double> &stressdiv_component,
                                         std::vector<double> &Je_theta_total_fluid_fullmomentum,
                                         std::vector<double> &Je_theta_PIC);

    void write_time_averaged_azimuthal(const double dt,
                                       const std::vector<double> &number_density,
                                       const std::vector<double> &vel_azimuthal,
                                       const std::vector<double> &momentum_azimuthal,
                                       const std::vector<double> &stress_azimuthal,
                                       const std::vector<double> &ExB_component,
                                       const std::vector<double> &diamagnetic_component,
                                       const std::vector<double> &inertial_component,
                                       const std::vector<double> &stressdiv_component,
                                       const std::vector<double> &Je_theta_total_fluid_fullmomentum,
                                       const std::vector<double> &Je_theta_PIC);
};

#endif
