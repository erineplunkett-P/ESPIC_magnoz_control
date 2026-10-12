#include "es_PIC.h"

void eval_EEDF( 
    int nspec,
    const long *nparts,
    const int ngpt[3],
    const double (*part)[MAXNPARTPROPS],
    const double *particleWeight, 
    const double EEDF_setup[3],              // [E_min, E_max, nbins]
    const int EEDF_region[4],                // [start_x, start_y, end_x, end_y]
    std::vector<double> &EEDF, 
    std::vector<double> &EEPF, 
    std::vector<double> &EEDF_centers
)
{
    // ----------------------------------
    // INITIALIZATION : dE, H, E_center 
    // ----------------------------------
    long npartsTot = nparts[MAXNSPECTYPES];
    double partX[3], partV[3];

    const double EMIN = EEDF_setup[0]; const double EMAX = EEDF_setup[1]; const int NBINS = (int)EEDF_setup[2]; 
    const int start_x = EEDF_region[0]; const int start_y = EEDF_region[1]; const int end_x = EEDF_region[2]; const int end_y = EEDF_region[3]; 

    double dE = (EMAX-EMIN) / (double)NBINS; 
    int region_area = (end_x - start_x) * (end_y - start_y); 
    const double J_to_eV = 1.0 / 1.602176634e-19; 

    double *internal_H = new double[NBINS]; 
    double *internal_E_center = new double[NBINS]; 
    double *EEDF_hist = new double[NBINS]; 
    double *EEPF_hist = new double[NBINS]; 

    memset(internal_H, 0, NBINS * sizeof(double)); 
    memset(internal_E_center, 0, NBINS * sizeof(double)); 
    memset(EEDF_hist, 0, NBINS * sizeof(double)); 
    memset(EEPF_hist, 0, NBINS * sizeof(double)); 

    for(int i = 0; i < NBINS; i++) {
        internal_H[i] = 0.0; 
        internal_E_center[i] = EMIN + (i + 0.5) * dE; 
    }

    int ix, iy, bin; 
    double xpm, ypm, me, vx, vy, vz, v2, ke_norm; 

    // ---------------------------------------------------------
    // PARTICLE LOOP : compute partX, partV, and kinetic energy 
    // ---------------------------------------------------------
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
        {
            partX[0] = domain[1] - smallValue;
        }
        if (partX[1] >= domain[3])
        {
            partX[1] = domain[3] - smallValue;
        }
        if (partX[2] >= domain[5])
        {
            partX[2] = domain[5] - smallValue;
        }
        if (partX[0] <= domain[0])
        {
            partX[0] = domain[0] + smallValue;
        }
        if (partX[1] <= domain[2])
        {
            partX[1] = domain[2] + smallValue;
        }
        if (partX[2] <= domain[4])
        {
            partX[2] = domain[4] + smallValue;
        }

        xpm = (partX[0] - domain[0]) * dxinv[0]; // z
        ypm = (partX[1] - domain[2]) * dxinv[1]; // r
        ix = (int)xpm;
        iy = (int)ypm;

        if(itype != 0) continue; // itype is either 0 or 1
        if(ix < start_x || ix > end_x || iy < start_y || iy > end_y) continue; 

        me = spec[itype][specmass]; 
        vx = partV[0]; vy = partV[1]; vz = partV[2]; 

        v2 = vx*vx + vy*vy + vz*vz; 
        ke_norm = 0.5 * me * v2; 

        if(ke_norm < EMIN || ke_norm > EMAX) continue; 
        bin = (int)floor((ke_norm - EMIN) / dE); 
        if(bin >= 0 && bin < NBINS) {
            internal_H[bin] += particleWeight[itype] * sqrt(internal_E_center[bin]); 
        }
        
    } // END PARTICLE LOOP 

    // -------------------------------------------------------------------
    // BUILD EEDF AND EEPF 
    // -------------------------------------------------------------------
    double sumH = 0.0; 
    for(int i = 0; i < NBINS; i++) sumH += internal_H[i]; 

    // build EEDF: normalized per region area and energy (normalized energy units) 
    for(int i = 0; i < NBINS; i++) {
        EEDF_hist[i] = internal_H[i] / (region_area * dE); // particles per cell per kB T_e0
    }

    // build EEPF 
    for(int i = 0; i < NBINS; i++) {
        if(internal_E_center[i] > 0.0) {
            EEPF_hist[i] = EEDF_hist[i] / sqrt(internal_E_center[i]); 
        }
        else EEPF_hist[i] = 0.0; 
    }

    // ----------------------------------------------
    // PACK RESULTS TO VECTORS 
    // ----------------------------------------------
    EEDF.clear(); 
    EEPF.clear(); 
    EEDF_centers.clear(); 

    EEDF.reserve(NBINS);
    EEPF.reserve(NBINS); 
    EEDF_centers.reserve(NBINS); 

    for(int i = 0; i < NBINS; i++) {
        EEDF.push_back(EEDF_hist[i]); 
        EEPF.push_back(EEPF_hist[i]); 
        EEDF_centers.push_back(internal_E_center[i]); // normalized energy units 
    }

    delete[] internal_H; 
    delete[] internal_E_center; 
    delete[] EEDF_hist; 
    delete[] EEPF_hist; 

    return; 
}