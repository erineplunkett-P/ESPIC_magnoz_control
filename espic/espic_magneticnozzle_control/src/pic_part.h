#ifndef PICPART_H
#define PICPART_H

#include <cstddef>
#include "objects_in_domain.h"

class Sources;

bool initPart(const int flag_read_restart, long *nparts, double *partWeight);

void partMove(const double half_dt,
              const int ngpt[3],
              long *nparts,
              double (*part)[MAXNPARTPROPS],
              double (*E)[3], double (*B)[3],
              const ObjectList &objlist,
              Sources &sourcelist,
              const int nthrds);

void part2meshChg(int nspec, long *nparts,
                  const int ngpt[3],
                  double *partWeight,
                  double **chgDenst,
                  const int nthrds);

void adjustPartMem();

extern int *vacMap;

#endif
