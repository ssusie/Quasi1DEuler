#ifndef __ZFEMELEMS
#define __ZFEMELEMS

#include "zfem_structs.hpp"
#include <iostream>

void ShapeFunctions(integration * quad, int num, double *N, double **dN);

void ConstitutiveLaw(physics* phys, bool lin, double **H, double **P, double **S);

#endif
