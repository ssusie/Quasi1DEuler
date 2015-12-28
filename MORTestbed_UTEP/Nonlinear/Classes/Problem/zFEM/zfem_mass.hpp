#ifndef __ZFEMMASS
#define __ZFEMMASS

#include "zfem_structs.hpp"
#include "array.hpp"
#include <iostream>

//Mass Matrix stored directly in sparse format 
void MassMatrix(mesh* msh, physics* phys, integration* quad, array<int> M_irow_cr, array<int> M_jcol, array<double> M_val, int nnz);

//Mass Matrix times DBC vector
void MassTimesDbcVec(mesh* msh, physics* phys, integration* quad, DirichletBCs* dbc, array<double> a, array<double> Ma);



//(reduced) Mass Matrix (np x np)
void MassMatrix(mesh* msh, physics* phys, integration* quad, array<double> M);

//(reduced) Mass Matrix times a state vector
void MassTimesStateVec(mesh* msh, DirichletBCs* dbc, array<double> M, double* a, array<double> Ma);

//(reduced) Mass Matrix times DBC vector
void MassTimesDbcVec(mesh* msh, DirichletBCs* dbc, array<double> M, array<double> Ma);

//(reduced) Mass Matrix plus (full) matrix
void MassPlusMatrix(mesh* msh, DirichletBCs* dbc, array<double> M, array<double> K);

#endif
