#ifndef __ZFEMFORCES
#define __ZFEMFORCES

#include "zfem_structs.hpp"
#include "array.hpp"

//Compute internal force only
void InternalForce(double* U, mesh* msh, physics* phys, DirichletBCs* dbc, integration* quad, array<double> f);

//Compute internal force and tangent stiffness matrix (full)
void InternalForce(double* U, mesh* msh, physics* phys, DirichletBCs* dbc, integration* quad, array<double> f, array<double> df);

//Compute internal force and tangent stiffness matrix (sparse)
void InternalForce(double* U, mesh* msh, physics* phys, DirichletBCs* dbc, integration* quad, array<double> f, array<int> df_irow_cr, array<int> df_jcol, array<double> df_val, int nnz);

void BodyForce(double* b, mesh* msh, physics* phys, DirichletBCs* dbc, integration* quad, array<double> f);

void  NodalStress(double* U, mesh* msh, physics* phys, DirichletBCs* dbc, array<double> sigma);
#endif
