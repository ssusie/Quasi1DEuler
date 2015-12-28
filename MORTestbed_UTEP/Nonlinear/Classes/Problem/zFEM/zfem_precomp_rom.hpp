#ifndef __ZFEMPRECOMP
#define __ZFEMPRECOMP

#include "zfem_structs.hpp"
#include "array.hpp"
#include <iostream>

void  PrecompRomQuant(array<double> V, array<double> ub, int nY, mesh* msh, physics* phys, DirichletBCs* dbc, integration* quad, array<double> alpha, array<double> beta, array<double> gamma, array<double> omega);

void  PrecompRomQuant(array<double> V, array<double> ub, int nY, mesh* msh, physics* phys, DirichletBCs* dbc, integration* quad, array<double> alpha_lam, array<double> beta_lam, array<double> gamma_lam, array<double> omega_lam, array<double> alpha_mu, array<double> beta_mu, array<double> gamma_mu, array<double> omega_mu);

void  PrecompRomQuant_lam(array<double> V, array<double> ub, int nY, mesh* msh, physics* phys, DirichletBCs* dbc, integration* quad, array<double> alpha_lam, array<double> beta_lam, array<double> gamma_lam, array<double> omega_lam);

void  PrecompRomQuant_mu(array<double> V, array<double> ub, int nY, mesh* msh, physics* phys, DirichletBCs* dbc, integration* quad, array<double> alpha_mu, array<double> beta_mu, array<double> gamma_mu, array<double> omega_mu);
#endif
