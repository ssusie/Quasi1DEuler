#include "zfem_structs.hpp"
#include "convertML2Cstruct.hpp"
#include "zfem_elems.hpp"
#include "zfem_precomp_rom.hpp"
#include "matrix.h"

void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[])
{
    //prhs[0]: V - reduced basis
    //prhs[1]: ub - basis offset vector
    //prhs[2]: nY - ROB size
    //prhs[3]: msh structure
    //prhs[4]: phys structure 
    //prhs[5]: integration structure
    //prhs[6]: DirichletBCs structure
    // if nrhs = 8 && nlhs = 4 (indicates only to compute lambda OR mu term, not both)
    //prhs[7]: ind (= 1 to compute lambda, = 2 to compute mu)

    //plhs[0]: alpha ( alpha_lam if nlhs = 8)
    //plhs[1]: beta  ( beta_lam  if nlhs = 8)
    //plhs[2]: gamma ( gamma_lam if nlhs = 8)
    //plhs[3]: omega ( omega_lam if nlhs = 8)
    // if nlhs = 8
    //plhs[4]: alpha_mu
    //plhs[5]: beta_mu
    //plhs[6]: gamma_mu
    //plhs[7]: omega_mu

    //Unpack structures
    mesh msh(prhs[3]);
    physics* phys = new physics[msh.nummat];
    for (int i=0; i<msh.nummat; ++i)
       phys[i].setProperties(prhs[4],i);
    //physics phys(prhs[4]);
    integration quad(prhs[5],msh);
    DirichletBCs dbc(prhs[6],msh);

    //Get pointer to ROB and offset
    double* V_ptr  = (double*)mxGetPr(prhs[0]);
    double* ub_ptr = (double*)mxGetPr(prhs[1]);

    //Get nY
    int nY = (int) mxGetScalar(prhs[2]);

    //Create input arrays 
    size_t dims[4];   dims[0]=msh.ndof; dims[1]=nY;
    array<double> V, ub;
    V.setarraysize(2, dims);  V.setreference(V_ptr);
    ub.setarraysize(1, dims); ub.setreference(ub_ptr);

    //Create output mxArrays
    int idims[4];  idims[0]=nY; idims[1]=nY; idims[2]=nY; idims[3]=nY;
    plhs[0] = mxCreateDoubleMatrix(nY,1,mxREAL);
    plhs[1] = mxCreateDoubleMatrix(nY,nY,mxREAL);
    plhs[2] = mxCreateNumericArray(3,idims,mxDOUBLE_CLASS,mxREAL);
    plhs[3] = mxCreateNumericArray(4,idims,mxDOUBLE_CLASS,mxREAL);

    //Pointers to output data
    double* alpha_ptr = (double*) mxGetPr(plhs[0]);
    double* beta_ptr  = (double*) mxGetPr(plhs[1]);
    double* gamma_ptr = (double*) mxGetPr(plhs[2]);
    double* omega_ptr = (double*) mxGetPr(plhs[3]);

    array<double> alpha, beta, gamma, omega;
    dims[0]=nY; dims[1]=nY; dims[2]=nY; dims[3]=nY;
    alpha.setarraysize(1,dims); alpha.setreference(alpha_ptr);
    beta.setarraysize(2,dims);  beta.setreference(beta_ptr);
    gamma.setarraysize(3,dims); gamma.setreference(gamma_ptr);
    omega.setarraysize(4,dims); omega.setreference(omega_ptr);

    if (nlhs == 4) {
       if (nrhs == 7)
         PrecompRomQuant(V, ub, nY, &msh, phys, &dbc, &quad,  alpha, beta, gamma, omega);
       
       if (nrhs == 8) {
          int ind = (int) mxGetScalar(prhs[7]);
          if (ind == 1)
             PrecompRomQuant_lam(V, ub, nY, &msh, phys, &dbc, &quad,  alpha, beta, gamma, omega);
          else if (ind == 2)
             PrecompRomQuant_mu(V, ub, nY, &msh, phys, &dbc, &quad,  alpha, beta, gamma, omega);
       }
    } else if (nlhs == 8) {

       //Create output mxArrays
       plhs[4] = mxCreateDoubleMatrix(nY,1,mxREAL);
       plhs[5] = mxCreateDoubleMatrix(nY,nY,mxREAL);
       plhs[6] = mxCreateNumericArray(3,idims,mxDOUBLE_CLASS,mxREAL);
       plhs[7] = mxCreateNumericArray(4,idims,mxDOUBLE_CLASS,mxREAL);

       double* alpha_mu_ptr = (double*) mxGetPr(plhs[4]);
       double* beta_mu_ptr  = (double*) mxGetPr(plhs[5]);
       double* gamma_mu_ptr = (double*) mxGetPr(plhs[6]);
       double* omega_mu_ptr = (double*) mxGetPr(plhs[7]);

       array<double> alpha_mu, beta_mu, gamma_mu, omega_mu;
       dims[0]=nY; dims[1]=nY; dims[2]=nY; dims[3]=nY;
       alpha_mu.setarraysize(1,dims); alpha_mu.setreference(alpha_mu_ptr);
       beta_mu.setarraysize(2,dims);  beta_mu.setreference(beta_mu_ptr);
       gamma_mu.setarraysize(3,dims); gamma_mu.setreference(gamma_mu_ptr);
       omega_mu.setarraysize(4,dims); omega_mu.setreference(omega_mu_ptr);

       alpha_mu.zeroout(); beta_mu.zeroout(); gamma_mu.zeroout(); omega_mu.zeroout();

       PrecompRomQuant(V, ub, nY, &msh, phys, &dbc, &quad, alpha, beta, gamma, omega, alpha_mu, beta_mu, gamma_mu, omega_mu);
    }
    delete[] phys;
}

/*    std::cout << alpha._size << std::endl;
    std::cout << alpha._dims[0] << std::endl << std::endl;

    std::cout << beta._size << std::endl;
    std::cout << beta._dims[0] << "  " << beta._dims[1] << std::endl << std::endl;

    std::cout << gamma._size << std::endl;
    std::cout << gamma._dims[0] << "  " << gamma._dims[1] << "  " << gamma._dims[2] << std::endl << std::endl;

    std::cout << omega._size << std::endl;
    std::cout << omega._dims[0] << "  " << omega._dims[1] << "  " << omega._dims[2] << "  " << omega._dims[3] << std::endl << std::endl;


    std::cout << "nY = " << nY << std::endl;
    for (int i=0; i < msh.ndof; ++i) {
       for (int j = 0; j < nY; ++j)
          std::cout << V(i,j) << " ";
       std::cout << ub(i) <<  std::endl;
    }

    msh.LM.printFormatted();
    msh.IEN.printFormatted();*/
