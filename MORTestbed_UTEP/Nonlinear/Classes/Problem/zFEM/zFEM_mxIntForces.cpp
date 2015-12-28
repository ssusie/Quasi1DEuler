#include "mex.h"
#include "zfem_structs.hpp"
#include "convertML2Cstruct.hpp"
#include "zfem_forces.hpp"
#include "zfem_elems.hpp"
#include "matrix.h"

void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[]) 
{
    //prhs[0]: U - current state vector (ndof x 1)
    //prhs[1]: t - scalar indicating current time
    //prhs[2]: msh structure
    //prhs[3]: phys structure 
    //prhs[4]: integration structure
    //prhs[5]: DirichletBCs structure
    //prhs[6]: irow_sr
    //prhs[7]: jcol
    //prhs[8]: nnz in jacobian

    //plhs[0]: R = fint(U) - fext (ndof x 1)
    //plhs[1]: J = dR/dU (ndof x ndof)
 
    //Extract state vector
    double* U = (double*)mxGetPr(prhs[0]);

    //Unpack structures
    mesh msh(prhs[2]);
    physics* phys = new physics[msh.nummat];
    for (int i=0; i<msh.nummat; ++i)
       phys[i].setProperties(prhs[3],i);
    //physics phys(prhs[3]);
    integration quad(prhs[4],msh);
    DirichletBCs dbc(prhs[5],msh);

    //std::cout <<  phys[0].lam << std::endl;
    //std::cout <<  phys[1].lam << std::endl;
    //std::cout <<  phys[0].mu << std::endl;
    //std::cout <<  phys[1].mu << std::endl;
    //std::cout <<  phys[0].rho0 << std::endl;
    //std::cout <<  phys[1].rho0 << std::endl;

    //Output force
    plhs[0] = mxCreateDoubleMatrix(msh.ndof,1,mxREAL);
    double* fptr = (double*) mxGetPr(plhs[0]);
   
    //Allocate force vector 
    size_t dims[1];
    array<double> f;
    dims[0]=msh.ndof;
    f.setarraysize(1, dims); f.setreference(fptr);

    if (nlhs == 1)
       InternalForce(U,&msh,phys,&dbc,&quad,f);
    else {
       int nnzeros = (int) mxGetScalar(prhs[8]);
       int* irow_ptr = (int*) mxGetData(prhs[6]);
       int* jcol_ptr = (int*) mxGetData(prhs[7]);
       array<int> irow;
       array<int> jcol;

       dims[0]=msh.ndof; dims[1]=1;
       irow.setarraysize(1, dims); irow.setreference(irow_ptr);
       dims[0]=nnzeros; dims[1]=1;
       jcol.setarraysize(1, dims); jcol.setreference(jcol_ptr);

       //Output sparse stiffness
       plhs[1] = mxCreateDoubleMatrix(nnzeros,1,mxREAL);
       double* spdf_ptr = mxGetPr(plhs[1]);
       array<double> spdf;

       dims[0]=nnzeros;
       spdf.setarraysize(1,dims); spdf.setreference(spdf_ptr);

       InternalForce(U,&msh,phys,&dbc,&quad,f,irow,jcol,spdf,nnzeros);
   }
   delete[] phys;
}
