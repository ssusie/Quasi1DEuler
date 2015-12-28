#include "mex.h"
#include "zfem_structs.hpp"
#include "convertML2Cstruct.hpp"
#include "zfem_mass.hpp"
#include "zfem_elems.hpp"
#include "matrix.h"

void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[]) 
{
    //prhs[0]: msh structure
    //prhs[1]: phys structure 
    //prhs[2]: integration structure
    //prhs[3]: irow_sr
    //prhs[4]: jcol
    //prhs[5]: nnz in jacobian

    //plhs[0]: M = mass matrix

    mesh msh(prhs[0]);
    physics* phys = new physics[msh.nummat];
    for (int i=0; i<msh.nummat; ++i)
       phys[i].setProperties(prhs[1],i);
    //physics phys(prhs[1]);
    integration quad(prhs[2],msh);

    int nnzeros = (int) mxGetScalar(prhs[5]);

    size_t dims[2];
    int* irow_ptr = (int*) mxGetData(prhs[3]);
    int* jcol_ptr = (int*) mxGetData(prhs[4]);
    array<int> irow;
    array<int> jcol;

    dims[0]=msh.ndof; dims[1]=1;
    irow.setarraysize(1, dims); irow.setreference(irow_ptr);
    dims[0]=nnzeros; dims[1]=1;
    jcol.setarraysize(1, dims); jcol.setreference(jcol_ptr);

    //Output sparse stiffness
    plhs[0] = mxCreateDoubleMatrix(nnzeros,1,mxREAL);
    double* spM_ptr = mxGetPr(plhs[0]);

    array<double> spM;

    dims[0]=nnzeros; dims[1]=1;
    spM.setarraysize(1,dims); spM.setreference(spM_ptr);

    MassMatrix(&msh,phys,&quad,irow,jcol,spM,nnzeros);

    delete[] phys;
}
