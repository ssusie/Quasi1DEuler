#include "mex.h"
#include "zfem_structs.hpp"
#include "convertML2Cstruct.hpp"
#include "zfem_forces.hpp"
#include "zfem_elems.hpp"
#include "matrix.h"

void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[]) 
{
    //prhs[0]: b - body forces (ndof x 1)
    //prhs[1]: t - scalar indicating current time
    //prhs[2]: msh structure
    //prhs[3]: phys structure 
    //prhs[4]: integration structure
    //prhs[5]: DirichletBCs structure

    //plhs[0]: fext (ndof x 1)
 
    //Extract state vector
    double* b = mxGetPr(prhs[0]);

    //Unpack structures
    mesh msh(prhs[2]);

    physics* phys = new physics[msh.nummat];
    for (int i=0; i<msh.nummat; ++i)
       phys[i].setProperties(prhs[3],i);
    //physics phys(prhs[3]);
    
    integration quad(prhs[4],msh);
    DirichletBCs dbc(prhs[5],msh);

    //Output force
    plhs[0] = mxCreateDoubleMatrix(msh.ndof,1,mxREAL);
    double* fptr = mxGetPr(plhs[0]);
   
    //Allocate force vector 
    size_t dims[1];
    array<double> f;
    dims[0]=msh.ndof;
    f.setarraysize(1, dims); f.setreference(fptr);

    BodyForce(b,&msh,phys,&dbc,&quad,f);
    delete[] phys;
}
