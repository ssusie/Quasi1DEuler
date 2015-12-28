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
    //prhs[3]: DirichletBCs structure
    //prhs[4]: BC acceleration

    //plhs[0]: M*a

    mesh msh(prhs[0]);
    physics* phys = new physics[msh.nummat];
    for (int i=0; i<msh.nummat; ++i)
       phys[i].setProperties(prhs[1],i);
    //physics phys(prhs[1]);
    integration quad(prhs[2],msh);
    DirichletBCs dbc(prhs[3],msh);

    size_t dims[1];
    dims[0]=dbc.ndbc;

    array<double> a;
    double* a_ptr = (double*) mxGetData(prhs[4]);
    a.setarraysize(1, dims); a.setreference(a_ptr);

    //Output M*a_dbc
    plhs[0] = mxCreateDoubleMatrix(msh.ndof,1,mxREAL);
    double* Ma_ptr = mxGetPr(plhs[0]);

    array<double> Ma;
    dims[0]=msh.ndof;
    Ma.setarraysize(1,dims); Ma.setreference(Ma_ptr);

    MassTimesDbcVec(&msh,phys,&quad,&dbc,a,Ma);
    delete[] phys;
}
