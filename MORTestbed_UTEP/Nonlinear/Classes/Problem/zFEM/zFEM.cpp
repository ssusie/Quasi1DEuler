#include "mex.h"
#include "zfem_structs.hpp"
#include "convertML2Cstruct.hpp"
#include "zfem_forces.hpp"
#include "zfem_mass.hpp"
#include "zfem_elems.hpp"
#include "zfem_precomp_rom.hpp"
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
    //prhs[9]: ROB
    //prhs[10]: nY
    //prhs[11]: ub

    //plhs[0]: R = fint(U) - fext (ndof x 1)
    //plhs[1]: J = dR/dU (ndof x ndof)
    //plhs[2]: M = mass matrix

    mexPrintf("Number of inputs: %i \n",nrhs);
    mexPrintf("Number of outputs: %i \n",nlhs);

    mesh msh(prhs[2]);
    physics phys(prhs[3]);
    integration quad(prhs[4],msh);
    DirichletBCs dbc(prhs[5],msh);

    std::cout << "Z = [" << (quad.Z)(0,0) << ", " << (quad.Z)(0,1) << ", " << (quad.Z)(0,2) << ", " << (quad.Z)(0,3) << "]" <<std::endl;
    std::cout << "Z = [" << (quad.Z)(1,0) << ", " << (quad.Z)(1,1) << ", " << (quad.Z)(1,2) << ", " << (quad.Z)(1,3) << "]" <<std::endl;
    std::cout << "Z = [" << (quad.Z)(2,0) << ", " << (quad.Z)(2,1) << ", " << (quad.Z)(2,2) << ", " << (quad.Z)(2,3) << "]" <<std::endl;
    std::cout << "Z = [" << (quad.Z)(3,0) << ", " << (quad.Z)(3,1) << ", " << (quad.Z)(3,2) << ", " << (quad.Z)(3,3) << "]" <<std::endl;
    std::cout << "Z = [" << (quad.Z)(4,0) << ", " << (quad.Z)(4,1) << ", " << (quad.Z)(4,2) << ", " << (quad.Z)(4,3) << "]" <<std::endl;
    

    mexPrintf("np = %i\n",msh.np);
    mexPrintf("nel = %i\n",msh.nel);
    mexPrintf("nsd = %i\n",msh.nsd);
    mexPrintf("nen = %i\n",msh.nen);
    mexPrintf("ndof = %i\n",msh.ndof);

    mexPrintf("lam = %f\n",phys.lam);
    mexPrintf("mu = %f\n",phys.mu);
    mexPrintf("rho0 = %f\n",phys.rho0);

    mexPrintf("nint = %i\n",quad.nint);
    mexPrintf("ndbc = %i\n",dbc.ndbc);
 
    mexPrintf("LM size = [%i, %i]\n",msh.LM._dims[0],msh.LM._dims[1]);
    mexPrintf("LM size = [%i]\n",msh.LM._size);
    mexPrintf("LM data = %i\n",msh.LM(10,10));
    mexPrintf("LM alloc = [%i]\n",msh.LM._alloc);
    //msh.LM.printFormatted();
    //msh.IEN.printFormatted();
    //msh.ID.printFormatted();
    //msh.X.printFormatted();
    //dbc.loc.printFormatted();

    double* U = mxGetPr(prhs[0]);

    int nnzeros = (int) mxGetScalar(prhs[8]);
    int nY = (int) mxGetScalar(prhs[10]);
 
    std::cout << "nnzeros = " << nnzeros << std::endl;
    std::cout << "nY = " << nY << std::endl;

    size_t dims[4];
    int* irow_ptr = (int*) mxGetData(prhs[6]);
    int* jcol_ptr = (int*) mxGetData(prhs[7]);
    double* ub_ptr = mxGetPr(prhs[11]);
    array<int> irow;
    array<int> jcol;
    array<double> ub;

    dims[0]=msh.ndof; dims[1]=1; dims[2]=1; dims[3]=1;
    irow.setarraysize(1, dims); irow.setreference(irow_ptr);
    dims[0]=nnzeros; dims[1]=1;
    jcol.setarraysize(1, dims); jcol.setreference(jcol_ptr);
    dims[0]=msh.ndof; dims[1]=1;
    ub.setarraysize(1, dims); ub.setreference(ub_ptr);

    double* V_ptr = mxGetPr(prhs[9]);
    dims[0]=msh.ndof; dims[1]=nY;
    array<double> V;
    V.setarraysize(2, dims); V.setreference(V_ptr); 

    //Output force
    plhs[0] = mxCreateDoubleMatrix(msh.ndof,1,mxREAL);
    double* fptr = mxGetPr(plhs[0]);
    
    //Output Jacobian
    plhs[1] = mxCreateDoubleMatrix(msh.ndof,msh.ndof,mxREAL);
    double* dfptr = mxGetPr(plhs[1]);

    //Output Mass
    plhs[2] = mxCreateDoubleMatrix(msh.np,msh.np,mxREAL);
    double* Mptr = mxGetPr(plhs[2]);

    //Output Mass times vector
    plhs[3] = mxCreateDoubleMatrix(msh.ndof,1,mxREAL);
    double* Maptr = mxGetPr(plhs[3]);

    //Output Mass times vector
    plhs[4] = mxCreateDoubleMatrix(msh.ndof,1,mxREAL);
    double* Madbcptr = mxGetPr(plhs[4]);

    //Output Mass times vector
    plhs[5] = mxCreateDoubleMatrix(msh.ndof,msh.ndof,mxREAL);
    double* MpKptr = mxGetPr(plhs[5]);

    //Output sparse stiffness
    plhs[6] = mxCreateDoubleMatrix(nnzeros,1,mxREAL);
    double* spdf_ptr = mxGetPr(plhs[6]);

    //Output sparse stiffness
    plhs[7] = mxCreateDoubleMatrix(nnzeros,1,mxREAL);
    double* spM_ptr = mxGetPr(plhs[7]);

    //Output sparse stiffness
    plhs[8] = mxCreateDoubleMatrix(nY,1,mxREAL);
    double* alpha_ptr = mxGetPr(plhs[8]);

    //Output sparse stiffness
    plhs[9] = mxCreateDoubleMatrix(nY,nY,mxREAL);
    double* beta_ptr = mxGetPr(plhs[9]);

    //Output sparse stiffness
    int idims[4];
    idims[0]=nY; idims[1]=nY; idims[2]=nY; idims[3]=nY;
    plhs[10] = mxCreateNumericArray(3,idims,mxDOUBLE_CLASS,mxREAL);
    double* gamma_ptr = mxGetPr(plhs[10]);

    //Output sparse stiffness
    plhs[11] = mxCreateNumericArray(4,idims,mxDOUBLE_CLASS,mxREAL);
    double* omega_ptr = mxGetPr(plhs[11]);

    array<double> f;
    array<double> df;
    array<double> M;
    array<double> Ma;
    array<double> Madbc;
    array<double> MpK;
    array<double> spdf;
    array<double> spM;
    array<double> alpha;
    array<double> beta;
    array<double> gamma;
    array<double> omega;

    dims[0]=msh.ndof; dims[1]=msh.ndof;
    f.setarraysize(1, dims); f.setreference(fptr);
    df.setarraysize(2,dims); df.setreference(dfptr);
    dims[0]=msh.np; dims[1]=msh.np;
    M.setarraysize(2,dims); M.setreference(Mptr);
    dims[0]=msh.ndof; dims[1]=msh.ndof;
    Ma.setarraysize(1,dims); Ma.setreference(Maptr);
    Madbc.setarraysize(1,dims); Madbc.setreference(Madbcptr);
    MpK.setarraysize(2,dims); MpK.setreference(MpKptr);
    dims[0]=nnzeros; dims[1]=1;
    spdf.setarraysize(1,dims); spdf.setreference(spdf_ptr);
    spM.setarraysize(1,dims); spM.setreference(spM_ptr);

    dims[0]=nY; dims[1]=nY; dims[2]=nY; dims[3]=nY;
    alpha.setarraysize(1,dims); alpha.setreference(alpha_ptr);
    beta.setarraysize(2,dims);  beta.setreference(beta_ptr);
    gamma.setarraysize(3,dims); gamma.setreference(gamma_ptr);
    omega.setarraysize(4,dims); omega.setreference(omega_ptr);

    InternalForce(U,&msh,&phys,&dbc,&quad,f);
    MassMatrix(&msh,&phys,&quad,M);
    MassTimesStateVec(&msh,&dbc,M,U,Ma);
    MassTimesDbcVec(&msh,&dbc,M,Madbc);

    for (int i=0; i<msh.ndof; ++i)
       for (int j=0; j < msh.ndof; ++j)
          MpK(i,j) = df(i,j);
    MassPlusMatrix(&msh,&dbc,M,MpK);

    InternalForce(U,&msh,&phys,&dbc,&quad,f,irow,jcol,spdf,nnzeros);
    MassMatrix(&msh,&phys,&quad,irow,jcol,spM,nnzeros);

    std::cout << alpha._size << std::endl;
    std::cout << alpha._dims[0] << std::endl << std::endl;

    std::cout << beta._size << std::endl;
    std::cout << beta._dims[0] << "  " << beta._dims[1] << std::endl << std::endl;

    std::cout << gamma._size << std::endl;
    std::cout << gamma._dims[0] << "  " << gamma._dims[1] << "  " << gamma._dims[2] << std::endl << std::endl;

    std::cout << omega._size << std::endl;
    std::cout << omega._dims[0] << "  " << omega._dims[1] << "  " << omega._dims[2] << "  " << omega._dims[3] << std::endl << std::endl;

    /*std::cout << "alpha = " << std::endl;
    for (int i=0; i < nY; ++i)
       std::cout << alpha(i) << std::endl;

    std::cout << "beta = " << std::endl;
    for (int i=0; i < nY; ++i)
       for (int j=0; j < nY; ++j)
          std::cout << beta(i,j) << std::endl;

    std::cout << "gamma = " << std::endl;
    for (int i=0; i < nY; ++i)
       for (int j=0; j < nY; ++j)
          for (int k=0; k < nY; ++k)
             std::cout << gamma(i,j,k) << std::endl;

    std::cout << "omega = " << std::endl;
    for (int i=0; i < nY; ++i)
       for (int j=0; j < nY; ++j)
          for (int k=0; k < nY; ++k)
             for (int l=0; l < nY; ++l)
                std::cout << omega(i,j,k,l) << std::endl;*/

    return;
    for (int i=0; i < msh.ndof; ++i) {
       for (int j = 0; j < nY; ++j)
          std::cout << V(i,j) << " ";
       std::cout << ub(i) <<  std::endl;
    }

    PrecompRomQuant(V, ub, nY, &msh, &phys, &dbc, &quad,  alpha, beta, gamma, omega);


    //InternalForce(U,msh,phys,dbc,quad,f,df);
    //BodyForce(b,msh,dbc,quad,f);
    //also need to add nodal forces
}
