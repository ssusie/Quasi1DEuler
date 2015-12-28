#ifndef __ZFEMSTRUCTS
#define __ZFEMSTRUCTS

#include "array.hpp"
#include "convertML2Cstruct.cpp"
#include <iostream>

class mesh {
 public:
    int np, nel, nsd, nen, ndof, nummat;
    array<int> mat;
    array<int> LM;  
    array<int> IEN;
    array<int> ID;
    array<double> X;

 public:
   explicit mesh(mxArray const* mlmesh) {

      getscalarsafe(mlmesh,"np", np);
      getscalarsafe(mlmesh,"nel", nel);
      getscalarsafe(mlmesh,"nsd", nsd);
      getscalarsafe(mlmesh,"nen", nen);
      getscalarsafe(mlmesh,"ndof", ndof);
      getscalarsafe(mlmesh,"nummat", nummat);

      size_t dims[2];
      dims[0] = nel; dims[1] = 1;
      mat.setarraysize(1, dims);
      getdatasafe<int>(mlmesh,"mat",mat);

      dims[0] = nel; dims[1] = nsd*nen;
      LM.setarraysize(2, dims);
      getdatasafe<int>(mlmesh,"LM",LM);

      dims[0] = nel; dims[1] = nen;
      IEN.setarraysize(2, dims);
      getdatasafe<int>(mlmesh,"IEN",IEN);

      dims[0] = np; dims[1] = nsd;
      ID.setarraysize(2, dims);
      getdatasafe<int>(mlmesh,"ID",ID);

      dims[0] = np; dims[1] = nsd;
      X.setarraysize(2, dims);
      getdatasafe<double>(mlmesh,"X",X);
   }
};

class physics {
 public:
    double lam, mu, rho0;

 public:
   explicit physics() {lam=0.0; mu=0.0; rho0=0.0;}
   explicit physics(mxArray const* mlphys) {
      getscalarsafe(mlphys,"lam", lam);
      getscalarsafe(mlphys,"mu", mu);
      getscalarsafe(mlphys,"rho0", rho0);
   }
   void setProperties(mxArray const* mlphys, int index){
      getscalarsafe(mlphys,"lam", lam, index);
      getscalarsafe(mlphys,"mu", mu, index);
      getscalarsafe(mlphys,"rho0", rho0, index);
   }
};

class integration {
 public:
    int nint;
    array<double> W;  
    array<double> Z;

 public:
   explicit integration(mesh msh) {

      nint = msh.nen;
      
      size_t dims[2];
      dims[0]=nint; dims[1]=msh.nsd+1;
      W.setarraysize(1, dims); 
      if (!W._alloc)
         W.allocmem(); 
      else
         W.zeroout();

      Z.setarraysize(2, dims);
      if (!Z._alloc)
         Z.allocmem(); 
      else
         Z.zeroout();

      //This ensures compatibility with node number in ShapeFunctions (in zfem_elems.cpp)
      Z(1,1)=1.0;
      Z(2,2)=1.0;
      Z(3,3)=1.0;
   }

   explicit integration(mxArray const* mlquad, mesh msh) {

      getscalarsafe(mlquad,"nint", nint);

      size_t dims[2];
      dims[0] = nint;
      W.setarraysize(1, dims);
      getdatasafe<double>(mlquad,"W",W);

      dims[0] = nint; dims[1] = msh.nsd+1;
      Z.setarraysize(2, dims);
      getdatasafe<double>(mlquad,"Z",Z);
   }
};

class DirichletBCs {
 public:
    int ndbc;
    array<bool> loc;  
    array<double> udbc;
    array<double> bdbc;

 public:
   explicit DirichletBCs(mxArray const* mldbc, mesh msh) {

      getscalarsafe(mldbc,"ndbc", ndbc);

      size_t dims[2];
      dims[0] = msh.np; dims[1] = msh.nsd;
      loc.setarraysize(2, dims);
      getdatasafe(mldbc,"loc",loc);

      dims[0]=ndbc; dims[1] = 0;
      udbc.setarraysize(1, dims);
      getdatasafe<double>(mldbc,"udbc",udbc);

      bdbc.setarraysize(1, dims);
      getdatasafe<double>(mldbc,"bdbc",bdbc);
   }
};

/*struct mesh
   {  int np, nel, nsd, nen, ndof;
      int** LM;
      size_t** IEN;
      double **X; };

struct physics
   {   double lam, mu, rho0;  };

struct integration
   {   int nint;
       double*  W;
       double** Z;  };

struct DirichletBCs
   {   int ndbc;
       bool** loc;
       double* udbc;
       double* bdbc;  };
*/

#endif
