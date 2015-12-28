#include "zfem_forces.hpp"
#include "zfem_elems.hpp"
#include "matops.hpp"

#include <assert.h>
#include <stdio.h>

void  InternalForce(double* U, mesh* msh, physics* phys, DirichletBCs* dbc, integration* quad, array<double> f) {

    for (int i=0; i < msh->ndof; ++i)
       f(i)=0.0;

    //Make arrays for elemental positions (reference) and displacements
    double* N   = new double [msh->nen];
    double** dN = new double* [msh->nen];
    double** fe = new double* [msh->nsd];
    double** Xe = new double* [msh->nsd];
    double** dNdX = new double* [msh->nen];
    double** dXdZ = new double* [msh->nsd];
    double** dZdX = new double* [msh->nsd];
    double** ue = new double* [msh->nsd];
    double** P  = new double* [msh->nsd];
    double** H  = new double* [msh->nsd];
    double** S  = new double* [msh->nsd];

    for (int i=0; i < msh->nsd; ++i)
       { Xe[i] = new double [msh->nen];
         ue[i] = new double [msh->nen];
         fe[i] = new double [msh->nen]; 
         P[i]  = new double [msh->nsd];
         H[i]  = new double [msh->nsd];
         S[i]  = new double [msh->nsd];
         dXdZ[i] =  new double [msh->nsd];
         dZdX[i] =  new double [msh->nsd]; }
         
    for (int i=0; i < msh->nen; ++i) {
       dN[i] = new double [msh->nsd];
       dNdX[i] = new double [msh->nsd]; }

    double J;

    for (int e=0; e < msh->nel; ++e){
        //Extract element positions (reference) and displacements,
        //taking DBCs into account
        for (int i=0; i < msh->nsd; ++i) {
           for (int j=0; j < msh->nen; ++j) {
              Xe[i][j]=(msh->X)( msh->IEN(e, j), i);
              if ((msh->LM)(e, i+msh->nsd*j) > -1)
                 ue[i][j]=U[(msh->LM)( e, i+msh->nsd*j )];
              else
                 ue[i][j]=(dbc->udbc)(-(msh->LM)( e, i+msh->nsd*j )-1);
           }
        }

        //Zero out fe and ke
        for (int i=0; i < msh->nsd; ++i)
           for (int j=0; j < msh->nen; ++j)
              fe[i][j]=0.0;

        for (int k=0; k< quad->nint; ++k){
           //Compute shape functions
           ShapeFunctions(quad, k, N, dN);

           //Compute jacobian of map      
           for (int i=0; i < msh->nsd; ++i) {
              for (int j=0; j < msh->nsd; ++j) {
                 dXdZ[i][j] = 0.0;
                 for (int p=0; p < msh->nen; ++p)
                    dXdZ[i][j]+=Xe[i][p]*dN[p][j];
              }
           }
           
           //Compute dZdX = inv(dXdZ)
           inverse3x3(dXdZ,dZdX);
           //Compute J = det(dXdZ)
           J=determinant3x3(dXdZ);

           //Compute dNdX
           for (int i=0; i < msh->nen; ++i) {
              for (int j=0; j < msh->nsd; ++j) {
                 dNdX[i][j]=0.0;
                 for (int p=0; p < msh->nsd; ++p)
                    dNdX[i][j]+=dN[i][p]*dZdX[p][j];
              }
           }

          //Compute H = ue*dNdX
           for (int i=0; i < msh->nsd; ++i) {
              for (int j=0; j < msh->nsd; ++j) {
                 H[i][j]=0.0;
                 for (int p=0; p < msh->nen; ++p)
                    H[i][j]+=ue[i][p]*dNdX[p][j];
              }
           }
           //Compute first Piola Kirchoff stress from constitutive law
           ConstitutiveLaw( &(phys[msh->mat(e)]) , false, H , P , S);

           //Elemental contribution to force
           for (int i=0; i < msh->nsd; ++i)
              for (int j=0; j < msh->nen; ++j)
                 for (int p=0; p < msh->nsd; ++p)
                    fe[i][j]+=(1.0/6.0)*((quad->W)(k))*J*dNdX[j][p]*P[p][i];
        }

        //Fill appropriate locations of force vector
        for (int i=0; i < (msh->nsd)*(msh->nen); ++i) { 
           //i = row + nsd*col (row < 3)
           //col = floor(i/nsd)
           //row = i = 3*col
 
           int col = i/(msh->nsd);
           int row = i - (msh->nsd)*col;

           if ((msh->LM)(e, i) > -1)
              f((msh->LM)(e,i))+=fe[row][col];
        }
    }

    for (int i=0; i < msh->nsd; ++i)
       { delete[] Xe[i];
         delete[] ue[i]; 
         delete[] fe[i];
         delete[] P[i];
         delete[] S[i];
         delete[] H[i];
         delete[] dXdZ[i];
         delete[] dZdX[i]; }
    for (int i=0; i < msh->nen; ++i)
       { delete[] dN[i];
         delete[] dNdX[i]; }
 
    delete[] Xe;
    delete[] ue;
    delete[] fe;
    delete[] P;
    delete[] H;
    delete[] S;
    delete[] dXdZ;
    delete[] dZdX;
    delete[] N;
    delete[] dN;
    delete[] dNdX;
}

void  InternalForce(double* U, mesh* msh, physics* phys, DirichletBCs* dbc, integration* quad, array<double> f, array<double> df) {

    for (int i=0; i < msh->ndof; ++i)
       f(i)=0.0;
    for (int i=0; i < msh->ndof; ++i)
       for (int j=0; j < msh->ndof; ++j)
          df(i,j)=0.0;

    //Make arrays for elemental positions (reference) and displacements
    double* N   = new double [msh->nen];
    double** fe = new double* [msh->nsd];
    double** ke = new double* [msh->nsd*msh->nen];
    double** dN = new double* [msh->nen];
    double** Xe = new double* [msh->nsd];
    double** dNdX = new double* [msh->nen];
    double** dXdZ = new double* [msh->nsd];
    double** dZdX = new double* [msh->nsd];
    double** ue = new double* [msh->nsd];
    double** P  = new double* [msh->nsd];
    double** H  = new double* [msh->nsd];
    double** S  = new double* [msh->nsd];
    double** dNdNt  = new double* [msh->nen];
    double** FFt  = new double* [msh->nsd];
    double** dNSdNt= new double* [msh->nen];
    double** FdNt  = new double* [msh->nsd];
    for (int i=0; i < (msh->nsd)*(msh->nen); ++i)
       ke[i] = new double [(msh->nen)*(msh->nsd)];

    for (int i=0; i < msh->nsd; ++i)
       { Xe[i] = new double [msh->nen];
         ue[i] = new double [msh->nen];
         fe[i] = new double [msh->nen]; 
         P[i]  = new double [msh->nsd];
         H[i]  = new double [msh->nsd];
         S[i]  = new double [msh->nsd];
         FdNt[i]  = new double [msh->nen];
         FFt[i]  = new double [msh->nsd];
         dXdZ[i] =  new double [msh->nsd];
         dZdX[i] =  new double [msh->nsd]; }
         
    for (int i=0; i < msh->nen; ++i) {
       dN[i] = new double [msh->nsd];
       dNdX[i] = new double [msh->nsd];
       dNdNt[i] = new double [msh->nen];
       dNSdNt[i] = new double [msh->nen]; }

    double J;

    for (int e=0; e < msh->nel; ++e){
        //Extract element positions (reference) and displacements,
        //taking DBCs into account
        for (int i=0; i < msh->nsd; ++i) {
           for (int j=0; j < msh->nen; ++j) {
              Xe[i][j]=(msh->X)( msh->IEN(e, j), i);
              if ((msh->LM)(e, i+msh->nsd*j) > -1)
                 ue[i][j]=U[(msh->LM)( e, i+msh->nsd*j )];
              else
                 ue[i][j]=(dbc->udbc)(-(msh->LM)( e, i+msh->nsd*j )-1);
           }
        }

        //Zero out fe and ke
        for (int i=0; i < msh->nsd; ++i)
           for (int j=0; j < msh->nen; ++j)
              fe[i][j]=0.0;

        for (int i=0; i < msh->nsd*msh->nen; ++i)
           for (int j=0; j < msh->nsd*msh->nen; ++j)
              ke[i][j]=0.0;

        for (int k=0; k< quad->nint; ++k){
           //Compute shape functions
           ShapeFunctions(quad, k, N, dN);

           //Compute jacobian of map      
           for (int i=0; i < msh->nsd; ++i) {
              for (int j=0; j < msh->nsd; ++j) {
                 dXdZ[i][j] = 0.0;
                 for (int p=0; p < msh->nen; ++p)
                    dXdZ[i][j]+=Xe[i][p]*dN[p][j];
              }
           }
           
           //Compute dZdX = inv(dXdZ)
           inverse3x3(dXdZ,dZdX);
           //Compute J = det(dXdZ)
           J=determinant3x3(dXdZ);

           //Compute dNdX
           for (int i=0; i < msh->nen; ++i) {
              for (int j=0; j < msh->nsd; ++j) {
                 dNdX[i][j]=0.0;
                 for (int p=0; p < msh->nsd; ++p)
                    dNdX[i][j]+=dN[i][p]*dZdX[p][j];
              }
           }

          //Compute H = ue*dNdX
           for (int i=0; i < msh->nsd; ++i) {
              for (int j=0; j < msh->nsd; ++j) {
                 H[i][j]=0.0;
                 for (int p=0; p < msh->nen; ++p)
                    H[i][j]+=ue[i][p]*dNdX[p][j];
              }
           }
           //Compute first Piola Kirchoff stress from constitutive law
           //ConstitutiveLaw( phys , H , P , S);
           ConstitutiveLaw( &(phys[msh->mat(e)]) , false, H , P , S);

           //Compute F*dNdX
           for (int i=0; i < msh->nsd; ++i)
              for (int j=0; j < msh->nen; ++j)
                 FdNt[i][j] = 0.0;

           for (int i=0; i < msh->nsd; ++i) {
              for (int j=0; j < msh->nen; ++j) {
                 for (int p=0; p < msh->nsd; ++p){
                    if (i == p)
                       FdNt[i][j]+=(H[i][p]+1.0)*dNdX[j][p]; 
                    else
                       FdNt[i][j]+=H[i][p]*dNdX[j][p];
                 }
              }
           }

           //Compute dNdX*dNdX^T
           for (int i =0; i < msh->nen; ++i)
              for (int j =0; j < msh->nen; ++j)
                 dNdNt[i][j] = 0.0;

           for (int i =0; i < msh->nen; ++i)
              for (int j =0; j < msh->nen; ++j)
                 for (int p=0; p < msh->nsd; ++p)
                    dNdNt[i][j] += dNdX[i][p]*dNdX[j][p];

           //Compute F*F^T
           for (int i=0; i < msh->nsd; ++i) {
              for (int j=0; j < msh->nsd; ++j) {
                 FFt[i][j]=H[i][j]+H[j][i];
                 if (i == j)
                    FFt[i][j]+=1.0;
              }
           }

           for (int i=0; i < msh->nsd; ++i)
              for (int j=0; j < msh->nsd; ++j) 
                 for (int p=0; p < msh->nsd; ++p)
                   FFt[i][j]+=H[i][p]*H[j][p];

           //Compute dNdX*S*dNdX^T
           for (int i =0; i < msh->nen; ++i)
              for (int j =0; j < msh->nen; ++j)
                 dNSdNt[i][j] = 0.0;

           for (int i =0; i < msh->nen; ++i)
              for (int j =0; j < msh->nen; ++j)
                 for (int k =0; k < msh->nsd; ++k)
                    for (int l =0; l < msh->nsd; ++l)
                       dNSdNt[i][j]+=dNdX[i][k]*S[k][l]*dNdX[j][l];
  
           //Elemental contribution to force
           for (int i=0; i < msh->nsd; ++i)
              for (int j=0; j < msh->nen; ++j)
                 for (int p=0; p < msh->nsd; ++p)
                    fe[i][j]+=(1.0/6.0)*((quad->W)(k))*J*dNdX[j][p]*P[p][i];

           //Elemental contribution to stiffness
           for (int ir=0; ir < msh->nsd; ++ir) {
              for (int jr=0; jr < msh->nen; ++jr) {
                 for (int ic=0; ic < msh->nsd; ++ic) {
                    for (int jc=0; jc < msh->nen; ++jc) {
                       ke[ir+jr*(msh->nsd)][ic+jc*(msh->nsd)] += (1.0/6.0)*((quad->W)(k))*J*( (phys[msh->mat(e)].lam)*FdNt[ir][jr]*FdNt[ic][jc] + (phys[msh->mat(e)].mu)*FFt[ir][ic]*dNdNt[jc][jr] + (phys[msh->mat(e)].mu)*FdNt[ir][jc]*FdNt[ic][jr]);
                       if (ir == ic)
                         ke[ir+jr*(msh->nsd)][ic+jc*(msh->nsd)] += (1.0/6.0)*((quad->W)(k))*J*dNSdNt[jr][jc]; 
                    }
                 }
              }
           }
        }

        /*for (int i=0; i < (msh->nsd)*(msh->nen); ++i){
           for (int j=0; j < (msh->nsd)*(msh->nen); ++j){
              std::cout << ke[i][j] << "  ";
           }
           std::cout << std::endl;
        }*/

        //Fill appropriate locations of force vector
        for (int i=0; i < (msh->nsd)*(msh->nen); ++i) { 
           //i = row + nsd*col (row < 3)
           //col = floor(i/nsd)
           //row = i = 3*col
 
           int col = i/(msh->nsd);
           int row = i - (msh->nsd)*col;
           //std::cout << (msh->LM)[i][e] << std::endl;

           if ((msh->LM)(e, i) > -1)
              f((msh->LM)(e,i))+=fe[row][col];
        }

        
        for (int i=0; i < (msh->nsd)*(msh->nen); ++i) { 
           for (int j=0; j < (msh->nsd)*(msh->nen); ++j) {
              if ((msh->LM)(e,i) > -1 && (msh->LM)(e,j) > -1)
                 df( (msh->LM)(e,i) , (msh->LM)(e,j) )+=ke[i][j];
           }
        }
    }

    for (int i=0; i < msh->nsd; ++i)
       { delete[] Xe[i];
         delete[] ue[i]; 
         delete[] fe[i];
         delete[] P[i];
         delete[] S[i];
         delete[] H[i];
         delete[] FdNt[i];
         delete[] FFt[i];
         delete[] dXdZ[i];
         delete[] dZdX[i]; }
    for (int i=0; i < msh->nen; ++i)
       { delete dN[i];
         delete dNdX[i];
         delete dNdNt[i];
         delete dNSdNt[i]; } 
    for (int i=0; i < msh->nen*msh->nsd; ++i)
       delete[] ke[i];
 
    delete[] Xe;
    delete[] ue;
    delete[] fe;
    delete[] P;
    delete[] H;
    delete[] S;
    delete[] FdNt;
    delete[] FFt;
    delete[] dNdNt;
    delete[] dNSdNt;
    delete[] dXdZ;
    delete[] dZdX;
    delete[] N;
    delete[] dN;
    delete[] dNdX;
    delete[] ke;
}

void InternalForce(double* U, mesh* msh, physics* phys, DirichletBCs* dbc, integration* quad, array<double> f, array<int> df_irow_cr, array<int> df_jcol, array<double> df_val, int nnz)
{ 
    for (int i=0; i < msh->ndof; ++i)
       f(i)=0.0;
    for (int i=0; i < nnz; ++i)
       df_val(i)=0.0;

    //Make arrays for elemental positions (reference) and displacements
    double* N   = new double [msh->nen];
    double** fe = new double* [msh->nsd];
    double** ke = new double* [msh->nsd*msh->nen];
    double** dN = new double* [msh->nen];
    double** Xe = new double* [msh->nsd];
    double** dNdX = new double* [msh->nen];
    double** dXdZ = new double* [msh->nsd];
    double** dZdX = new double* [msh->nsd];
    double** ue = new double* [msh->nsd];
    double** P  = new double* [msh->nsd];
    double** H  = new double* [msh->nsd];
    double** S  = new double* [msh->nsd];
    double** dNdNt  = new double* [msh->nen];
    double** FFt  = new double* [msh->nsd];
    double** dNSdNt= new double* [msh->nen];
    double** FdNt  = new double* [msh->nsd];
    for (int i=0; i < (msh->nsd)*(msh->nen); ++i)
       ke[i] = new double [(msh->nen)*(msh->nsd)];

    for (int i=0; i < msh->nsd; ++i)
       { Xe[i] = new double [msh->nen];
         ue[i] = new double [msh->nen];
         fe[i] = new double [msh->nen];
         P[i]  = new double [msh->nsd];
         H[i]  = new double [msh->nsd];
         S[i]  = new double [msh->nsd];
         FdNt[i]  = new double [msh->nen];
         FFt[i]  = new double [msh->nsd];
         dXdZ[i] =  new double [msh->nsd];
         dZdX[i] =  new double [msh->nsd]; }

    for (int i=0; i < msh->nen; ++i) {
       dN[i] = new double [msh->nsd];
       dNdX[i] = new double [msh->nsd];
       dNdNt[i] = new double [msh->nen];
       dNSdNt[i] = new double [msh->nen]; }

    double J;

    for (int e=0; e < msh->nel; ++e){
        //Extract element positions (reference) and displacements,
        //taking DBCs into account
        for (int i=0; i < msh->nsd; ++i) {
           for (int j=0; j < msh->nen; ++j) {
              Xe[i][j]=(msh->X)( msh->IEN(e, j), i);
              if ((msh->LM)(e, i+msh->nsd*j) > -1)
                 ue[i][j]=U[(msh->LM)( e, i+msh->nsd*j )];
              else
                 ue[i][j]=(dbc->udbc)(-(msh->LM)( e, i+msh->nsd*j )-1);
           }
        }

        //Zero out fe and ke
        for (int i=0; i < msh->nsd; ++i)
           for (int j=0; j < msh->nen; ++j)
              fe[i][j]=0.0;

        for (int i=0; i < msh->nsd*msh->nen; ++i)
           for (int j=0; j < msh->nsd*msh->nen; ++j)
              ke[i][j]=0.0;

        for (int k=0; k< quad->nint; ++k){
           //Compute shape functions
           ShapeFunctions(quad, k, N, dN);

           //Compute jacobian of map      
           for (int i=0; i < msh->nsd; ++i) {
              for (int j=0; j < msh->nsd; ++j) {
                 dXdZ[i][j] = 0.0;
                 for (int p=0; p < msh->nen; ++p)
                    dXdZ[i][j]+=Xe[i][p]*dN[p][j];
              }
           }

           //Compute dZdX = inv(dXdZ)
           inverse3x3(dXdZ,dZdX);
           //Compute J = det(dXdZ)
           J=determinant3x3(dXdZ);

           //Compute dNdX
           for (int i=0; i < msh->nen; ++i) {
              for (int j=0; j < msh->nsd; ++j) {
                 dNdX[i][j]=0.0;
                 for (int p=0; p < msh->nsd; ++p)
                    dNdX[i][j]+=dN[i][p]*dZdX[p][j];
              }
           }

          //Compute H = ue*dNdX
           for (int i=0; i < msh->nsd; ++i) {
              for (int j=0; j < msh->nsd; ++j) {
                 H[i][j]=0.0;
                 for (int p=0; p < msh->nen; ++p)
                    H[i][j]+=ue[i][p]*dNdX[p][j];
              }
           }
           //Compute first Piola Kirchoff stress from constitutive law
           //ConstitutiveLaw( phys , H , P , S);
           ConstitutiveLaw( &(phys[msh->mat(e)]) , false, H , P , S);

           //Compute F*dNdX
           for (int i=0; i < msh->nsd; ++i)
              for (int j=0; j < msh->nen; ++j)
                 FdNt[i][j] = 0.0;

           for (int i=0; i < msh->nsd; ++i) {
              for (int j=0; j < msh->nen; ++j) {
                 for (int p=0; p < msh->nsd; ++p){
                    if (i == p)
                       FdNt[i][j]+=(H[i][p]+1.0)*dNdX[j][p];
                    else
                       FdNt[i][j]+=H[i][p]*dNdX[j][p];
                 }
              }
           }

           //Compute dNdX*dNdX^T
           for (int i =0; i < msh->nen; ++i)
              for (int j =0; j < msh->nen; ++j)
                 dNdNt[i][j] = 0.0;

           for (int i =0; i < msh->nen; ++i)
              for (int j =0; j < msh->nen; ++j)
                 for (int p=0; p < msh->nsd; ++p)
                    dNdNt[i][j] += dNdX[i][p]*dNdX[j][p];

           //Compute F*F^T
           for (int i=0; i < msh->nsd; ++i) {
              for (int j=0; j < msh->nsd; ++j) {
                 FFt[i][j]=H[i][j]+H[j][i];
                 if (i == j)
                    FFt[i][j]+=1.0;
              }
           }

           for (int i=0; i < msh->nsd; ++i)
              for (int j=0; j < msh->nsd; ++j)
                 for (int p=0; p < msh->nsd; ++p)
                   FFt[i][j]+=H[i][p]*H[j][p];

           //Compute dNdX*S*dNdX^T
           for (int i =0; i < msh->nen; ++i)
              for (int j =0; j < msh->nen; ++j)
                 dNSdNt[i][j] = 0.0;

           for (int i =0; i < msh->nen; ++i)
              for (int j =0; j < msh->nen; ++j)
                 for (int k =0; k < msh->nsd; ++k)
                   for (int l =0; l < msh->nsd; ++l)
                       dNSdNt[i][j]+=dNdX[i][k]*S[k][l]*dNdX[j][l];

           //Elemental contribution to force
           for (int i=0; i < msh->nsd; ++i)
              for (int j=0; j < msh->nen; ++j)
                 for (int p=0; p < msh->nsd; ++p)
                    fe[i][j]+=(1.0/6.0)*((quad->W)(k))*J*dNdX[j][p]*P[p][i];

           //Elemental contribution to stiffness
           for (int ir=0; ir < msh->nsd; ++ir) {
              for (int jr=0; jr < msh->nen; ++jr) {
                 for (int ic=0; ic < msh->nsd; ++ic) {
                    for (int jc=0; jc < msh->nen; ++jc) {
                       ke[ir+jr*(msh->nsd)][ic+jc*(msh->nsd)] += (1.0/6.0)*((quad->W)(k))*J*( (phys[msh->mat(e)].lam)*FdNt[ir][jr]*FdNt[ic][jc] + (phys[msh->mat(e)].mu)*FFt[ir][ic]*dNdNt[jc][jr] + (phys[msh->mat(e)].mu)*FdNt[ir][jc]*FdNt[ic][jr]);
                       if (ir == ic)
                         ke[ir+jr*(msh->nsd)][ic+jc*(msh->nsd)] += (1.0/6.0)*((quad->W)(k))*J*dNSdNt[jr][jc];
                    }
                 }
              }
           }
        }

        //Fill appropriate locations of force vector
        for (int i=0; i < (msh->nsd)*(msh->nen); ++i) {
           //i = row + nsd*col (row < 3)
           //col = floor(i/nsd)
           //row = i = 3*col

           int col = i/(msh->nsd);
           int row = i - (msh->nsd)*col;
           //std::cout << (msh->LM)[i][e] << std::endl;

           if ((msh->LM)(e, i) > -1)
              f((msh->LM)(e,i))+=fe[row][col];
        }

        for (int i=0; i < (msh->nsd)*(msh->nen); ++i) {
           int row = (msh->LM)(e,i);
           if (row < 0)
              continue;

           int index;
           int lowbnd = df_irow_cr(row);
           int uppbnd = df_irow_cr(row+1);             

           for (int j=0; j < (msh->nsd)*(msh->nen); ++j) {
              int col = (msh->LM)(e,j);
              if (col < 0)
                 continue;

              for (int p=lowbnd; p< uppbnd; ++p)
              {   if ((df_jcol(p)-1) == col)
                    { index = p; break;}  }
                 
              df_val( index )+=ke[i][j];
           }
        }
    }

    for (int i=0; i < msh->nsd; ++i)
       { delete[] Xe[i];
         delete[] ue[i];
         delete[] fe[i];
         delete[] P[i];
         delete[] S[i];
         delete[] H[i];
         delete[] FdNt[i];
         delete[] FFt[i];
         delete[] dXdZ[i];
         delete[] dZdX[i]; }
    for (int i=0; i < msh->nen; ++i)
       { delete dN[i];
         delete dNdX[i];
         delete dNdNt[i];
         delete dNSdNt[i]; }
    for (int i=0; i < msh->nen*msh->nsd; ++i)
       delete[] ke[i];

    delete[] Xe;
    delete[] ue;
    delete[] fe;
    delete[] P;
    delete[] H;
    delete[] S;
    delete[] FdNt;
    delete[] FFt;
    delete[] dNdNt;
    delete[] dNSdNt;
    delete[] dXdZ;
    delete[] dZdX;
    delete[] N;
    delete[] dN;
    delete[] dNdX;
    delete[] ke;
}

void BodyForce(double* b, mesh* msh, physics* phys, DirichletBCs* dbc, integration* quad, array<double> f) {

    for (int i=0; i < msh->ndof; ++i)
       f(i)=0.0;

    //Make arrays for elemental positions (reference) and displacements
    double* N   = new double [msh->nen];
    double* b_interp = new double [msh->nsd];
    double** fe = new double* [msh->nsd];
    double** dN = new double* [msh->nen];
    double** Xe = new double* [msh->nsd];
    double** dXdZ = new double* [msh->nsd];
    double** be = new double* [msh->nsd];
    for (int i=0; i < msh->nsd; ++i)
       { Xe[i] = new double [msh->nen];
         be[i] = new double [msh->nen];
         fe[i] = new double [msh->nen];
         dXdZ[i] =  new double [msh->nsd]; }

    for (int i=0; i < msh->nen; ++i) 
       dN[i] = new double [msh->nsd];

    double J;

    for (int e=0; e < msh->nel; ++e){
        //Extract element positions (reference) and displacements,
        //taking DBCs into account
        for (int i=0; i < msh->nsd; ++i) {
           for (int j=0; j < msh->nen; ++j) {
              Xe[i][j]=(msh->X)(msh->IEN(e,j), i);
              if ((msh->LM)(e,i+msh->nsd*j) > -1)
                 be[i][j]=b[(msh->LM)(e,i+msh->nsd*j)];
              else
                 be[i][j]=(dbc->bdbc)(-(msh->LM)(e,i+msh->nsd*j)-1);
           }
        }

        //Zero out fe
        for (int i=0; i < msh->nsd; ++i)
           for (int j=0; j < msh->nen; ++j)
              fe[i][j]=0.0;

        for (int k=0; k< quad->nint; ++k){
           //Compute shape functions
           ShapeFunctions(quad, k, N, dN);

           //Compute jacobian of map      
           for (int i=0; i < msh->nsd; ++i) {
              for (int j=0; j < msh->nsd; ++j) {
                 dXdZ[i][j] = 0.0;
                 for (int p=0; p < msh->nen; ++p)
                    dXdZ[i][j]+=Xe[i][p]*dN[p][j];
              }
           }

           //Compute J = det(dXdZ)
           J=determinant3x3(dXdZ);

          //Compute b_interp = be*N
           for (int i=0; i < msh->nsd; ++i) {
              b_interp[i]=0.0;
              for (int j=0; j < msh->nen; ++j) {
                 b_interp[i]+=be[i][j]*N[j];
              }
           }

           //Elemental c ntribution to force
           for (int i=0; i < msh->nsd; ++i)
              for (int j=0; j < msh->nen; ++j)
                 fe[i][j]+=(1.0/6.0)*((quad->W)(k))*J*b_interp[i]*N[j]*(phys[msh->mat(e)].rho0);
        }

        //Fill appropriate locations of force vector
        for (int i=0; i < (msh->nsd)*(msh->nen); ++i) {
           //i = row + nsd*col (row < 3)
           //col = floor(i/nsd)
           //row = i = 3*col

           int col = i/(msh->nsd);
           int row = i - (msh->nsd)*col;

           if ((msh->LM)(e,i) > -1)
              f((msh->LM)(e,i))+=fe[row][col];
        }
    }

    for (int i=0; i < msh->nsd; ++i)
       { delete[] Xe[i];
         delete[] be[i];
         delete[] fe[i];
         delete[] dXdZ[i];}
    for (int i=0; i < msh->nen; ++i)
       delete dN[i];

    delete[] Xe;
    delete[] be;
    delete[] fe;
    delete[] dXdZ;
    delete[] N;
    delete[] b_interp;
    delete[] dN;
}

void  NodalStress(double* U, mesh* msh, physics* phys, DirichletBCs* dbc, array<double> sigma) {

    for (int i=0; i < msh->np; ++i)
       sigma(i)=0.0;

    
    integration quad(*msh);

    //Make arrays for elemental positions (reference) and displacements
    double* N   = new double [msh->nen];
    double** dN = new double* [msh->nen];
    double* se  = new double [msh->nen];
    double* prinstr  = new double [msh->nsd];
    double** Xe = new double* [msh->nsd];
    double** dNdX = new double* [msh->nen];
    double** dXdZ = new double* [msh->nsd];
    double** dZdX = new double* [msh->nsd];
    double** ue = new double* [msh->nsd];
    double** P  = new double* [msh->nsd];
    double** H  = new double* [msh->nsd];
    double** S  = new double* [msh->nsd];

    for (int i=0; i < msh->nsd; ++i)
       { Xe[i] = new double [msh->nen];
         ue[i] = new double [msh->nen];
         P[i]  = new double [msh->nsd];
         H[i]  = new double [msh->nsd];
         S[i]  = new double [msh->nsd];
         dXdZ[i] =  new double [msh->nsd];
         dZdX[i] =  new double [msh->nsd]; }
         
    for (int i=0; i < msh->nen; ++i) {
       dN[i] = new double [msh->nsd];
       dNdX[i] = new double [msh->nsd]; }

    for (int e=0; e < msh->nel; ++e){
        //Extract element positions (reference) and displacements,
        //taking DBCs into account
        for (int i=0; i < msh->nsd; ++i) {
           for (int j=0; j < msh->nen; ++j) {
              Xe[i][j]=(msh->X)( msh->IEN(e, j), i);
              if ((msh->LM)(e, i+msh->nsd*j) > -1)
                 ue[i][j]=U[(msh->LM)( e, i+msh->nsd*j )];
              else
                 ue[i][j]=(dbc->udbc)(-(msh->LM)( e, i+msh->nsd*j )-1);
           }
        }

        //Zero out se
        for (int i=0; i < msh->nen; ++i)
           se[i]=0.0;

        for (int k=0; k< msh->nen; ++k){
           ShapeFunctions(&quad, k, N, dN);

           //Compute jacobian of map      
           for (int i=0; i < msh->nsd; ++i) {
              for (int j=0; j < msh->nsd; ++j) {
                 dXdZ[i][j] = 0.0;
                 for (int p=0; p < msh->nen; ++p)
                    dXdZ[i][j]+=Xe[i][p]*dN[p][j];
              }
           }
           
           //Compute dZdX = inv(dXdZ)
           inverse3x3(dXdZ,dZdX);

           //Compute dNdX
           for (int i=0; i < msh->nen; ++i) {
              for (int j=0; j < msh->nsd; ++j) {
                 dNdX[i][j]=0.0;
                 for (int p=0; p < msh->nsd; ++p)
                    dNdX[i][j]+=dN[i][p]*dZdX[p][j];
              }
           }

          //Compute H = ue*dNdX
           for (int i=0; i < msh->nsd; ++i) {
              for (int j=0; j < msh->nsd; ++j) {
                 H[i][j]=0.0;
                 for (int p=0; p < msh->nen; ++p)
                    H[i][j]+=ue[i][p]*dNdX[p][j];
              }
           }
           //Compute first Piola Kirchoff stress from constitutive law
           //ConstitutiveLaw( phys , H , P , S);
           ConstitutiveLaw( &(phys[msh->mat(e)]) , false, H , P , S);

           //Get principle stresses
           eig3x3sym(S, prinstr);

           //Assume we want 1st principle stress for now
           se[k] = prinstr[0];
        }

        //Fill appropriate locations of force vector
        for (int i=0; i < msh->nen; ++i)
              sigma((msh->IEN)(e,i))+=se[i];
    }

    for (int i=0; i < msh->nsd; ++i)
       { delete[] Xe[i];
         delete[] ue[i]; 
         delete[] P[i];
         delete[] S[i];
         delete[] H[i];
         delete[] dXdZ[i];
         delete[] dZdX[i]; }
    for (int i=0; i < msh->nen; ++i)
       { delete[] dN[i];
         delete[] dNdX[i]; }
 
    delete[] Xe;
    delete[] ue;
    delete[] se;
    delete[] P;
    delete[] H;
    delete[] S;
    delete[] prinstr;
    delete[] dXdZ;
    delete[] dZdX;
    delete[] N;
    delete[] dN;
    delete[] dNdX;
}
