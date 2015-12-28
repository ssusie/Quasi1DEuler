#include "zfem_mass.hpp"
#include "zfem_elems.hpp"
#include "matops.hpp"
#include <assert.h>
#include <stdio.h>

void MassMatrix(mesh* msh, physics* phys, integration* quad, array<int> M_irow_cr, array<int> M_jcol, array<double> M_val, int nnz)
{
    double J;
    //Zero out mass matrix
    for (int i=0; i < nnz; ++i)
       M_val(i) = 0.0;

    double* N   = new double [msh->nen];
    double** dN = new double* [msh->nen];
    double** Xe = new double* [msh->nsd];
    double** dXdZ = new double* [msh->nsd];
    double** me = new double* [msh->nen];

    for (int i=0; i < msh->nsd; ++i)
       { Xe[i] = new double [msh->nen];
         dXdZ[i] = new double [msh->nsd]; }

    for (int i=0; i < msh->nen; ++i)
       { dN[i] = new double [msh->nsd];
         me[i] = new double [msh->nen]; }

    for (int e=0; e < msh->nel; ++e){
        //Extract element positions (reference) and displacements,
        //taking DBCs into account
        for (int i=0; i < msh->nsd; ++i)
           for (int j=0; j < msh->nen; ++j)
              Xe[i][j]=(msh->X)( msh->IEN(e, j), i);

        //Zero out me
        for (int i=0; i < msh->nen; ++i)
           for (int j=0; j < msh->nen; ++j)
              me[i][j]=0.0;

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

           //Elemental contribution to mass
           for (int i=0; i < msh->nen; ++i)
              for (int j=0; j < msh->nen; ++j)
                  me[i][j]+=(1.0/6.0)*((quad->W)(k))*(phys[msh->mat(e)].rho0)*J*N[i]*N[j];
        }

        for (int i=0; i < (msh->nen); ++i) {
           for (int s=0; s < (msh->nsd); ++s) {
              int row = (msh->LM)(e,s+(msh->nsd)*i);
              if (row < 0)
                 continue;

              int index;
              int lowbnd = M_irow_cr(row);
              int uppbnd = M_irow_cr(row+1);

              for (int j=0; j < (msh->nen); ++j) {
                 for (int t=0; t < (msh->nsd); ++t) {
                    int col = (msh->LM)(e,t+(msh->nsd)*j);
                    if (col < 0 || (s!=t))
                       continue;

                    for (int p=lowbnd; p< uppbnd; ++p)
                    {   if ((M_jcol(p)-1) == col)
                          { index = p; break;}  }

                    M_val( index )+=me[i][j];
                 }
              }
           }
        }
    }
    for (int i=0; i < msh->nsd; ++i)
       { delete[] Xe[i];
         delete[] dXdZ[i]; }

    for (int i=0; i < msh->nen; ++i)
       { delete[] dN[i];
         delete[] me[i]; }
    delete[] N;
    delete[] dN;
    delete[] Xe;
    delete[] dXdZ;
    delete[] me;
}

void MassTimesDbcVec(mesh* msh, physics* phys, integration* quad, DirichletBCs* dbc, array<double> a, array<double> Ma) {

    double J;

    //Zero out a
    for (int i=0; i < msh->ndof; ++i)
       Ma(i) = 0.0;

    double* N   = new double [msh->nen];
    double** dN = new double* [msh->nen];
    double** Xe = new double* [msh->nsd];
    double** dXdZ = new double* [msh->nsd];
    double** me = new double* [msh->nen];

    for (int i=0; i < msh->nsd; ++i)
       { Xe[i] = new double [msh->nen];
         dXdZ[i] = new double [msh->nsd]; }

    for (int i=0; i < msh->nen; ++i)
       { dN[i] = new double [msh->nsd];
         me[i] = new double [msh->nen]; }

    for (int e=0; e < msh->nel; ++e){
        //Extract element positions (reference) and displacements,
        //taking DBCs into account
        for (int i=0; i < msh->nsd; ++i)
           for (int j=0; j < msh->nen; ++j)
              Xe[i][j]=(msh->X)( msh->IEN(e, j), i);

        //Zero out me
        for (int i=0; i < msh->nen; ++i)
           for (int j=0; j < msh->nen; ++j)
              me[i][j]=0.0;

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

           //Elemental contribution to mass
           for (int i=0; i < msh->nen; ++i)
              for (int j=0; j < msh->nen; ++j)
                  me[i][j]+=(1.0/6.0)*((quad->W)(k))*(phys[msh->mat(e)].rho0)*J*N[i]*N[j];
        }

        for (int i=0; i < msh->nen*msh->nsd; ++i) {
           int row=i/msh->nsd;
           int idof=i-(row+1)*msh->nsd;
           for (int j=0; j < msh->nen*msh->nsd; ++j) {
              int col=j/msh->nsd;
              int jdof=j-(col+1)*msh->nsd;
              if ( (msh->LM)(e,i)>-1 && (msh->LM)(e,j)<0 && ( idof == jdof ))
                 Ma( (msh->LM)(e,i) ) += me[row][col]*a( -(msh->LM)(e,j)-1 );
           }
        }
    }

    for (int i=0; i < msh->nsd; ++i)
       { delete[] Xe[i];
         delete[] dXdZ[i]; }

    for (int i=0; i < msh->nen; ++i)
       { delete[] dN[i];
         delete[] me[i]; }
    delete[] N;
    delete[] dN;
    delete[] Xe;
    delete[] dXdZ;
    delete[] me;
}





void MassMatrix(mesh* msh, physics* phys, integration* quad, array<double> M) {

    double J;

    //Zero out mass matrix
    for (int i=0; i < msh->np; ++i)
       for (int j=0; j < msh->np; ++j)
          M(i,j) = 0.0;

    double* N   = new double [msh->nen];
    double** dN = new double* [msh->nen];
    double** Xe = new double* [msh->nsd];
    double** dXdZ = new double* [msh->nsd];
    double** me = new double* [msh->nen];

    for (int i=0; i < msh->nsd; ++i)
       { Xe[i] = new double [msh->nen];
         dXdZ[i] = new double [msh->nsd]; }

    for (int i=0; i < msh->nen; ++i)
       { dN[i] = new double [msh->nsd];
         me[i] = new double [msh->nen]; }

    for (int e=0; e < msh->nel; ++e){
        //Extract element positions (reference) and displacements,
        //taking DBCs into account
        for (int i=0; i < msh->nsd; ++i)
           for (int j=0; j < msh->nen; ++j)
              Xe[i][j]=(msh->X)( msh->IEN(e, j), i);

        //Zero out me
        for (int i=0; i < msh->nen; ++i)
           for (int j=0; j < msh->nen; ++j)
              me[i][j]=0.0;

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

           //Elemental contribution to mass
           for (int i=0; i < msh->nen; ++i)
              for (int j=0; j < msh->nen; ++j)
                  me[i][j]+=(1.0/6.0)*((quad->W)(k))*(phys[msh->mat(e)].rho0)*J*N[i]*N[j];
        }

        for (int i=0; i < msh->nen; ++i)
           for (int j=0; j < msh->nen; ++j)
              M( (msh->IEN)(e,i) , (msh->IEN)(e,j) )+=me[i][j];
    }

    for (int i=0; i < msh->nsd; ++i)
       { delete[] Xe[i];
         delete[] dXdZ[i]; }

    for (int i=0; i < msh->nen; ++i)
       { delete[] dN[i];
         delete[] me[i]; }
    delete[] N;
    delete[] dN;
    delete[] Xe;
    delete[] dXdZ;
    delete[] me;
}

void MassTimesStateVec(mesh* msh, DirichletBCs* dbc, array<double> M, double* a, array<double> Ma) {
 //Mass matrix times state vector

    for (int i=0; i<(msh->ndof); ++i)
       Ma(i) = 0.0;

    for (int i=0; i<(msh->np); ++i) {
       for (int k=0; k<(msh->nsd); ++k){
          for (int j=0; j<(msh->np); ++j){
             if (!dbc->loc(i,k) && !dbc->loc(j,k))
               Ma( msh->ID(i,k) ) += M( i , j )*a[ msh->ID(j,k) ];
          }
       }
    } 
}

void MassTimesDbcVec(mesh* msh, DirichletBCs* dbc, array<double> M, array<double> Ma) {
 //Mass matrix times DBC vector

    for (int i=0; i<(msh->ndof); ++i)
       Ma(i) = 0.0;

    for (int i=0; i<(msh->np); ++i) {
       for (int k=0; k<(msh->nsd); ++k){
          for (int j=0; j<(msh->np); ++j){
             if (!dbc->loc(i,k) && dbc->loc(j,k))
               Ma( msh->ID(i,k) ) += M( i , j )*dbc->udbc( -msh->ID(j,k) - 1);
          }
       }
    } 
}

void MassPlusMatrix(mesh* msh, DirichletBCs* dbc, array<double> M, array<double> K) {
 //(full) Mass Matrix plus full matrix K (result in K)

    for (int i=0; i<(msh->np); ++i) {
       for (int k=0; k<(msh->nsd); ++k){
          for (int j=0; j<(msh->np); ++j){
             if (!dbc->loc(i,k) && !dbc->loc(j,k))
                K( msh->ID(i,k) , msh->ID(j,k) ) += M( i , j );
          }
       }
    }
}
