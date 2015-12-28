#include "zfem_elems.hpp"
#include <assert.h>
#include <stdio.h>

void ShapeFunctions(integration * quad, int num, double *N, double **dN) {

    //Linear shape functions for tetrahedra
    N[0] = 1-(quad->Z)(num,0)-(quad->Z)(num,1)-(quad->Z)(num,2); //Equivalent to: N[0]=(quad->Z)[num][3]
    N[1] = (quad->Z)(num,0);
    N[2] = (quad->Z)(num,1);
    N[3] = (quad->Z)(num,2);

    dN[0][0]=-1; dN[0][1]=-1; dN[0][2]=-1;
    dN[1][0]=1;  dN[1][1]=0;  dN[1][2]=0;
    dN[2][0]=0;  dN[2][1]=1;  dN[2][2]=0;
    dN[3][0]=0;  dN[3][1]=0;  dN[3][2]=1;
}

void ConstitutiveLaw(physics* phys, bool lin, double **H, double **P, double **S) {

    //assert(msh->nsd == 3);

    double E[3][3];
    //double S[3][3];
    double trE=0.0;

    //Green-Lagrange Stress Tensor
    for (int i=0; i<3; ++i) {
       for (int j=0; j<3; ++j) {
         E[i][j]=0.5*(H[i][j]+H[j][i]);
         for (int k=0; k<3; ++k)
           E[i][j]+=0.5*H[k][i]*H[k][j];
       }
       trE+=E[i][i];
    }

    //Second Piola Kirchoff Stress 
    for (int i=0; i<3; ++i) {
       for (int j=0; j<3; ++j) {
          S[i][j] = 2*(phys->mu)*E[i][j];
          if (i == j)
             S[i][i] += (phys->lam)*trE;
       }
    }

    //First Piola Kirchoff Stress Tensor
    if (lin) {
       for (int i=0; i<3; ++i)
          for (int j=0; j<3; ++j)
            P[i][j] = S[i][j];
    } else { 
       for (int i=0; i<3; ++i) {
          for (int j=0; j<3; ++j) {
             P[i][j] = S[i][j];
                for (int k=0; k<3; ++k)
                   P[i][j]+= S[i][k]*H[j][k];
          }
       }
    }
}
